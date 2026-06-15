library(tidyverse)
library(dplyr)
library(readxl)
library(saccades)
library(arrow)
library(data.table)
library(this.path)
library(runner)

### Read in DMS data

# Load in data
setwd(this.dir())
non_aug <- read_feather("./data/non_aug_data.feather")

# Replace repeated button-press observations with zeros
non_aug <- non_aug %>%
  group_by(rleid(Task_Available)) %>%
  mutate(
    # Correct button press is first frame pressed within task window
    Correct_Left_Button_Press = case_when(
      Left_Button_Press == 1 &
        !duplicated(Left_Button_Press) &
        Task_Available %in% c("Left", "Right") ~ 1,
      TRUE ~ 0
    ),
    # Left button press is any press that is not holding an ongoing press
    Left_Button_Press = case_when(
      (Left_Button_Press == lag(Left_Button_Press)) &
        (Left_Button_Press == 1) |
        (Left_Button_Press == 2) ~ 0,
      TRUE ~ Left_Button_Press
    )
  ) %>%
  ungroup()

# Replace coded 99 values from python script with NA values
non_aug[non_aug == 99] <- NA

# Use 1-base indexing for X and Frame_Num (match row number)
non_aug <- non_aug %>%
  mutate(X = 1:nrow(non_aug), .before = 1) %>%
  group_by(DaqName) %>%
  mutate(Frame_Num = row_number(), .after = X) %>%
  ungroup()

# Add average pupil size of both eyes
non_aug$Avg_Pupil_Diameter <- (non_aug$RPupil_Diameter +
  non_aug$LPupil_Diameter) / 2

# Rolling SD of vehicle lateral deviation
non_aug <- non_aug %>%
  group_by(DaqName) %>%
  mutate(Vehicle_Lat_Dev_rolling_sd3 = runner(
    Vehicle_Lat_Dev, sd,
    k = 180, na_pad = FALSE, na.rm = TRUE
  )) %>%
  ungroup()

# Get `ST_*` columns

output <- read_excel("./data/output_IIHS_withDMS.xls") %>%
  rename("ST_Right_N" = "ST_Rights_N") %>%
  select(matches("ST_(Right|Left).*"), DaqName, Subject, Drive) %>%
  relocate(DaqName, .after = Drive) %>%
  drop_na()

output$DaqName <- output$DaqName %>% str_sub(0, -5)

# Add counted hits from data set
output <- output %>% left_join(
  non_aug %>%
    filter(Task_Available %in% c("Left", "Right")) %>%
    group_by(DaqName, Task_Available) %>%
    summarize(
      Counted_Hits = sum(Correct_Left_Button_Press),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Task_Available,
      values_from = Counted_Hits,
      names_prefix = "Counted_Hits_",
      values_fill = 0
    ),
  by = "DaqName"
)

### Add and modify predictors

non_aug$Blinking <- non_aug$LPupil_Diameter == 0 &
  non_aug$RPupil_Diameter == 0

# Get saccades, fixations, and blinks
saccade_sens <- 15 # Detection threshold parameter for `saccades` event detection
non_saccades <- detect.fixations(
  non_aug %>%
    select(
      DaqName,
      Gaze_Yaw,
      Gaze_Pitch,
      X
    ) %>%
    rename(
      trial = DaqName,
      x = Gaze_Yaw,
      y = Gaze_Pitch,
      time = X
    ),
  lambda = saccade_sens
)

# Add column to `non_aug` for eye event type
non_aug$Eye_Event <- factor(
  "saccade",
  levels = c("saccade", "blink", "fixation", "too short")
)
indices <- mapply(seq, non_saccades$start, non_saccades$end)
flattened_indices <- unlist(indices)
non_aug$Eye_Event[flattened_indices] <- rep(
  as.character(non_saccades$event),
  times = lengths(indices)
)
non_aug[is.na(non_aug$Eye_Event), "Eye_Event"] <- "saccade"
rm(non_saccades, indices, flattened_indices)

### Aggregate matrix

non_aug <- non_aug %>%
  group_by(consecutive_id(Task_Available)) %>%
  mutate(Task_Start_X = first(X)) %>%
  ungroup()

# Long-form wrt pre/post pairs

# Sensitivity of lateral deviation
lat_dev_sens <- 4

# Cleaned button press frames
press_events <- non_aug %>%
  filter(Correct_Left_Button_Press == 1) %>%
  select(
    DaqName,
    Subject,
    Drive,
    X,
    Frame_Num,
    Left_Button_Press,
    Task_Start_X,
  )

# Window of frames around each press
task_windows <- press_events %>%
  rowwise() %>%
  reframe(
    DaqName = DaqName,
    Subject = Subject,
    Drive = Drive,
    # Want the frame number within the run of the simulator
    Press_Frame = Frame_Num,
    Frame_Index = X,
    Frame = c(
      Task_Start_X:(X - 1),
      X:(2 * X - Task_Start_X - 1)
    ), # Include frame of press since press occurred before recording
    Phase = c(
      rep("Pre_Press", X - Task_Start_X),
      rep("Post_Press", X - Task_Start_X)
    )
  )

# Find the boundaries of every task block per Daq run
task_boundaries <- task_windows %>%
  group_by(DaqName, Subject, Drive, Press_Frame, Frame_Index) %>%
  summarise(
    Task_Start_X = min(Frame),
    Segment_Length = n() / 2,
    .groups = "drop"
  ) %>%
  arrange(DaqName, Task_Start_X) %>%
  group_by(DaqName) %>%
  mutate(
    Control_Start = Task_Start_X - Segment_Length,
    Control_End = Task_Start_X - 1
  ) %>%
  ungroup() %>%
  filter(!is.na(Control_Start) & Control_End >= Control_Start)

# Expand these boundaries into actual individual frame rows
control_windows <- task_boundaries %>%
  rowwise() %>%
  reframe(
    DaqName = DaqName,
    Subject = Subject,
    Drive = Drive,
    Press_Frame = Press_Frame, # Linking it to the upcoming press
    Frame_Index = Frame_Index,
    Frame = Control_Start:Control_End,
    Phase = "Control"
  )
rm(task_boundaries)

# Combine with original Pre/Post task windows with the new control windows
all_windows <- bind_rows(task_windows, control_windows) %>%
  arrange(DaqName, Frame) %>%
  mutate(Phase = factor(Phase, levels = c(
    "Control",
    "Pre_Press",
    "Post_Press"
  )))
rm(control_windows, task_windows)

# Pull in raw eye data for ALL windows (Tasks and Controls)
task_windows_with_control <- all_windows %>%
  inner_join(
    select(non_aug, -DaqName),
    by = c("Subject", "Drive", "Frame" = "X")
  )
rm(all_windows)

# Aggregate into one row per window with the control phases
task_matrix <- task_windows_with_control %>%
  group_by(DaqName, Subject, Drive, Press_Frame, Frame_Index, Phase) %>%
  summarise(
    LPupil_Diameter_mean = mean(LPupil_Diameter, na.rm = TRUE),
    RPupil_Diameter_mean = mean(RPupil_Diameter, na.rm = TRUE),
    Avg_Pupil_Diameter_mean = mean(Avg_Pupil_Diameter, na.rm = TRUE),
    Gaze_Pitch_mean = mean(Gaze_Pitch, na.rm = TRUE),
    Gaze_Yaw_mean = mean(Gaze_Yaw, na.rm = TRUE),
    Deviation_frames_prop = sum(
      abs(Vehicle_Lat_Dev) > Road_Width / lat_dev_sens
    ) / n(),
    BAC = Start_BAC[1],
    KSS = KSS_Score[1],
    Blink_Count = max(Blink_Counter) - min(Blink_Counter),
    # Per close is the proportion that **both eyes are closed simultaneously**
    Per_Close = mean(Blinking, na.rm = TRUE),
    Lat_Dev_SD = sd(Vehicle_Lat_Dev, na.rm = TRUE),
    Speed_SD = sd(Vehicle_Speed, na.rm = TRUE),
    Braking_Events = sum(lag(Brake_Pedal_Force) == 0 & Brake_Pedal_Force > 0, na.rm = TRUE),
    Reaction_Frames = n(),
    Saccade_Count = sum(Eye_Event == "saccade" & lag(Eye_Event) != "saccade", na.rm = TRUE),
    Saccade_Prop = mean(Eye_Event == "saccade", na.rm = TRUE),
    .groups = "drop"
  )
rm(task_windows_with_control)

write.csv(task_matrix, "./data/task_matrix.csv", row.names = FALSE)

# Make wide-form version for paired differences
wideform_task_matrix <- task_matrix %>%
  pivot_wider(
    names_from = Phase,
    values_from = setdiff(
      tail(colnames(task_matrix), -5),
      c("BAC", "KSS")
    )
  ) %>%
  mutate(
    Avg_Pupil_Diameter_mean_diff = (
      Avg_Pupil_Diameter_mean_Post_Press - Avg_Pupil_Diameter_mean_Pre_Press
    ),
    RPupil_Diameter_mean_diff = (
      RPupil_Diameter_mean_Post_Press - RPupil_Diameter_mean_Pre_Press
    ),
    LPupil_Diameter_mean_diff = (
      LPupil_Diameter_mean_Post_Press - LPupil_Diameter_mean_Pre_Press
    ),
    Gaze_Yaw_mean_diff = (
      Gaze_Yaw_mean_Post_Press - Gaze_Yaw_mean_Pre_Press
    ),
    Gaze_Pitch_mean_diff = (
      Gaze_Pitch_mean_Post_Press - Gaze_Pitch_mean_Pre_Press
    ),
    Deviation_Frames_Prop_Diff = (
      Deviation_frames_prop_Post_Press - Deviation_frames_prop_Pre_Press
    ),
    Lat_Dev_SD_diff = (
      Lat_Dev_SD_Post_Press - Lat_Dev_SD_Pre_Press
    ),
    Speed_SD_diff = (
      Speed_SD_Post_Press - Speed_SD_Pre_Press
    ),
    Braking_Events_diff = (
      Braking_Events_Post_Press - Braking_Events_Pre_Press
    ),
  )

### Write environment
save.image("./data/processed_data.RData")
