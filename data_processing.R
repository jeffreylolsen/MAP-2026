library(tidyverse)
library(dplyr)
library(readxl)
library(saccades)
library(arrow)
library(data.table)
library(this.path)
library(runner)
library(qs2)

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
non_aug$`rleid(Task_Available)` <- NULL

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

non_aug <- non_aug %>% mutate(
  Blinking = Left_Eyelid_Closed & Right_Eyelid_Closed
)

# Get saccades, fixations, and blinks
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
    )
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
non_aug$`consecutive_id(Task_Available)` <- NULL

# Sensitivity of lateral deviation
lat_dev_sens <- 4

# Number of frames to add to task segments on top of reaction time
frame_length <- 300

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
    Frame = Task_Start_X:(X + frame_length - 1),
    Phase = rep("Task", X - Task_Start_X + frame_length)
  )

# Find the boundaries of every task block per Daq run
task_boundaries <- task_windows %>%
  group_by(DaqName, Subject, Drive, Press_Frame, Frame_Index) %>%
  summarize(
    Task_Start_X = min(Frame),
    Segment_Length = n(),
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
non_aug$Task_Start_X <- NULL

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
    "Task"
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
  summarize(
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
    Braking_Events = sum(lag(Brake_Pedal_Force) == 0 & Brake_Pedal_Force > 0,
      na.rm = TRUE
    ),
    Brake_Force_mean = mean(Brake_Pedal_Force, na.rm = TRUE),
    Saccade_Count = sum(Eye_Event == "saccade" & lag(Eye_Event) != "saccade",
      na.rm = TRUE
    ),
    Saccade_Prop = mean(Eye_Event == "saccade", na.rm = TRUE),
    Frame_Count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    BAC_cent = BAC - median(BAC, na.rm = TRUE),
    KSS_cent = KSS - median(KSS, na.rm = TRUE)
  ) %>%
  group_by(Frame_Index) %>%
  mutate(Reaction_Frames = Frame_Count[1] - frame_length) %>%
  ungroup()
rm(task_windows_with_control)

write.csv(task_matrix, "./data/task_matrix.csv", row.names = FALSE)

# Make wide-form version for paired differences
wideform_task_matrix <- task_matrix %>%
  pivot_wider(
    names_from = Phase,
    values_from = setdiff(
      tail(colnames(task_matrix), -5),
      c("BAC", "KSS", "BAC_cent", "KSS_cent", "Reaction_Frames")
    )
  ) %>%
  mutate(
    Avg_Pupil_Diameter_mean_diff = (
      Avg_Pupil_Diameter_mean_Task - Avg_Pupil_Diameter_mean_Control
    ),
    RPupil_Diameter_mean_diff = (
      RPupil_Diameter_mean_Task - RPupil_Diameter_mean_Control
    ),
    LPupil_Diameter_mean_diff = (
      LPupil_Diameter_mean_Task - LPupil_Diameter_mean_Control
    ),
    Gaze_Yaw_mean_diff = (
      Gaze_Yaw_mean_Task - Gaze_Yaw_mean_Control
    ),
    Gaze_Pitch_mean_diff = (
      Gaze_Pitch_mean_Task - Gaze_Pitch_mean_Control
    ),
    Deviation_Frames_Prop_Diff = (
      Deviation_frames_prop_Task - Deviation_frames_prop_Control
    ),
    Lat_Dev_SD_diff = (
      Lat_Dev_SD_Task - Lat_Dev_SD_Control
    ),
    Speed_SD_diff = (
      Speed_SD_Task - Speed_SD_Control
    ),
    Braking_Events_diff = (
      Braking_Events_Task - Braking_Events_Control
    ),
  )

### Write environment
qs_save(as.list(environment()), "./data/processed_data.qs")
