library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(glue)
library(mgcv)
library(data.table)
library(this.path)
library(glmmTMB)
library(sjPlot)
library(performance)
library(visreg)
library(modelsummary)
library(emmeans)
library(plotly)
library(qs2)

# Read in processed data
setwd(this.dir())

if (!file.exists("./data/processed_data.qs") |
  file.mtime("./data/processed_data.qs") <
    file.mtime("./data_processing.R")) {
  local(source("./data_processing.R", local = TRUE))
}
qs_readm("./data/processed_data.qs")

models <- list()

for (frame_length in as.character(c(180, 300, 600))) {
  for (var in c("task_matrix", "wideform_task_matrix")) {
    assign(var, get(glue("{var}_{frame_length}")))
  }
  rm(var)

  # Preliminary models
  models[["pup"]][[frame_length]] <- lmer(
    Avg_Pupil_Diameter_mean_diff ~
      BAC +
      KSS_cent +
      (1 | Subject),
    data = wideform_task_matrix
  )

  models[["sdlp"]][[frame_length]] <- lmer(
    Lat_Dev_SD_diff ~
      BAC +
      KSS_cent +
      (1 | Subject),
    data = wideform_task_matrix
  )

  models[["speed"]][[frame_length]] <- lmer(
    Speed_SD_diff ~ BAC +
      KSS_cent +
      (1 | Subject),
    data = wideform_task_matrix
  )

  models[["dev"]][[frame_length]] <- lmer(
    Deviation_Frames_prop_diff ~ BAC +
      KSS_cent +
      (1 | Subject),
    data = wideform_task_matrix
  )

  # Pupil diameter difference
  models[["pupil_diameter_diff"]][[frame_length]] <- lmer(
    Avg_Pupil_Diameter_mean_diff ~
      BAC +
      KSS_cent +
      as.numeric(Drive) +
      (1 | Subject),
    data = wideform_task_matrix
  )

  # Saccades ZIP model
  models[["saccade_zip"]][[frame_length]] <- glmmTMB(
    Saccade_Count ~ Phase + BAC + KSS_cent + Road_Surface +
      (1 | Subject) + (1 | Frame_Index),
    ziformula = ~1, # What predicts being in the "zero" group?
    family = poisson,
    data = task_matrix
  )

  # Saccades logit model
  models[["saccade_logit"]][[frame_length]] <- glmer(
    data = task_matrix,
    formula =
      (Saccade_Count > 0) ~
        BAC +
        KSS_cent +
        Road_Surface +
        BAC * KSS_cent +
        # BAC * Phase +
        as.numeric(Drive) +
        Phase +
        (1 | Subject) + (1 | Frame_Index),
    family = "binomial"
  )

  # Reaction time models
  models[["rt_driving"]][[frame_length]] <- lmer(
    (Reaction_Frames * (1000 / 60)) ~
      BAC +
      KSS_cent +
      Road_Surface +
      as.numeric(Drive) + # Drive term to adjust for learning effect
      (1 | Subject),
    data = wideform_task_matrix
  )
  models[["reaction_vs_latdev_by_drive"]][[frame_length]] <-
    lmer(
      data = wideform_task_matrix,
      formula = (
        (Reaction_Frames * (1000 / 60)) ~
          Lat_Dev_SD_Control +
          KSS_cent +
          BAC +
          Road_Surface +
          as.numeric(Drive) + # To adjust for learning effect
          Lat_Dev_SD_Control * BAC +
          (1 | Subject)
      )
    )

  # SDLP model
  models[["lat_dev_phase"]][[frame_length]] <- lmer(
    data = task_matrix,
    formula = Lat_Dev_SD ~
      Lat_Dev_abs_init +
      KSS_cent +
      BAC +
      BAC * KSS_cent +
      BAC * Phase +
      Road_Surface +
      as.numeric(Drive) +
      Phase +
      (1 | Subject) + (1 | Frame_Index)
  )
  models[["lat_dev_phase_visreg"]][[frame_length]] <- visreg(
    models[["lat_dev_phase"]][[frame_length]],
    "BAC",
    by = "KSS_cent",
    plot = FALSE
  )

  # Lane departure model
  models[["lane_departure"]][[frame_length]] <- glmer(
    data = task_matrix,
    formula = (
      Lane_Departure ~
        KSS_cent +
        BAC +
        BAC * KSS_cent +
        BAC * Phase +
        Road_Surface +
        as.numeric(Drive) +
        Phase +
        (1 | Subject) + (1 | Frame_Index)
    ),
    family = "binomial"
  )

  # Speed models
  models[["speed_sd"]][[frame_length]] <- lmer(
    data = task_matrix,
    formula = Speed_SD ~
      KSS_cent +
      BAC +
      BAC * KSS_cent +
      BAC * Phase +
      Road_Surface +
      as.numeric(Drive) +
      Phase +
      (1 | Subject) + (1 | Frame_Index)
  )
  models[["speed_mean"]][[frame_length]] <- lmer(
    data = task_matrix,
    formula = Speed_mean ~
      KSS_cent +
      BAC +
      BAC * KSS_cent +
      BAC * Phase +
      Road_Surface +
      as.numeric(Drive) +
      Phase +
      (1 | Subject) + (1 | Frame_Index)
  )

  # Per Close model
  models[["per_close"]][[frame_length]] <- glmmTMB(
    # A single observation of Per_Close = 1,
    # not likely to be seen frequently in population data
    data = filter(task_matrix, Per_Close != 1),
    formula = (
      Per_Close ~
        KSS_cent +
        BAC +
        Phase +
        BAC * KSS_cent +
        BAC * Phase +
        Road_Surface +
        as.numeric(Drive) +
        (1 | Subject) + (1 | Frame_Index)
    ),
    ziformula = ~1,
    family = beta_family()
  )

  # Brake force model
  models[["brake_force"]][[frame_length]] <- lmer(
    data = task_matrix,
    formula = (
      Brake_Force_mean ~
        Phase +
        KSS_cent +
        BAC +
        BAC * Phase +
        BAC * KSS_cent +
        Road_Surface +
        as.numeric(Drive) +
        (1 | Subject) + (1 | Frame_Index)
    )
  )
}
rm(frame_length, task_matrix, wideform_task_matrix)

# Write environment
modeling_path <- this.path()
qs_save(as.list(environment()), "./data/final_data.qs")
