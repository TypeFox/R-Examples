## ----results='hide'------------------------------------------------------
set.seed(42)

library("Matrix")
library("lme4")
library("ggplot2")
library("eyetrackingR")

data("word_recognition")
data <- make_eyetrackingr_data(word_recognition, 
                       participant_column = "ParticipantName",
                       trial_column = "Trial",
                       time_column = "TimeFromTrialOnset",
                       trackloss_column = "TrackLoss",
                       aoi_columns = c('Animate','Inanimate'),
                       treat_non_aoi_looks_as_missing = TRUE
)

## ---- eval=FALSE---------------------------------------------------------
#  animate_aoi <- read.csv("./interest_areas_for_animate_aoi.csv")
#  
#  #            Trial Left Top Right Bottom
#  # 1   FamiliarBird  500 100   900    500
#  # 2 FamiliarBottle  400 200   800    600
#  # 3    FamiliarCow  500 300   900    700
#  # 4    FamiliarDog  300 100   700    500
#  # 5  FamiliarHorse  500 200   900    600
#  # 6  FamiliarSpoon  350 300   750    700
#  
#  data <- add_aoi(data = data, aoi_dataframe = animate_aoi,
#                 x_col = "GazeX", y_col = "GazeY",
#                 aoi_name = "Animate",
#                 x_min_col = "Left", x_max_col = "Right", y_min_col = "Top", y_max_col = "Bottom")

## ------------------------------------------------------------------------
table(data$Animate)
table(is.na(data$Animate)) # if all TRUE, then something went wrong.

## ---- echo=FALSE---------------------------------------------------------
data$Message <- with(data, ifelse(TimeFromTrialOnset==0, "TrialStart", "."))
data$ResponseWindowStart <- 15500 

## ------------------------------------------------------------------------
data <- subset_by_window(data, window_start_msg = "TrialStart", msg_col = "Message", rezero= TRUE)

## ------------------------------------------------------------------------
response_window <- subset_by_window(data, window_start_col = "ResponseWindowStart", rezero= FALSE, remove= TRUE)

## ------------------------------------------------------------------------
response_window <- subset_by_window(response_window, window_end_time = 21000, rezero= FALSE, remove= TRUE)

## ---- warning=FALSE------------------------------------------------------
# analyze amount of trackloss by subjects and trials
(trackloss <- trackloss_analysis(data = response_window))

response_window_clean <- clean_by_trackloss(data = response_window, trial_prop_thresh = .25)

## ---- warning=FALSE------------------------------------------------------
trackloss_clean <- trackloss_analysis(data = response_window_clean)

(trackloss_clean_subjects <- unique(trackloss_clean[, c('ParticipantName','TracklossForParticipant')]))

## ---- warning=FALSE------------------------------------------------------
# get mean samples contributed per trials, with SD
mean(1 - trackloss_clean_subjects$TracklossForParticipant)
sd(1- trackloss_clean_subjects$TracklossForParticipant)

## ---- warning=FALSE------------------------------------------------------
# look at the NumTrials column
(final_summary <- describe_data(response_window_clean, 'Animate', 'ParticipantName'))

## ---- warning=FALSE------------------------------------------------------
mean(final_summary$NumTrials)
sd(final_summary$NumTrials)

## ---- warning=FALSE------------------------------------------------------
response_window_clean$Target <- as.factor( ifelse(test = grepl('(Spoon|Bottle)', response_window_clean$Trial), 
                                       yes = 'Inanimate', 
                                       no  = 'Animate') )

