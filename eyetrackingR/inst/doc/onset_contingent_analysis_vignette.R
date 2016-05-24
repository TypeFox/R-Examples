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

# subset to response window post word-onset
response_window <- subset_by_window(data, 
                                    window_start_time = 15500, 
                                    window_end_time = 21000, 
                                    rezero = FALSE)

# analyze amount of trackloss by subjects and trials
(trackloss <- trackloss_analysis(data = response_window))

# remove trials with > 25% of trackloss
response_window_clean <- clean_by_trackloss(data = response_window, trial_prop_thresh = .25)

# create Target condition column
response_window_clean$Target <- as.factor( ifelse(test = grepl('(Spoon|Bottle)', response_window_clean$Trial), 
                                       yes = 'Inanimate', 
                                       no  = 'Animate') )

## ---- warning=FALSE------------------------------------------------------
# recode AOIs to target & distractor
response_window_clean$TrialTarget <- ifelse(test = response_window_clean$Target == 'Animate', 
                                            yes = response_window_clean$Animate, 
                                            no = response_window_clean$Inanimate)
response_window_clean$TrialDistractor <- ifelse(test = response_window_clean$Target == 'Animate', 
                                                yes = response_window_clean$Inanimate, 
                                                no = response_window_clean$Animate)

## ---- warning=FALSE------------------------------------------------------
onsets <- make_onset_data(response_window_clean, onset_time = 15500, 
                          fixation_window_length = 1, target_aoi='TrialTarget')
# participants' ability to orient to the trial target overall:
plot(onsets) + theme(legend.text=element_text(size=5))

## ---- warning=FALSE------------------------------------------------------
# participants' ability to orient to the trial target, split by which target:
plot(onsets, predictor_columns = "Target") + theme(legend.text=element_text(size=6))

## ---- warning=FALSE------------------------------------------------------
# we can also visualize numeric predictors:
plot(onsets, predictor_columns = "MCDI_Total") + theme(legend.text=element_text(size=6))

## ---- warning= FALSE-----------------------------------------------------
onset_switches <- make_switch_data(onsets, predictor_columns = "Target")

# visualize subject's switch times
plot(onset_switches, predictor_columns = c("Target"))

# center predictor:
onset_switches$FirstAOIC <- ifelse(onset_switches$FirstAOI == 'TrialTarget', .5, -.5)
onset_switches$FirstAOIC <- scale(onset_switches$FirstAOIC, center=TRUE, scale=FALSE) 
onset_switches$TargetC <- ifelse(onset_switches$Target == 'Animate', .5, -.5)
onset_switches$TargetC <- scale(onset_switches$TargetC, center=TRUE, scale=FALSE) 

# build model:
model_switches <- lmer(FirstSwitch ~ FirstAOIC*TargetC + 
                (1 | Trial) + (1 | ParticipantName), data=onset_switches, REML=FALSE)
# cleanly show important parts of model (see `summary()` for more)
broom::tidy(model_switches, effects="fixed")
drop1(model_switches,~.,test="Chi")

