## ----results='hide'------------------------------------------------------
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
# aggregate across trials within subjects in time analysis
response_time <- make_time_sequence_data(response_window_clean, time_bin_size = 100, 
                                 predictor_columns = c("Target"),
                                 aois = "Animate"
                            )

# visualize time results
plot(response_time, predictor_column = "Target") + 
  theme_light() +
  coord_cartesian(ylim = c(0,1))

## ---- warning=FALSE------------------------------------------------------
# generate dataframe summarizing values of each vector
timecodes <- unique(response_time[, c('ot1','ot2','ot3','ot4','ot5','ot6','ot7')])
timecodes$num <- 1:nrow(timecodes)

ggplot(timecodes, aes(x=num, y=ot1)) +
               geom_line() +
               geom_line(aes(y=ot2), color='red') +    # quadratic
               geom_line(aes(y=ot3), color='blue') +   # cubic
               geom_line(aes(y=ot4), color='green') +  # quartic
               geom_line(aes(y=ot5), color='purple') + # quintic 
               geom_line(aes(y=ot6), color='yellow') + # sextic
               geom_line(aes(y=ot7), color='pink') +   # septic
               scale_x_continuous(name="") +
               scale_y_continuous(name="")

round(cor(timecodes[, c(1:7)]),5)

## ---- warning=FALSE------------------------------------------------------
# sum-code and center our predictor:
response_time$TargetC <- ifelse(response_time$Target == 'Animate', .5, -.5)
response_time$TargetC <- as.numeric(scale(response_time$TargetC, center=TRUE, scale=FALSE))

# Construct model
model_time_sequence <- lmer(Elog ~ TargetC*(ot1) + (1 + ot1 | Trial) + (1 + ot1 | ParticipantName), 
              data = response_time, REML = FALSE)

# cleanly show important parts of model (see `summary()` for more)
broom::tidy(model_time_sequence, effects = "fixed")
drop1(model_time_sequence, ~., test="Chi")

## ---- warning=FALSE------------------------------------------------------
plot(response_time, predictor_column = "Target", dv = "Elog", model = model_time_sequence) +
  theme_light()

## ---- warning=FALSE------------------------------------------------------
model_time_sequence <- lmer(Elog ~ TargetC*(ot1 + ot2 + ot3 + ot4) + (1 | Trial) + (1 | ParticipantName), 
                            data = response_time, REML = FALSE)

# commented out because this many random slopes takes too long to fit for a vignette
# model_time_sequence <- lmer(Elog ~ TargetC*(ot1 + ot2 + ot3 + ot4) + 
#                               (1 + ot1 + ot2 + ot3 + ot4 | Trial) + (1 + ot1 + ot2 + ot3 + ot4 | ParticipantName), 
#                             data = response_time, REML = FALSE)

# cleanly show important parts of model (see `summary()` for more)
broom::tidy(model_time_sequence, effects = "fixed")
drop1(model_time_sequence, ~., test="Chi")

plot(response_time, predictor_column = "Target", dv = "Elog", model = model_time_sequence) +
  theme_light()

