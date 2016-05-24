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

# remove trials with > 25% of trackloss
response_window_clean <- clean_by_trackloss(data = response_window,
                                            trial_prop_thresh = .25)

# create Target condition column
response_window_clean$Target <- as.factor( ifelse(test = grepl('(Spoon|Bottle)', response_window_clean$Trial), 
                                       yes = 'Inanimate', 
                                       no  = 'Animate') )

## ---- warning=FALSE------------------------------------------------------
(data_summary <- describe_data(response_window_clean, 
                               describe_column='Animate', group_columns=c('Target','ParticipantName')))
plot(data_summary)

## ---- warning=FALSE------------------------------------------------------
# aggregate by subject across the response window
response_window_agg_by_sub <- make_time_window_data(response_window_clean, 
                                             aois='Animate',
                                             predictor_columns=c('Target','Age','MCDI_Total'),
                                             summarize_by = "ParticipantName")

# take a quick peek at data
plot(response_window_agg_by_sub, predictor_columns="Target", dv = "ArcSin")

# show condition means
describe_data(response_window_agg_by_sub, describe_column = "ArcSin", group_columns = "Target")

# simple paired t-test between conditions
t.test(ArcSin ~ Target, data= response_window_agg_by_sub, paired=TRUE)

## ---- warning=FALSE------------------------------------------------------
# you should almost always sum-code and center your predictors when performing regression analyses
response_window_agg_by_sub$AgeC <- response_window_agg_by_sub$Age - mean(response_window_agg_by_sub$Age)
response_window_agg_by_sub$MCDI_TotalC <- response_window_agg_by_sub$MCDI_Total - mean(response_window_agg_by_sub$MCDI_Total)

model <- lm(ArcSin ~ Target*AgeC*MCDI_TotalC, data=response_window_agg_by_sub)
broom::tidy(model)

## ---- warning=FALSE------------------------------------------------------
response_window_agg <- make_time_window_data(response_window_clean, 
                                             aois='Animate', 
                                             predictor_columns=c('Target','Age','MCDI_Total'))

# sum-code and center predictors
response_window_agg$TargetC <- ifelse(response_window_agg$Target == 'Animate', .5, -.5)
response_window_agg$TargetC <- as.numeric(scale(response_window_agg$TargetC, center=TRUE, scale=FALSE))

# mixed-effects linear model on subject*trial data
model_time_window <- lmer(Elog ~ TargetC + (1 + TargetC | Trial) + (1 | ParticipantName), 
                          data = response_window_agg, REML = FALSE)
# cleanly show important parts of model (see `summary()` for more)
(est <- broom::tidy(model_time_window, effects="fixed"))

# use model comparison to attain p-values
drop1(model_time_window,~.,test="Chi")

## ---- warning=FALSE------------------------------------------------------
condition_estimate <- with(est, 
                           c(estimate[term=="(Intercept)"] + estimate[term=="TargetC"] / 2,
                             estimate[term=="(Intercept)"] - estimate[term=="TargetC"] / 2))

## ---- warning=FALSE------------------------------------------------------
exp(condition_estimate)/(1+exp(condition_estimate))

## ---- warning=FALSE------------------------------------------------------
plot(model_time_window)

## ---- warning=FALSE------------------------------------------------------
model_time_window_logit <- lmer(LogitAdjusted ~ TargetC + (1 + TargetC | Trial) + (1 | ParticipantName), 
                          data = response_window_agg, REML = FALSE)
plot(model_time_window_logit)
drop1(model_time_window_logit,~.,test="Chi")
est_logit <- broom::tidy(model_time_window_logit, effects="fixed")
condition_estimate_logit <- with(est_logit, 
                           c(estimate[term=="(Intercept)"] + estimate[term=="TargetC"] / 2,
                             estimate[term=="(Intercept)"] - estimate[term=="TargetC"] / 2))
exp(condition_estimate_logit)/(1+exp(condition_estimate_logit))

## ---- warning=FALSE------------------------------------------------------
response_window_agg$AgeC <- response_window_agg$Age - mean(response_window_agg$Age)
response_window_agg$MCDI_TotalC <- response_window_agg$MCDI_Total - mean(response_window_agg$MCDI_Total)

model_time_window_add_predictors <- lmer(Elog ~ TargetC*AgeC*MCDI_TotalC + (1 + TargetC | Trial) + (1 | ParticipantName), 
              data = response_window_agg, REML = FALSE)
# cleanly show important parts of model (see `summary()` for more)
broom::tidy(model_time_window_add_predictors, effects="fixed")

# use model comparison to attain p-values
drop1(model_time_window_add_predictors,~.,test="Chi")

