## ----results="hide"------------------------------------------------------
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
response_window_clean <- clean_by_trackloss(data = response_window,
                                            trial_prop_thresh = .25)

# create Target condition column
response_window_clean$Target <- as.factor( ifelse(test = grepl('(Spoon|Bottle)', response_window_clean$Trial), 
                                       yes = 'Inanimate', 
                                       no  = 'Animate') )

## ---- warning=FALSE------------------------------------------------------
response_time <- make_time_sequence_data(response_window_clean,
                                  time_bin_size = 100, 
                                  predictor_columns = c("Target"),
                                  aois = "Animate",
                                  summarize_by = "ParticipantName" )

# visualize timecourse
plot(response_time, predictor_column = "Target") + 
  theme_light() +
  coord_cartesian(ylim = c(0,1))

## ---- warning=FALSE------------------------------------------------------
tb_analysis <- analyze_time_bins(data = response_time, predictor_column = "Target", test = "t.test", alpha = .05)
plot(tb_analysis, type = "estimate") + theme_light()
summary(tb_analysis)

## ------------------------------------------------------------------------
alpha <- .05
num_time_bins <- nrow(tb_analysis)
(prob_no_false_alarm_per_bin <- 1-alpha)
(prob_no_false_alarm_any_bin <- prob_no_false_alarm_per_bin^num_time_bins)
(prob_at_least_one_false_alarm <- 1-prob_no_false_alarm_any_bin)

## ------------------------------------------------------------------------
alpha <- .05 / num_time_bins
(prob_no_false_alarm_per_bin <- 1-alpha)
(prob_no_false_alarm_any_bin <- prob_no_false_alarm_per_bin^num_time_bins)
(prob_at_least_one_false_alarm <- 1-prob_no_false_alarm_any_bin)

## ---- warning=FALSE------------------------------------------------------
tb_analysis_bonf <- analyze_time_bins(data = response_time, predictor_column = "Target", test = "t.test", alpha = .05,
                                 p_adjust_method = "bonferroni")
plot(tb_analysis_bonf) + theme_light()
summary(tb_analysis_bonf)

## ---- warning=FALSE------------------------------------------------------
tb_analysis_holm <- analyze_time_bins(data = response_time, predictor_column = "Target", test = "t.test", alpha = .05,
                                 p_adjust_method = "holm")
plot(tb_analysis_holm) + theme_light()
summary(tb_analysis_holm)

## ---- warning=FALSE------------------------------------------------------
tb_bootstrap <- analyze_time_bins(response_time, predictor_column = 'Target', test= 'boot_splines', 
                                  within_subj = TRUE, bs_samples = 1000, alpha = .05)
plot(tb_bootstrap) + theme_light()
summary(tb_bootstrap)

## ---- warning=FALSE------------------------------------------------------
tb_bootstrap_bonf <- analyze_time_bins(response_time, predictor_column = 'Target', test= 'boot_splines', 
                                  within_subj = TRUE, alpha = .05/num_time_bins)
plot(tb_bootstrap_bonf) + theme_light()
summary(tb_bootstrap_bonf)

## ---- warning=FALSE------------------------------------------------------
num_sub = length(unique((response_window_clean$ParticipantName)))
threshold_t = qt(p = 1 - .05/2, 
                 df = num_sub-1) # pick threshold t based on alpha = .05 two tailed

## ---- warning=FALSE------------------------------------------------------
df_timeclust <- make_time_cluster_data(response_time, 
                                      test= "t.test", paired=TRUE,
                                      predictor_column = "Target", 
                                      threshold = threshold_t) 
plot(df_timeclust) +  ylab("T-Statistic") + theme_light()
summary(df_timeclust)

## ---- warning=FALSE------------------------------------------------------
clust_analysis <- analyze_time_clusters(df_timeclust, within_subj=TRUE, paired=TRUE,
                                        samples=150) # in practice, you should use a lot more
plot(clust_analysis) + theme_light()

## ---- warning=FALSE------------------------------------------------------
summary(clust_analysis)

## ---- warning=FALSE------------------------------------------------------
response_time_between <- make_time_sequence_data(response_window_clean,
                                  time_bin_size = 100, 
                                  predictor_columns = c("Sex", "MCDI_Total"),
                                  aois = "Animate",
                                  summarize_by = "ParticipantName" )

df_timeclust_between <- make_time_cluster_data(response_time_between, 
                                      test= "lm",
                                      predictor_column = "MCDI_Total", 
                                      threshold = threshold_t) 
plot(df_timeclust_between) +  ylab("T-Statistic") + theme_light()
summary(df_timeclust_between)

## ------------------------------------------------------------------------
clust_analysis_between <- analyze_time_clusters(df_timeclust_between, within_subj = FALSE, 
                                        samples=100) # in practice, you should use a lot more
plot(clust_analysis_between) + theme_light()
summary(clust_analysis_between)

