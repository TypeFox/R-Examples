library(Greg)
library(magrittr)
data("melanoma", package = "boot")
melanoma$status %<>% 
  factor(
    levels=1:3,
    labels=c("melanoma-specific death", "alive", "other death"))
melanoma$ulcer %<>% 
  factor(
    levels=0:1,
    labels=c("Absent", "Present")
  )
melanoma$time <-
  melanoma$time/365.25

library(survival)
regular_model <- coxph(Surv(time, status == "melanoma-specific death") ~
                         age + sex + year + thickness + ulcer,
                       data = melanoma,
                       x = TRUE, y = TRUE)
spl_melanoma <-
  melanoma %>% 
  timeSplitter(by = .1,
               event_var = "status",
               event_start_status = "alive",
               time_var = "time",
               time_related_vars = c("age", "year"))

interval_model <-
  update(regular_model, 
         Surv(Start_time, Stop_time, status == "melanoma-specific death") ~ .,
         data = spl_melanoma)
mismatch <- abs(sum(coef(interval_model) - coef(regular_model)))
if (mismatch > 10^-10)
  stop("Failed to match interval with regular cox model.",
       " Total coefficient difference = ", mismatch,
       "\n Regular: ", paste(txtRound(coef(regular_model), 3),
                             collapse = ", "),
       "\n Interval: ", paste(txtRound(coef(interval_model), 3),
                              collapse = ", "))
