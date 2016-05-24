library(Greg)
library(magrittr)
data("melanoma", package = "boot", envir = environment())

library(dplyr)
library(magrittr)
melanoma %<>% 
  mutate(status = factor(status,
                         levels = 1:3,
                         labels = c("Died from melanoma", 
                                    "Alive", 
                                    "Died from other causes")),
         ulcer = factor(ulcer,
                        levels = 0:1,
                        labels = c("Absent", "Present")),
         time = time/365.25, # All variables should be in the same time unit
         sex = factor(sex,
                      levels = 0:1,
                      labels = c("Female", "Male")))

library(survival)
model <- coxph(Surv(time, status == "Died from melanoma") ~ sex + age,
               data = melanoma)

nl_model <- addNonlinearity(model, "age", 
                            spline_fn = "pspline", 
                            verbal = TRUE,
                            workers = FALSE)
# Note that there is no support for nonlinearity in this case
