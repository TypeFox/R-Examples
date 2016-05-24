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

library(splines)
nl_model <- addNonlinearity(model, "age", "ns", flex_param = 2:7,
                            workers = FALSE)  
if(!all.equal(model, nl_model))
  stop("Failed check coxph with ns")

nl_model <- addNonlinearity(model, "age", "ns", flex_param = 2:7, sig_level = .7,
                            workers = FALSE)  
if(length(all.equal(model, nl_model)) == 1)
  stop("Failed check coxph with ns with sensitivity increased")

nl_model <- addNonlinearity(model, "age", "pspline", flex_param = "Asdasdsadasda",
                            workers = FALSE)  
if(!all.equal(model, nl_model))
  stop("Failed check coxph with pspline")

nl_model <- addNonlinearity(model, "age", "pspline", flex_param = "Asdasdsadasda",sig_level = .7,
                            workers = FALSE)  
if(length(all.equal(model, nl_model)) == 1)
  stop("Failed check coxph with pspline with sensitivity increased")
