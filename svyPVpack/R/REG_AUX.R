########################################################################
# REG AUX FUNCTIONS
########################################################################

glm.dummy.no.std <- function(formula, 
                             design_dummy_4_glm.dummy.no.std,
                             family_dummy_4_glm.dummy.no.std){
  design_dummy_4_glm.dummy.no.std <<- design_dummy_4_glm.dummy.no.std
  family_dummy_4_glm.dummy.no.std <<- family_dummy_4_glm.dummy.no.std
  # estimate mod
  mod <- svyglm(formula, 
                design_dummy_4_glm.dummy.no.std, 
                family=family_dummy_4_glm.dummy.no.std)
  res <- summary(mod)
  # goodness of fit
  test.term <- as.character(formula[3])
  #.e <- environment()
  test <- regTermTest(mod, as.formula(paste("~", test.term)),
                      method="LRT")
  rm(list=c("design_dummy_4_glm.dummy.no.std", 
            "family_dummy_4_glm.dummy.no.std"),
     envir = .GlobalEnv)
  return(list(coef = res$coef, 
              mod.fit = cbind(test$chi, t(test$lambda), test$ddf)))
}



