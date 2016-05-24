Prentice <- function(Dataset, Surr, True, Treat, Pat.ID, Alpha=.05){
  
  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Trial.ID <- rep(1, nrow(Dataset))
  Pat.ID <- Dataset[,paste(substitute(Pat.ID))]
  
  Data.Proc <- .Data.Processing(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat, Trial.ID=Trial.ID, Pat.ID=Pat.ID, Min.Trial.Size=0)
  wide <- Data.Proc$wide
  dataS <- Data.Proc$dataS
  dataT <- Data.Proc$dataT
  Data.analyze <- Data.Proc$Data.analyze
  N.total <- Data.Proc$N.total
  
  # Prentice
  model12 <- lm(cbind(wide$Surr, wide$True)~ wide$Treat, data=wide) 
  alpha <- model12$coefficients[2,1]
  beta <- model12$coefficients[2,2]
  model3 <- lm(wide$True ~ wide$Surr, data=wide)
  gamma <- model3$coefficients[2]
  model4 <- lm(wide$True ~ wide$Treat + wide$Surr, data=wide)
  beta_s <-  model4$coefficients[2]
  
  P_model1 <- summary(model12)$"Response wide$Surr"$coefficients
  rownames(P_model1) <- c("Intercept", "Treatment")
  P_model2 <- summary(model12)$"Response wide$True"$coefficients
  rownames(P_model2) <- c("Intercept", "Treatment")
  P_model3 <- summary(model3)$coefficients
  rownames(P_model3) <- c("Intercept", "Surrogate")
  P_model4 <- summary(model4)$coefficients
  rownames(P_model4) <- c("Intercept", "Treatment", "Surrogate")
  
  P_crit1 <- summary(model12)$"Response wide$Surr"$coefficients[2, 4] 
  P_crit2 <- summary(model12)$"Response wide$True"$coefficients[2, 4] 
  P_crit3 <- summary(model3)$coefficients[2,4] 
  P_crit4 <- data.frame(summary(model4)$coefficients)[2, 4]  
  
  if ((P_crit1 < Alpha & P_crit2 < Alpha & P_crit3 < Alpha & P_crit4 > Alpha)==TRUE) {Prentice.Passed <- TRUE}
  if ((P_crit1 < Alpha & P_crit2 < Alpha & P_crit3 < Alpha & P_crit4 > Alpha)==FALSE) {Prentice.Passed <- FALSE}
  
  fit <-   
    list(Prentice.Model.1=P_model1, Prentice.Model.2=P_model2, Prentice.Model.3=P_model3, Prentice.Model.4=P_model4, 
         Prentice.Passed=Prentice.Passed, Call=match.call())  
  
  class(fit) <- "Prentice"
  fit
  
}


summary.Prentice <- function(object, ..., Object){
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Estimated parameters models 1-4")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("\nalpha = ", round(x = Object$Prentice.Model.1[2,1], digits=4), ",  p-value = ", Object$Prentice.Model.1[2,4], sep="")
  cat("\nbeta = ", round(x = Object$Prentice.Model.2[2,1], digits=4), ",  p-value = ", Object$Prentice.Model.2[2,4], sep="")
  cat("\ngamma = ", round(x = Object$Prentice.Model.3[2,1], digits=4), ",  p-value = ", Object$Prentice.Model.3[2,4], sep="")
  cat("\nbeta_S = ", round(x = Object$Prentice.Model.4[2,1], digits=4), ",  p-value = ", Object$Prentice.Model.4[2,4], sep="")

  cat("\n\nPrentice criteria passed: ")
  if (Object$Prentice.Passed == FALSE){cat("No\n\n")}
  if (Object$Prentice.Passed == TRUE){cat("Yes\n\n")}
  
}
