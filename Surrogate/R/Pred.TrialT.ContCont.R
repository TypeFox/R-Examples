Pred.TrialT.ContCont <- function(Object, mu_S0, alpha_0, alpha.CI=0.05){

  if (class(Object)=="BimixedContCont") {Dmat <- Object$D}
  if (class(Object)=="UnimixedContCont"){Dmat <- Object$D.Equiv}
  if (class(Object)=="BifixedContCont") {Dmat <- Object$D.Equiv}
  if (class(Object)=="UnifixedContCont"){Dmat <- Object$D.Equiv}
  
  if ((dim(Dmat)[1])==4){Model <- "Full"}
  if ((dim(Dmat)[1])==2){Model <- "Reduced"}
  
  # Full
  if (Model == "Full"){  
    
  if (class(Object)=="BimixedContCont"){
      Dmat <- Object$D
      mu_S_mod <- Object$Fixed.Effect.Pars[1,1]
      alpha_mod <- Object$Fixed.Effect.Pars[2,1]
      beta_mod <- Object$Fixed.Effect.Pars[4,1]
    }
    
    
  if (class(Object)=="UnimixedContCont"){
    Dmat <- Object$D.Equiv
    mu_S_mod <- Object$Fixed.Effect.Pars[1,1]
    alpha_mod <- Object$Fixed.Effect.Pars[2,1]
    beta_mod <- Object$Fixed.Effect.Pars[4,1]
  }
  
  
  if (class(Object)=="BifixedContCont"){
    Dmat <- Object$D.Equiv
    mu_S_mod <- mean(Object$Results.Stage.1[,3])
    alpha_mod <- mean(Object$Results.Stage.1[,5])
    beta_mod <- mean(Object$Results.Stage.1[,6])
    }
  
  if (class(Object)=="UnifixedContCont"){
    Dmat <- Object$D.Equiv
    mu_S_mod <-  mean(Object$Results.Stage.1[,3])
    alpha_mod <- mean(Object$Results.Stage.1[,5])
    beta_mod <- mean(Object$Results.Stage.1[,6])
    
  }
  
  p_1 <- (matrix(data = c(Dmat[4,1], Dmat[4,3]), nrow = 2, byrow = T))
  p_2 <- (matrix(data=c(Dmat[1,1], Dmat[3,1], Dmat[1,3], Dmat[3,3]), 
                 nrow=2, byrow=TRUE))
  p_3 <- matrix(data=c((mu_S0 - mu_S_mod), (alpha_0 - alpha_mod)), nrow = 2)
  exp_val <- beta_mod + ((t(p_1)) %*% (solve(p_2)) %*% p_3)  
  var_val <- Dmat[4,4] - (t(p_1) %*% solve(p_2) %*% p_3)
  
  } # end full
  
  if (Model == "Reduced"){  
    
    if (class(Object)=="BimixedContCont"){
      Dmat <- Object$D
      alpha_mod <- Object$Fixed.Effect.Pars[2,1]
      beta_mod <- Object$Fixed.Effect.Pars[4,1]
    }
    
    
    if (class(Object)=="UnimixedContCont"){
      Dmat <- Object$D.Equiv
      alpha_mod <- Object$Fixed.Effect.Pars[2,1]
      beta_mod <- Object$Fixed.Effect.Pars[4,1]
    }
    
    
    if (class(Object)=="BifixedContCont"){
      Dmat <- Object$D.Equiv
      alpha_mod <- mean(Object$Results.Stage.1[,3])
      beta_mod <- mean(Object$Results.Stage.1[,4])
    }
    
    if (class(Object)=="UnifixedContCont"){
      Dmat <- Object$D.Equiv
      alpha_mod <- mean(Object$Results.Stage.1[,3])
      beta_mod <- mean(Object$Results.Stage.1[,4])
    }
    
    exp_val <- beta_mod + ((Dmat[1,2] / Dmat[1,1]) * (alpha_0 - alpha_mod))  
    var_val <- Dmat[2,2] - ((Dmat[1,2]**2) / Dmat[1,1])
    
  } # end reduced
    
  Est.lb <- exp_val + (qnorm(alpha.CI/2) * sqrt(var_val))
  Est.ub <- exp_val - (qnorm(alpha.CI/2) * sqrt(var_val))
  
  fit <- 
    list(Beta_0=exp_val, Variance=var_val, Lower=Est.lb, Upper=Est.ub, 
         alpha.CI=alpha.CI, Surr.Model= Object, 
         alpha_0=alpha_0, Call=match.call())   
  
  class(fit) <- "PredTrialTContCont"
  fit
}



summary.PredTrialTContCont <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("Function call:\n")
  print(Object$Call)
  cat("\nResults:\n")
  cat("--------\n")
  cat("Expected treatment effect (variance) on T (beta_0) = ", round(Object$Beta_0, digits = 4), 
      " (", round(Object$Variance, digits = 4), ")\n", sep="")
  cat((1-Object$alpha.CI)*100, "% CI: [", round(Object$Lower, digits = 4), 
      "; ", round(Object$Upper, digits = 4), "] \n", sep="")
}



plot.PredTrialTContCont <- function(x, Size.New.Trial=5, CI.Segment=1, ...){
  
  plot(x = x$Surr.Model, Indiv.Level = FALSE, col="grey", Main.Trial=expression(paste("Predicted ",beta[0])), ...)
  points(x$alpha_0, x$Beta_0, col="black", pch=1, cex=Size.New.Trial, lwd=2)
  segments(x0 = x$alpha_0, y0 = x$Lower, x1 = x$alpha_0, y1 = x$Upper, 
           lwd=2, lty=2)
  segments(x0 = -CI.Segment+x$alpha_0, y0 = x$Lower, x1 = CI.Segment+x$alpha_0, y1 = x$Lower, lwd=CI.Segment, lty=2)
  segments(x0 = -CI.Segment+x$alpha_0, y0 = x$Upper, x1 = CI.Segment+x$alpha_0, y1 = x$Upper, lwd=CI.Segment, lty=2)
}