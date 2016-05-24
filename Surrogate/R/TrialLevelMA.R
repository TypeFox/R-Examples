TrialLevelMA <- function(Alpha.Vector, Beta.Vector, N.Vector, Weighted=TRUE,  
                             Alpha=.05){
    
  # stage 2
  if (Weighted==FALSE){
      Results.Stage.2 <- lm(Beta.Vector ~ Alpha.Vector)
  }
  if (Weighted==TRUE){
      Results.Stage.2 <- lm(Beta.Vector ~ Alpha.Vector, weights=N.Vector) 
    }
  
  
  # Trial R2
  Trial.R2.value <- as.numeric(summary(Results.Stage.2)[c("r.squared")])
  Trial.R2.sd <- sqrt((4*Trial.R2.value*(1-Trial.R2.value)^2)/(length(N.Vector)-3))
  Trial.R2.lb <- max(0, Trial.R2.value + qnorm(Alpha/2) *(Trial.R2.sd))
  Trial.R2.ub <- min(1, Trial.R2.value + qnorm(1-Alpha/2)*(Trial.R2.sd))
  Trial.R2 <- data.frame(cbind(Trial.R2.value, Trial.R2.sd, Trial.R2.lb, Trial.R2.ub))
  colnames(Trial.R2) <- c("R2 Trial", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ") 
  
  # Trial R
  Trial.R.value <- sqrt(as.numeric(summary(Results.Stage.2)[c("r.squared")]))
  Z <- .5*log((1+Trial.R.value)/(1-Trial.R.value)) 
  Trial.R.lb <- max(0, (exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(length(N.Vector)-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(length(N.Vector)-3)))))+1))
  Trial.R.ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(length(N.Vector)-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(length(N.Vector)-3)))))+1))
  Trial.R.sd <- sqrt((1-Trial.R.value**2)/(length(N.Vector)-2))
  Trial.R <- data.frame(cbind(Trial.R.value, Trial.R.sd , Trial.R.lb, Trial.R.ub))
  colnames(Trial.R) <- c("R Trial", "Standard Error", "CI lower limit", "CI upper limit")
  row.names(Trial.R) <- c(" ") 
  
  
  fit <- 
    list(Alpha.Vector=Alpha.Vector, Beta.Vector=Beta.Vector, N.Vector=N.Vector, Trial.R2=Trial.R2, 
         Trial.R=Trial.R, 
         Model.2.Fit=summary(Results.Stage.2), Call=match.call())
  
  class(fit) <- "TrialLevelMA"
  fit
  
}
