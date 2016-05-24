TrialLevelIT <- function(Alpha.Vector, Mu_S.Vector=NULL, Beta.Vector, N.Trial,
                         Model="Reduced", Alpha=.05){
  
  # R2_ht          
   if (Model=="Full"){
      Model.1 <- glm(Beta.Vector ~ 1 + Mu_S.Vector + Alpha.Vector, family=gaussian)  
      L1 <- -2 * logLik(Model.1)[1]
    }
  
  if (Model==c("Reduced")){
      Model.1 <- glm(Beta.Vector ~ 1 + Alpha.Vector, family=gaussian)  
      L1 <- -2 * logLik(Model.1)[1]
    }
  
  Model.0 <- glm(Beta.Vector ~ 1, family=gaussian)      
  L0 <- -2 * logLik(Model.0)[1]
  g2 <- -(L1-L0)
  R2ht.value <- 1 - exp(-g2/N.Trial)
  k1 <- qchisq(Alpha, 1, g2)
  d1 <- qchisq((1-Alpha), 1, g2)
  R2ht.lb <- max(0, 1-exp(-k1/N.Trial))
  R2ht.ub <- min(1, 1-exp(-d1/N.Trial))
  R2ht <- data.frame(cbind(R2ht.value, R2ht.lb, R2ht.ub)) #output
  colnames(R2ht) <- c("R2ht", "CI lower limit", "CI upper limit")
  rownames(R2ht) <- c(" ")    
  
  fit <- 
    list(Alpha.Vector=Alpha.Vector, Beta.Vector=Beta.Vector, N.Trial=N.Trial, R2.ht=R2ht, 
         Call=match.call())
  
  class(fit) <- "TrialLevelIT"
  fit
  
}



