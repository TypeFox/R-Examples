### Note: Simulations in Panel R20.11 take a long time

###################################################
### code chunk: Chap20.4init
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)
date()

sessionInfo()

require(nlme)    #

data(armd, package="nlmeU")

lm3.form   <-                           # (12.9)
     formula(visual ~ visual0 + time + treat.f) 
fm16.5  <-     #update(fm16.4,          # M16.5 <- M16.4           
     lme (lm3.form, 
     random = list(subject = pdDiag(~time)),
     weights = varPower(form = ~time), # 
     data = armd)            
fm16.5ml <- update(fm16.5, method = "ML")


###################################################
### code chunk: R20.11
###################################################
library(nlmeU)
simY <- simulateY(fm16.5ml, nsim = 1000, seed = 1238917)  # Simulated y from M16.1
auxDt <- subset(armd,                    # Auxiliary data
   select = c(subject, visual, visual0, time, treat.f))
simYsumm <- 
   apply(simY, 
         MARGIN = 2,                     # Over columns
         FUN = function(y){    
            auxDt$visual <- y            # Dependent variable updated
            auxFit <-                    # Update M16.1 with new y
               update(fm16.5ml, data = auxDt)
            summ <- summary(auxFit)      # Summary
            beta <- fixef(summ)     
            list(beta = beta)
})
simYsumm[[1]]                            # beta for the 1st simulation


###################################################
### code chunk: R20.12
###################################################
betaE <- sapply(simYsumm,            # Matrix with betas
  FUN = function(x) x$beta)
rowMeans(betaE)                      # Empirical beta (see Panel R20.6b)
cov(t(betaE))                        # Emoirical Var(beta)

## sessionInfo
sessionInfo()
