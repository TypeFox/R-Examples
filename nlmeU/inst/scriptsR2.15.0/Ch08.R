
###################################################
### code chunk: Chap8init
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)
date()
packageVersion("nlmeU")
packageVersion("nlme")
sessionInfo()
data(armd, package = "nlmeU")

## lm1.form was defined in Chapter 4
lm1.form <- formula(visual ~ -1 + visual0 + time.f + treat.f:time.f )
library(nlme)

###################################################
### code chunk: R8.1
###################################################
(val <- c("12wks" = 0.5, "24wks" = 2))  # delta1 = 1, delta2 = 0.5, delta3 = 2     
(fix <- c("52wks" = 3))                 # delta4 = 3 (fixed)
frm  <- formula(~1|time.f)              # time.f is a stratifying factor
(vf0 <- 
   varIdent(value = val,                # Var. function object defined... 
            fixed = fix,
            form  = frm)) 
(vf0i <- Initialize(vf0, armd))         # ... and initialized


###################################################
### code chunk: R8.2a
###################################################
coef(vf0i, unconstrained = FALSE, allCoef = TRUE) # All delta coefs
coef(vf0i, unconstrained = FALSE, allCoef = FALSE)# Varying only


###################################################
### code chunk: R8.2b
###################################################
coef(vf0i, unconstrained = TRUE, allCoef = TRUE)  # All delta* coefs
coef(vf0i, unconstrained = TRUE, allCoef = FALSE) # Varying (default)
coef(vf0i) <- c(-0.6, 0.7)                        # New coefs assigned   
coef(vf0i, allCoef = TRUE)                        # All coefs printed


###################################################
### code chunk: R8.3
###################################################
summary(vf0i)               # Summary
formula(vf0i)               # Variance function formula
getCovariate(vf0i)          # Variance covariate
getGroupsFormula(vf0i)      # Formula for variance strata
length(stratum <-           # Length of stratum indicator
         getGroups(vf0i)) 
unique(stratum)             # Unique strata
stratum[1:6]                # First six observations
varWeights(vf0i)[3:6]       # Variance weights 1/lambdai:(7.8)
logLik(vf0i)                # Contribution to the log-likelihood

###### sessionInfo() with packages attached
sessionInfo()
detach(package:nlme)
