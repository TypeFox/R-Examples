###################################################
### code chunk: Chap18mer
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)
date()
packageVersion("nlmeU")
packageVersion("lme4.0")
packageVersion("lattice")

sessionInfo()
  
require(lattice) 

data(SIIdata, package = "nlmeU")


###################################################
### code chunk number 5: R18.17
###################################################
library(lme4.0)
fm18.6mer <-                                 
   lmer(mathgain ~ ses + minority + poly(mathkind, 3) + ses:minority + 
          (1|schoolid) + (1|schoolid:classid), # Syntax #1 (general)
        data = SIIdata, REML = FALSE)
summ  <- summary(fm18.6mer)
print(summ, corr = FALSE)

update(fm18.6mer, 
       mathgain ~ ses + minority + poly(mathkind, 3) + ses:minority + 
         (1|schoolid) + (1|classid))           # Syntax #2 

update(fm18.6mer, 
       mathgain ~ ses + minority + poly(mathkind, 3) + ses:minority + 
         (1|schoolid/classid))                 # Syntax #3


###################################################
### code chunk number 10: R18.18
###################################################
anova(fm18.6mer)	                 # Approximate F-test statistics
logLik(fm18.6mer)                  # ML value
unlist(VarCorr(fm18.6mer))         # d_1 and d_2
sigma(fm18.6mer)                   # sigma


###################################################
### code chunk number 11: R18.19a.
###################################################
rsd6 <- resid(fm18.6mer)
qqnorm(rsd6)            


###################################################
### code chunk number 12: R18.19b
###################################################
rnf6mer <- ranef(fm18.6mer)          # Random effects
rnf6qn  <- plot(rnf6mer, grid = TRUE)# Q-Q plot for random effects
update(rnf6qn[["schoolid:classid"]], # For classid  (see Fig. 18.9a)
       ylab = c("Random effects for classes")) 
update(rnf6qn[["schoolid"]],         # For schoolid (see Fig. 18.9b) 
       ylab = c("Random effects for schools")) 

#### sessionInfo with packages attached 
sessionInfo() 


