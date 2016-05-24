###################################################
### code chunk Chap16lmeinit
###################################################
options(width=65, digits=5,show.signif.stars = FALSE) 
date()
packageVersion("nlmeU")
packageVersion("nlme")
packageVersion("RLRsim")
packageVersion("lattice")
  
sessionInfo()

library(nlme)
library(lattice)
data(armd, package = "nlmeU")
rlims <- c(-4.5, 4.5)
xlims <- c(0, 90)
xlimc <- c("4", "12", "24", "52wks")

lm1.form <- 
    formula(visual ~ -1 + visual0 + time.f + treat.f:time.f)
vis.gls1 <- gls(lm1.form, data = armd)
vis.gls1.tpwr <- 
    gls(lm1.form, weights = varPower(form=~time), data=armd)


###################################################
### code chunk: R16.1
###################################################
lm2.form <-                               # (16.1)
   formula(visual ~ visual0 + time + treat.f + treat.f:time)
(fm16.1 <-                                # M16.1
    lme(lm2.form,                         
    random = ~1|subject, data = armd))    # b_0i:(16.5)
printCoefmat(summary(fm16.1)$tTable,      # Print fixed effects, etc.
             has.Pvalue = TRUE, P.values = TRUE) # ... with p-values


###################################################
### code chunk: R16.2
###################################################
getGroupsFormula(fm16.1)           # Grouping formula
str(grpF <- getGroups(fm16.1))     # Grouping factor 
grpF[1:17]              
levels(grpF)[1:5]
range(xtabs(~grpF))      # Min, Max no. of observations


###################################################
### code chunk: R16.3a
###################################################
getVarCov(fm16.1, individual = "2")   # d_11: (16.6)
VarCorr(fm16.1)                       # d_11, sigma^2


###################################################
### code chunk: R16.3b
###################################################
getVarCov(fm16.1,                     
          type = "conditional",       # R_i: (16.6)
          individual = "2")   

###################################################
### code chunk: R16.4
###################################################
(fm16.1cov <- 
    getVarCov(fm16.1,                       
              type = "marginal",      # V_i: (16.7)
              individual = "2")) 
(cov2cor(fm16.1cov[[1]]))             # Corr(V_i)      


###################################################
### code chunk: R16.5
###################################################
(fm16.2 <-                            # M16.2 <- M16.1
   update(fm16.1,                        
          weights = varPower(form = ~ time),  # (9.4)
          data = armd))


###################################################
### code chunk: R16.6
###################################################
VarCorr(fm16.2)                            # d_11: (16.6), sigma^2
getVarCov(fm16.2,                          # R_i: (16.8)
          type = "conditional",
          individual = "2")
(fm16.2cov <-                              # V_i: (16.9)
   getVarCov(fm16.2,
             type = "marginal",
             individual = "2"))  
cov2cor(fm16.2cov[[1]])                    # Corr(V_i)


###################################################
### code chunk: R16.7a
###################################################
 plot(fm16.2)                         # Fig. 16.1


###################################################
### code chunk: R16.7b
###################################################
 plot(fm16.2,                         # Figure not shown
      resid(., type = "pearson") ~ time | treat.f,
      id = 0.05)
 bwplot(resid(fm16.2, type = "p") ~ time.f | treat.f,  # Fig. 16.2  
        # panel = panel.bwxplot2,     # User-defined panel (will be provided at a later time)
        data = armd)



###################################################
### code chunk: R16.7c
###################################################
qqnorm(fm16.2, ~resid(.) | time.f)   # Fig. 16.3
qqnorm(fm16.2, ~ranef(.))            # Fig. 16.4

###################################################
### code chunk: R16.8
###################################################

id <- 0.05                              # Argument for qnorm()
outliers.idx <- 
   within(armd,
          {
          resid.p <- resid(fm16.2, type = "pearson") # Pearson residuals
          idx <- abs(resid.p) > -qnorm(id/2)         # Indicator vector
          })
outliers  <- subset(outliers.idx, idx)  # Data with outliers
nrow(outliers)                          # Number of outliers
outliers$subject                        # IDs of outliers

###################################################
### code chunk: R16.9
###################################################
aug.Pred <-                  # augPred for M16.2
   augPred(fm16.2,                             
           primary = ~time,  # Primary covariate
           level = 0:1,      # Marginal(0) and subj.-spec.(1)
           length.out = 2)    
plot(aug.Pred, layout = c(4, 4, 1),   # Fig. 16.5
     key = list(lines = list(lty = c(1,2)),
                text = list(c("Marginal", "Subject-specific")),
                columns = 2))


###################################################
### code chunk: R16.10
###################################################
fm16.3 <-                                # M16.3 <- M16.2
   update(fm16.2,
          random = ~1 + time | subject,  # D: (16.16)
          data = armd)
getVarCov(fm16.3, individual = "2")     
intervals(fm16.3, which = "var-cov")     # 95% CI for theta_d, delta (16.8), sigma

###################################################
### code chunk: R16.11
###################################################
fm16.4 <-                               # M16.4 <- M16.3
    update(fm16.3,                              
           random = list(subject = pdDiag(~time)),  # Diagonal D
           data = armd)
intervals(fm16.4)                       # 95% CI for beta, theta_d, delta, sigma  


###################################################
### code chunk: R16.12
###################################################
anova(fm16.4, fm16.3)              # H0: d_12=0 

###################################################
### code chunk: R16.13
###################################################
lm3.form <- formula(visual ~ visual0 + time + treat.f) # (12.9)
fm16.5  <-  
   update(fm16.4,                        
          lm3.form, data = armd)            
summary(fm16.5)$tTable                    # beta, se(beta), t-test
intervals(fm16.5, which = "var-cov")      # 95% CI for theta_d, delta, sigma 


###################################################
### code chunk: R16.14
###################################################
VarCorr(fm16.5)                            # D, (16.16)sigma
getVarCov(fm16.5,                          # Ri (16.8)
          type = "conditional", individual = "2")
(fm16.5cov <-                              # Vi (16.9) 
   getVarCov(fm16.5,
             type = "marginal",
             individual = "2"))  
cov2cor(fm16.5cov[[1]])                    # Corr(Vi) 




###################################################
### code chunk: R16.15a
###################################################
(fm16.6 <-                                       #  M16.6 <- M16.3 
    update(fm16.3, weights = varIdent(form = ~1 | time.f)))


###################################################
### code chunk: R16.15b
###################################################
anova(fm16.3, fm16.6)        # varPower (M16.3) nested in varIdent (M16.6)

###################################################
### code chunk: R16.16
###################################################
AIC(fm16.1, fm16.2,                     # M16.1, M16.2
    fm16.3, fm16.4)                     # M16.3, M16.4
fm16.4ml <- update(fm16.4, method = "ML")
fm16.5ml <- update(fm16.5, method = "ML")
anova(fm16.4ml, fm16.5ml)               # M16.4 nested in M16.5

###################################################
### code chunk: R16.17a
###################################################
vis.gls1a   <-                                # Null model
   gls(lm2.form, data = armd)              
(anova.res  <- anova(vis.gls1a, fm16.1))      # Null vs. M16.1
(anova.res[["p-value"]][2])/2                 


###################################################
### code chunk: R16.17b
###################################################
library(RLRsim)    
exactRLRT(fm16.1)                             # M16.1 (alternative model)


###################################################
### code chunk: R16.18a
###################################################
fm16.7 <-                                  # M16.7 <- M16.4
   update(fm16.4, weights = NULL,          # Constant resid. variance 
          data = armd)              


###################################################
### code chunk: R16.18b
###################################################
(an.res <-                                 # M16.1 (null)
    anova(fm16.1, fm16.7))                 # M16.7 (alternative)
(RLRT <- an.res[["L.Ratio"]][2])           # LR-test statistic
.5 * pchisq(RLRT, 1, lower.tail = FALSE) + # 0.5* chi(1)^2 + 0.5 * chi(2)^2
  .5 * pchisq(RLRT, 2, lower.tail = FALSE) 


###################################################
### code chunk: R16.18c
###################################################
mAux <-        # Auxiliary model with random slopes only. 
   update(fm16.1, random = ~0 + time|subject, 
          data = armd)         
exactRLRT(m = mAux,           # Auxiliary model
          m0 = fm16.1,        # M16.1 (null)
          mA = fm16.7)        # M16.7 (alternative)


###################################################
### code chunk: R16.18d
###################################################
vis.lme2.sim <- 
     simulate(fm16.1, m2 = fm16.7, nsim = 10000, seed = 654321)
## save(vis.lme2.sim , file = "ARMDsimLMM.dat")
## load("ARMDsimLMM.dat")            # vis.lme2.sim loaded
plot(vis.lme2.sim, df = c(1, 2),     # Fig. 16.12
     abline = c(0, 1, lty = 2))
               
### sessionInfo

sessionInfo()                      # Before detaching packages
 
detach(package:nlme)
detach(package:RLRsim)
detach(package:lattice)


