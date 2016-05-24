###################################################
### code chunk Chap16merinit
###################################################

 
        
options(width=65, digits=5,show.signif.stars = FALSE) 
date()
packageVersion("nlmeU")
packageVersion("reshape")
packageVersion("lme4.0")
packageVersion("RLRsim")
sessionInfo()


data(armd, package = "nlmeU")


###################################################
### code chunk: R16.19a
###################################################
library(lme4.0)

fm16.1mer  <-                       
   lmer(visual ~ visual0 + time * treat.f + (1|subject),
        data = armd)
print(fm16.1mer, corr = FALSE)               


###################################################
### code chunk: R16.19b
###################################################
vcovb <- vcov(fm16.1mer)                     #
corb <- cov2cor(vcovb)                       #  
nms <- abbreviate(names(fixef(fm16.1mer)), 5)
rownames(corb) <- nms
corb


###################################################
### code chunk: R16.20a
###################################################
VarCorr(fm16.1mer)                    # Estimates of D, Corr(D)       
(sgma <- sigma(fm16.1mer))            # Estimate of sigma   


###################################################
### code chunk: R16.20b
###################################################
A <- getME(fm16.1mer, "A")            # A
I.n <- Diagonal(ncol(A))              # I_N
V <- sgma^2 * (I.n + crossprod(A))    # V = sigma^2 (I_N + A'A) 

str(getME(fm16.1mer, "flist"))        # Grouping factor

# V[3:6, 3:6]                         # V_i commented out (see R16.4)


###################################################
### code chunk: R16.21a
###################################################
coefs <- coef(summary(fm16.1mer)) 
ddf <- c(631, 231, 631, 231, 631)                 # Denominator df
pT <- 2 * (1 - pt(abs(coefs[, "t value"]), ddf))  # p -value
tTable <- cbind(coefs, ddf, pT)
printCoefmat(tTable, P.values = TRUE, has.Pvalue = TRUE)


###################################################
### code chunk: R16.21b
###################################################
(dtaov <- anova(fm16.1mer))
ddf1 <- ddf[-1]                          # ddf for intercept omitted
within(dtaov,
      {
       `Pr(>F)` <- pf(`F value`, Df, ddf1, lower.tail = FALSE)
       denDf <- ddf1
       })


###################################################
### code chunk: R16.22a
###################################################
SeedValue <- 17432157                 # Check it out
set.seed(SeedValue)

merObject <- fm16.1mer                      # M16.1 fit
simD1 <- simulate(merObject, nsim = 1000)   # Simulated y 
SimD1summ <- apply(simD1, 
     2,                               # Over columns
  function(y){
    auxFit <- refit(merObject, y)     # Refit M16.1 with new y
    summ <- summary(auxFit)           # Summary
    beta <- fixef(summ)               # beta
    Sx <- getME(auxFit, "theta")      # S element
    sgma <- sigma(auxFit)     
    list(beta = beta, ST = Sx, sigma = sgma)
             })



###################################################
### code chunk: R16.22b
###################################################
betaE  <-                          # Matrix with beta esimates
   sapply(SimD1summ, FUN = function(x) x$beta)
STe <- sapply(SimD1summ, FUN = function(x) x$ST)
sigmaE <- sapply(SimD1summ, FUN = function(x) x$sigma)


###################################################
### code chunk: R16.23a
###################################################
betaEm <- apply(betaE, 1, mean)         # Means (for each row)
betaEq <-                               # Quantiles
    apply(betaE, 1,                        
          FUN = function(x) quantile(x, c(0.5, 0.025, 0.975)))
ptE <-                                     # \zx!p?-values
    apply(betaE, 1,                       
          FUN = function(x){
             prb <- mean(x > 0)
             2 * pmax(0.5/ncol(betaE), pmin(prb, 1 - prb)) 
                           })
cbind(betaEm, t(betaEq), ptE)         # Bind results columnwise 


###################################################
### code chunk: R16.23b
###################################################
d11E <- STe * sigmaE
rndE <- rbind(d11E, sigmaE)                # Matrix with two rows
rndEm <- rowMeans(rndE)                    # Means (for each row)
rndEq <- apply(rndE, 1,                    # Quantiles
   FUN = function(x) quantile(x, c(0.5, 0.025, 0.975)))
cbind(rndEm, t(rndEq))                     # Bind results 


###################################################
### code chunk: R16.24
###################################################
names(sigmaE) <- names(STe) <- NULL          # For vectors
parSimD1  <-                                 # Matrix
   rbind(betaE, ST1 = STe, sigma = sigmaE) 
parSimD1t <-                                 # Transposed 
    data.frame(t(parSimD1), check.names=FALSE)
parSimD1s <-                                 # Subset
    subset(parSimD1t, select = -`(Intercept)`) # Intercept omitted
require(reshape)                             # melt function needed
densityplot(~value | variable,               # Fig. 16.13
            data = melt(parSimD1s),          # Molten data 
            scales = list(relation = "free"), 
            plot.points = FALSE) 
detach(package:reshape)    




###################################################
### code chunk: R16.25a
###################################################
lm2.form  <- visual ~ visual0 + time + treat.f + treat.f:time
vis.lm2 <- lm(lm2.form, data = armd)           # The null model
(RLRTstat <-                                   # Compare to R16.7
   -2 * as.numeric(logLik(vis.lm2, REML=TRUE) 
   - logLik(fm16.1mer)))        # log-REML for M16.1 (alternative)
0.5 * pchisq(RLRTstat, 1, lower.tail = FALSE)  # p-value


###################################################
### code chunk: R16.25b
###################################################
require(RLRsim)
exactRLRT(fm16.1mer) # \ref{M:sec:ARMDLMM:Rint} (alternative)


###################################################
### code chunk: R16.25c
###################################################
SeedValue <- 1321719  # check it out!!
set.seed(SeedValue)

lm2sim <- simulate(vis.lm2, nsim = 100)  # y simulated from the null model 

RLRTstatSim <- apply(lm2sim, 
  2,                                           # For each column
  function(y){ 
   dfAux  <- within(armd, visual <- y)         # Auxiliary data 
   lm0    <- lm(formula(vis.lm2), data = dfAux)# The null model
   llik0  <- as.numeric(logLik(lm0, REML=TRUE))# log-REML, the null
   llikA  <- as.numeric(logLik(refit(fm16.1mer, y))) 
   RLRTstat<- -2 * (llik0 - llikA)             # LR-test statistics              
             })
mean(RLRTstat <= RLRTstatSim)       # Empirical p-value


###################################################
### code chunk: R16.26a
###################################################
fm16.2mer  <-                        # M16.7
   lmer(visual ~ visual0 + time + treat.f + treat.f:time +   
           (1|subject) + (0 + time|subject),         
        data = armd)
summ <- summary(fm16.2mer)
coef(summ)                                     # t-Table
unlist(VarCorr(fm16.2mer))           # D. Short printout
sigma(fm16.2mer)                     # sigma 


###################################################
### code chunk: R16.26b
###################################################
fm16.2aux  <-                                  # Model M16.7 with ...
   update(fm16.2mer, . ~ . - treat.f:time)     #...interaction omitted
anova(fm16.2aux, fm16.2mer)



###################################################
### code chunk: R16.27a
###################################################
RML0 <- logLik(fm16.1mer)            # log-REML, M16.1 (null)
RMLa <- logLik(fm16.2mer)            # log-REML, M16.7 (alternative)
(RLRTstat <- -2 * as.numeric(RML0 - RMLa))  
.5 * pchisq(RLRTstat, 1, lower.tail = FALSE) +    # \zx{p}-value
  .5 * pchisq(RLRTstat, 2, lower.tail = FALSE)     


###################################################
### code chunk: R16.27b
###################################################
require(RLRsim)
mAux  <- lmer(visual ~               # Auxiliary model with ...
           visual0 + time + treat.f + treat.f:time +  
           (0 + time| subject),        # ... random slopes only. 
         data = armd)          
exactRLRT(m = mAux,                  # The auxiliary model
          m0= fm16.1mer,             # M16.1 (null)
          mA= fm16.2mer)             # M16.7 (alternative)

#### sessionInfo ###

sessionInfo()           # with packages attached

detach(package:lme4.0)
detach(package:RLRsim)

