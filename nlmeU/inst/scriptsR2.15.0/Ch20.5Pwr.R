### Code for Chapter 20. Section 5

###################################################
### code chunk:  Chap.20.5init
###################################################

options(width=65, digits=5,show.signif.stars = FALSE)   
date()
packageVersion("nlmeU")
packageVersion("nlme")
sessionInfo()


##  M16.5 (best in Ch.16)
require(nlme)
data(armd, package = "nlmeU")
lm3.form <-   
    formula(visual ~ visual0 + time + treat.f) 

fm16.5 <- lme(lm3.form,    #  See R16.13         
        random = list(subject = pdDiag(~time)),       
        weights = varPower(form = ~time),
        data = armd)     


###################################################
### code chunk: R20.13a
###################################################
formula(fm16.5)                           # Recall formula
fixef(fm16.5)                             # beta


###################################################
### code chunk: R20.13b
###################################################
anova(fm16.5)                             # Default call 
anova(fm16.5, Terms = "treat.f")          # Terms argument
anova(fm16.5, L = c("treat.fActive" = 1)) # L argument


###################################################
### code chunk: R20.14a 
###################################################
alpha <- 0.05                              # alpha
df1 <- 1                                   # numDF
df2 <- 231                                 # denDF
Fvalue <- 5.5345                           # F-value (from R20.13b)
(Fcrit <-               # Critical value for the F-test under H_0
   qf(1 - alpha,     
      df1 = df1, df2 = df2, ncp =0))
nc <- Fvalue * df1                         # Non-centrality parameter
pf(Fcrit,                                  # Power
   df1 = df1, df2 = df2, 
   ncp = nc, lower.tail = FALSE)


###################################################
### code chunk: R20.14b
###################################################
library(nlmeU)
Pwr(fm16.5)                                # Default call 
Pwr(fm16.5,  L = c("treat.fActive" = 1))   # The L argument


###################################################
### code chunk:  cleanup
###################################################
rm(Fcrit, nc, df1, df2, Fvalue)



###################################################
### code chunk: R20.15a
###################################################
npg <- 20                                  # No of subjects per group
subject <- 1:(2*npg)                       # Subjects' ids
treat.f <- gl(2, npg, labels = c("Placebo", "Active"))
dts <- data.frame(subject, treat.f)        # Subject-level data

dtL <- 
   list(time = c(4, 12, 24, 52),
        subject = subject)
dtLong <- expand.grid(dtL)                 # "Long" format
mrgDt  <- merge(dtLong, dts, sort = FALSE) # Merged
exmpDt <- 
   within(mrgDt, 
          {    
           m0 <- 65 - 0.1 * time    # Under H0  
           mA <- 75 - 0.1 * time    # Under Ha. 85 changed to 75
           mA <- ifelse(treat.f %in% "Active", mA, m0) 
          })

###################################################
### code chunk: R20.15b
###################################################
selDt <- 
   with(exmpDt, 
        {
         lvls <- levels(treat.f)           # "Placebo", "Active"
         i <- match(lvls, treat.f)         # 1, 81
         subj <- subject[i]                # 1, 21
         subset(exmpDt, subject %in% subj) 
        })        
library(lattice)
xyplot(mA ~ time,                          # Fig. 20.5
       groups = treat.f,           
       data = selDt, 
       type = "l",
       auto.key = list(lines = TRUE, points = FALSE),
       grid = TRUE)


###################################################
### code chunk: R20.16a
###################################################
D0 <- diag(c(100, 0.09))                 # calligraphic D     
sgma  <- 5                               # sigma
(D  <- D0/(sgma*sgma))                   # D
(pd1 <- pdDiag(D, form = ~time, data = armd))
(vF <- varPower(form = ~time, fixed = 0.15))


###################################################
### code chunk: R20.16b
###################################################
cntrl <- 
   lmeControl(maxIter = 0, msMaxIter = 0, niterEM = 0, 
              returnObject = TRUE, opt = "optim")

fmA <- 
   lme(mA ~ time + treat.f,
       random = list(subject = pd1),
       weights = vF,
       data = exmpDt,
       control = cntrl)
fixef(fmA)                               # beta verified
sigma(fmA)                               # sigma approx. 0


###################################################
### code chunk: R20.17a
###################################################
Pwr(fmA, sigma = sgma, L = c("treat.fActive" = 1))


###################################################
### code chunk: R20.17b
###################################################
dif <- seq(1, 15, by = 0.1)                 # Delta
dim(dif) <- c(length(dif), 1)
colnames(dif) <- "treat.fActive"
dtF <-                                      # Data for Fig.20.6 
   Pwr(fmA, sigma = sgma,                  
       L = c("treat.fActive" = 1), altB = dif)   
dtF[ ,1:4]                                  # Four variables

###################################################
### code chunk: R20.17c
###################################################
xyplot(Power ~ treat.fActive,                # Fig.20.6
       data = dtF, type="l",
       auto.key = list(lines = TRUE, points = FALSE),
       grid = TRUE)




###################################################
### code chunk: R20.18
###################################################
simA <- simulateY(fmA, sigma = sgma, nsim = 1000, seed = 8917437) # Simulation
dt <- exmpDt                                      # Working copy 
simfmA <- 
   apply(simA,     
         2,                                       # Over columns
         function(y){
            dt$mA <- y                            # mA over-written
            auxFit <- update(fmA, data = dt)
            anova(auxFit)                         # ANOVA table
            })
simfmA[[1]]                                       # First ANOVA 

###################################################
### code chunk: R20.19
###################################################
FstatE <-                     # Empirical F-test statistics under HA 
   sapply(simfmA, function(x) x["treat.f", "F-value"])
summary(FstatE)
Fcrit <- qf(1- 0.05, 1, 38, ncp =0)
(nsim <- length(FstatE))
(powerE <- sum(FstatE > Fcrit)/nsim)             # Empirical power

### sessionInfo

sessionInfo()                                    # Before detaching packages

detach(package:nlme)
