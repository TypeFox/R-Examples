date()
## See R20.15a
npg <- 20 # No of subjects per group
subject <- 1:(2*npg) # Subjects' ids
treat.f <- gl(2, npg, labels = c("Placebo", "Active"))
dts <- data.frame(subject, treat.f) # Subject-level data
dtL <-
  list(time = c(4, 12, 24, 52),
  subject = subject)
  dtLong <- expand.grid(dtL) # "Long" format
  mrgDt <- merge(dtLong, dts, sort = FALSE) # Merged
  exmpDt <-
  within(mrgDt,
  {
  m0 <- 75 - 0.1 * time # Under H0                      # Changed from 65 to 75 (Feb. 2013)
  mA <- 85 - 0.1 * time # Under HA
  mA <- ifelse(treat.f %in% "Active", mA, m0)
  })

## See R20.16a
data(armd, package = "nlmeU")
library(nlme)
D0 <- diag(c(100, 0.09)) 
sgma <- 5 
(D <- D0/(sgma*sgma)) # D
(pd1 <- pdDiag(D, form = ~time, data = armd))
(vF <- varPower(form = ~time, fixed = 0.15))

## See R20.16b
cntrl <-
  lmeControl(maxIter = 0, msMaxIter = 0, niterEM = 0,
             returnObject = TRUE, opt = "optim")
fmA <-
   lme(mA ~ time + treat.f,
       random = list(subject = pd1),
       weights = vF,
       data = exmpDt,
       control = cntrl)
fixef(fmA)

detach(package:nlme)      

library(nlmeU)
sigma(fmA)
Pwr(fmA, sigma = sgma, L = c("treat.fActive" = 1)) 

dif <- 10
dim(dif) <- c(length(dif), 1)
colnames(dif) <- "treat.fActive"
dtF <- Pwr(fmA, sigma = sgma,
    L = c("treat.fActive" = 1), altB = dif)
dtF


packageVersion("nlme")
sessionInfo()
detach(package:nlmeU)



