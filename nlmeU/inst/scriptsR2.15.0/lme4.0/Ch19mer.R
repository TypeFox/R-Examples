###################################################
### code chunk: Chap19mer_init: Panels R19.1 - R19.7.
### code for panels R19.8 - R19.13 is stored in ../Ch19.R
###################################################
options(digits = 5, show.signif.stars = FALSE)

packageVersion("nlmeU")
sessionInfo()
require(nlme)    
require(lattice)

###################################################
### code chunk: R19.1
###################################################
data(fcat, package = "nlmeU")
opts <- options()                         # Global options saved
options(contrasts =                       # Default contrasts changed 
          c("contr.sum", "contr.poly"))
options("contrasts")                      # Changes verified 

fm19.1 <-                                 # M19.1: (19.1)
   lm(scorec ~ id + target, data = fcat) 
formula(fm19.1)
length(beta <- coef(fm19.1))
beta[c(1:20, 535:547)]   
options(opts)                             # Global options restored


###################################################
### code chunk: R19.2a
###################################################
fxd <- coef(fm19.1)
idx <- substr(names(fxd), 1, 2) == "id"          # Logical vector
names(fxi <- fxd[idx])
(fxd.id <- c(fxi, "id539" = -sum(fxi)))          # beta_2,s

###################################################
### code chunk: R19.2b
###################################################
idx <- substr(names(fxd), 1, 6) == "target"
names(fxi <- fxd[idx])
(fxd.trgt <- c(fxi, "target9" = -sum(fxi)))       # beta_1,t


###################################################
### code chunk: R19.3
###################################################
library(lme4.0)
system.time(
fm19.2mer <- lmer(scorec ~ (1|target) + (1|id), 
    data=fcat)
)
fm19.2mer                                           # M19.2: (19.2)


###################################################
### code chunk: R19.4
###################################################
summ.merFit <- summary(fm19.2mer)          # Summary of the model-fit
isREML(fm19.2mer)                          # REML used?
(cl <- getCall(summ.merFit))               # Function call
cl$data                                    # The name of data frame
formula(fm19.2mer)                         # Formula
fixef(fm19.2mer)                           # beta
coef(summ.merFit)                          # 
(VCorr <- unlist(VarCorr(fm19.2mer)))      # d_S, d_T
sigma(fm19.2mer)                           # sigma


###################################################
### code chunk: R19.5
###################################################
rnf <- ranef(fm19.2mer)         # ranef.mer-class object
names(rnf)
length(plx <- plot(rnf))        # Two Q-Q plots saved. 
plx[1]                          # Fig. 19.1a 
plx[2]                          # Fig. 19.1b  
plot(coef(fm19.2mer))           # Fig. 19.2


###################################################
### code chunk: R19.6a
###################################################
dpx <- dotplot(rnf)
# dpx[1]                      # Dotplot for id (not shown)
dpx[2]                        # Fig. 19.3a


###################################################
### code chunk: R19.6b
###################################################
rnf.pVar <- ranef(fm19.2mer, postVar = TRUE) # ranef.mer-class object
dpx.pVar <- dotplot(rnf.pVar)
# dpx.pVar[1]                 # Dotplot for id (not shown)
dpx.pVar[2]                   # Fig. 19.3b




###################################################
### code chunk: R19.7a
###################################################
(eVCorr <-  sapply(rnf, var))
VCorr                           
all(eVCorr < VCorr)


###################################################
### code chunk: R19.7b
###################################################
rnf.id    <- rnf$id               # Data frame with b_1,S
arnf.id   <- abs(rnf.id)          # abs(b_1,S)
afxd.id   <- abs(fxd.id)          # abs(beta_1,S)    
range(afxd.id - arnf.id)


###################################################
### code chunk: R19.7c
###################################################
rnf.trgt  <- rnf$target
arnf.trgt <- abs(rnf.trgt) 
afxd.trgt <- abs(fxd.trgt)
range(afxd.trgt - arnf.trgt)


###################################################
### code chunk: R19.7d
###################################################
names(dt  <- data.frame(afxd.id, arnf.id))
names(dt)[2] <- "arnf.id"

myPanel <- function(x, y, ...){
  panel.grid(h = -1, v = -1)
  panel.lines(c(0, 3), c(0, 3), lty = 2)
  panel.xyplot(x, y, ...)
}

xyplot(arnf.id ~ afxd.id,         # Fig. 19.4a
       data = dt, panel = myPanel)

sessionInfo()


