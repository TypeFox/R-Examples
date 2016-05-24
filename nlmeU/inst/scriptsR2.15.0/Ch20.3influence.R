###################################################
### code chunk: Chap20.3influence_init
###################################################

options(width = 65, digits = 5, show.signif.stars = FALSE)
date()
packageVersion("nlmeU")
packageVersion("nlme")
packageVersion("lattice")

sessionInfo()
require(nlme)    
require(lattice) 


data(armd, package="nlmeU")

## Model M16.5 
lm3.form <-                 # (12.9, 16.17)
        formula(visual ~ visual0 + time + treat.f) 
fm16.5  <-                  # R16.13
lme(lm3.form,              
        random = list(subject = pdDiag(~time)),       
        weights =varPower(form=~time),
        data = armd)     



###################################################
### code chunk: R20.6a
###################################################
fm16.5ml <- update(fm16.5, method = "ML") # ML estimation
formula(fm16.5ml)                         # Recall model formula.
fm16.5ml$call$data                        # Recall data name.
logLik(fm16.5ml)                          # Log-likelihood value


###################################################
### code chunk: R20.6b
###################################################
beta0  <- fixef(fm16.5ml)                 # beta
names(beta0)                              # Long names
names(beta0) <- abbreviate(names(beta0), minlength = 7) # Short names 
beta0                                     # beta printed.
vcovb  <- vcov(fm16.5ml)                  # vcovb 
colnames(vcovb) <- names(beta0)           # Short names
vcovb                                     # vcovb printed. 


###################################################
### code chunk: R20.7a
###################################################
require(nlmeU)  
df1 <- subset(armd, subject %in% "1")          # Data for subject "1" 
logLik1(fm16.5ml, df1)                         # logLik_i for subject "1" 

lLik.i <- by(armd, armd$subject,
   FUN = function(dfi) logLik1(fm16.5ml, dfi))
lLik.i <- as.vector(lLik.i)   # Coerse array to vector  
lLik.i[1:5]                   # logLik_i for the first five subjects
sum(lLik.i)                   # Sum logLik_i; compare to Panel 20.6a


###################################################
### code chunk: R20.7b
###################################################
nx <- by(armd, armd$subject, nrow)             # ni
lLik.n <- lLik.i/as.vector(nx)                 # logLiki
outL <- lLik.n < -6                            # TRUE for values < -6
lLik.n[outL]                                   # logLiki/ni <  -6
subject.c <- levels(armd$subject)
subject.x <- as.numeric(subject.c)

plot(lLik.n ~ subject.x, type = "h")           # Fig. 20.1
points(subject.x[outL], lLik.n[outL], type = "p", pch = 16)
text(subject.x[outL], lLik.n[outL], subject.c[outL])


###################################################
### code chunk: R20.8a
###################################################
lmeU <- function(cx) { 
    dfU <- subset(armd, subject != cx)          # LOO data 
    update(fm16.5ml, data = dfU)                # LOO fit 
}

lmeUall        <- lapply(subject.c, lmeU)       # List with LOO fits
names(lmeUall) <- subject.c                     # Names assigned          


###################################################
### code chunk: R20.8b
###################################################
names(lmeUall)[1:6]                             
dataU6 <- lmeUall[["6"]]$data     # LOO data for subject 6
dim(dataU6)                       # Number of rows is 863
unique(dataU6$subject)[1:6]       # Subject no. 6 omitted


###################################################
### code chunk: R20.9a
###################################################
lLik <- function(cx){                  
    lmeU   <- lmeUall[[cx]]           # LOO fit extracted 
    lLikU  <- logLik(lmeU, REML = FALSE)  # LOO log-likelihood
    df.s   <-                         # Data for subject cx...
       subset(armd, subject == cx)                 
    lLik.s <- logLik1(lmeU, df.s)       # ...and log-likelihood.
    return(lLikU + lLik.s)            # "Displaced" log-likelihood...
}
lLikUall <- sapply(subject.c, lLik)     # ...for all subjects.         

dif.2Lik <- 2*(logLik(fm16.5ml) - lLikUall) # Vector of LDi
summary(dif.2Lik)


###################################################
### code chunk: R20.9b
###################################################
names(dif.2Lik) <- subject.c              # Subjects' ids assigned
outL  <-  dif.2Lik > 0.5                  # Outlying LDi's
dif.2Lik[outL]
library(lattice)
subject.f <- factor(subject.c, levels = subject.c)
myPanel <- function(x, y, ...){
   x1 <- as.numeric(x)
   panel.xyplot(x1, y, ...)
   ltext(x1[outL], y[outL], subject.c[outL])  # Label outlying LDi
}

dtp <-                                        # Fig. 20.2
   dotplot(dif.2Lik ~ subject.f, panel = myPanel, type = "h")           
lxlims <- length(dtp$x.limits)         
update(dtp, xlim = rep("", lxlims), grid = "h") 



###################################################
### code chunk: R20.10a
###################################################
betaUall <- sapply(lmeUall, fixef)          # Matrix with beta(-i)
vb.inv <- solve(vcovb)                       
CookDfun <- function(betaU){  
   dbetaU <- betaU - beta0                  # beta(-i) - beta
   CookD.value <- t(dbetaU) %*% vb.inv %*% dbetaU 
}
CookD.num <- apply(betaUall, 2, CookDfun)
(n.fixeff <- length(beta0))                 # Number of fixed effects
rankX <- n.fixeff                           # Rank of matrix X
CookD <- CookD.num/rankX                    # Cook's distance Di


###################################################
### code chunk: R20.10b
###################################################
outD <- CookD > 0.03                        # Outlying Di's
subject.c[outD]                             # Subjects' ids 

plot(CookD ~ subject.c, 
     ylab = "Cook's D", type = "h")         # Fig. 20.3
text(as.numeric(subject.c[outD]),
     CookD[outD], subject.c[outD])          # Annotation 
points(subject.c[outD], CookD[outD])

sessionInfo()
