### R code from vignette source 'survRM2-vignette3-1.Rnw'

###################################################
### code chunk number 1: preamble
###################################################
options(width=80, continue=" ")
makefig <- function(file, top=1, right=1, left=4) {
    pdf(file, width=11, height=7, pointsize=18)
    par(mar=c(4, left, top, right) +.1)
    }
library(survival)


###################################################
### code chunk number 2: survRM2-vignette3-1.Rnw:34-36 (eval = FALSE)
###################################################
## library(survival)
## ?pbc


###################################################
### code chunk number 3: survRM2-vignette3-1.Rnw:41-46
###################################################
library(survRM2)

D=rmst2.sample.data()
nrow(D)
head(D[,1:3])


###################################################
### code chunk number 4: survRM2-vignette3-1.Rnw:52-54
###################################################
plot(survfit(Surv(time, status)~arm, data=D), col=c("blue","red"), lwd=2, mark.time=F, xlab="Years",ylab="Probability")
legend("bottomleft", c("Placebo (arm=0)","D-penicillamine (arm=1)"), col=c("blue","red"), lwd=2)


###################################################
### code chunk number 5: survRM2-vignette3-1.Rnw:76-88
###################################################
  fit=survfit(Surv(D$time[D$arm==1], D$status[D$arm==1])~1)
  tau=10
  tmp.xx=c(0, fit$time); tmp.yy=c(1, fit$surv) ;
  idx=tmp.xx<=tau
  y.tau = min(tmp.yy[idx])
  xx=c(tmp.xx[idx],   tau)
  yy=c(tmp.yy[idx], y.tau)  
  x.step=sort(c(0, tmp.xx, tmp.xx))
  y.step=rev(sort(c(1,1,tmp.yy, tmp.yy[-length(tmp.yy)])))

  rmst=summary(fit, rmean=10)$table[5]



###################################################
### code chunk number 6: survRM2-vignette3-1.Rnw:91-108
###################################################
  par(mfrow=c(1,2))

  #--------
  plot(fit, mark.time=F, xlab="Years",ylab="Probability",conf.int=F, lwd=2, main="Restricted mean survival time (RMST)", col="red", cex.main=0.8)
  for (i in 1: (length(xx)-1)){  
  polygon(c(xx[i], xx[i+1], xx[i+1], xx[i]), c(0, 0, yy[i+1], yy[i]), col="pink", density=80, angle=80, lwd=2)
  }
  lines(x.step, y.step, col="red", lwd=3) 
  text(5,0.4, paste(round(rmst, digits=2),"years"), cex=0.9)
  
  #--------
  plot(fit, mark.time=F, xlab="Years",ylab="Probability", conf.int=F, lwd=2, main="Restricted mean time lost (RMTL)", col="red",cex.main=0.8)
  for (i in 1: (length(xx)-1)){  
  polygon(c(xx[i], xx[i+1], xx[i+1], xx[i]), c(yy[i], yy[i+1], 1,1), col="orange", density=80, angle=80, lwd=2)
  }
  lines(x.step, y.step, col="red", lwd=3) 
  text(7,0.8, paste(round(tau-rmst, digits=2),"years"), cex=0.9)


###################################################
### code chunk number 7: survRM2-vignette3-1.Rnw:126-129
###################################################
time=D$time
status=D$status
arm=D$arm


###################################################
### code chunk number 8: survRM2-vignette3-1.Rnw:132-133 (eval = FALSE)
###################################################
## rmst2(time, status, arm, tau=10)


###################################################
### code chunk number 9: survRM2-vignette3-1.Rnw:142-143 (eval = FALSE)
###################################################
## rmst2(time, status, arm)


###################################################
### code chunk number 10: survRM2-vignette3-1.Rnw:153-155
###################################################
obj=rmst2(time, status, arm, tau=10)
print(obj)


###################################################
### code chunk number 11: survRM2-vignette3-1.Rnw:163-164
###################################################
plot(obj, xlab="Years", ylab="Probability")


###################################################
### code chunk number 12: survRM2-vignette3-1.Rnw:180-181 (eval = FALSE)
###################################################
## rmst2(time, status, arm, tau=10, covariates=x)


###################################################
### code chunk number 13: survRM2-vignette3-1.Rnw:187-189
###################################################
x=D[,c(4,6,7)]
head(x)


###################################################
### code chunk number 14: survRM2-vignette3-1.Rnw:206-207
###################################################
rmst2(time, status, arm, tau=10, covariates=x)


