####
#
# ECOLOGICAL BIAS TESTING
#
###

library(general)
library(ipdmeta2)
source(file="~/code/R/ipdmeta2/beta/make.ecological.bias.ipd.R")
source(file="~/code/R/coxmeta/R/mixed.meta.data.R")
source(file="~/code/R/coxmeta/R/make.meta.from.ipd.R")

#FIXED EFFECTS MODEL WITH BIAS
beta = array(c(0,-.5,-.2,0,-.4,0))
names(beta) <- c("int","trt","x","x.bar","trt and x","trt and x.bar")

#WITHOUT BIAS
beta.nobias = array(c(0,-.5,-.2,-.2,-.4,-.4))
names(beta.nobias) <- c("int","trt","x","x.bar","trt and x","trt and x.bar")

ipd <- make.ecological.bias.ipd(rep(200,10),beta=beta,diag(0),.1,.2)
mixed <- mixed.meta.data(ipd$ipd,.5)
mixed$meta$x.bar <- mixed$meta$x

###FIXED EFFECTS

fit <- coxmeta.fixed(
     Surv(time,event)~trt*I(x-x.bar)+trt*x.bar,surv~log(time)+trt*x.bar,
     mixed$ipd,
     mixed$meta,
     mixed$meta$sigma2,
     mixed$meta$sub.group,
     beta.ad=runif(5),
     beta.ipd=beta[-1]
     )

fit$coef

#INTERACTION DIFFERNCE?
C <- array(c(0,0,0,-1,1,0,0))
C%*%fit$coef/sqrt(t(C)%*%fit$var%*%C)


###MIXED EFFECTS

fit <- coxmeta.mixed(
     Surv(time,event)~trt*I(x-x.bar)+trt*x.bar,surv~log(time)+trt*x.bar,
     ~(1|group),
    mixed$ipd,
    mixed$meta,
    ipd.groups=10,meta.groups=10,
    mixed$meta$sigma2,
    mixed$meta$sub.group,
    max.iter=5,
    min.sample=200,
    est.delta=.05,mc.step=1.5,df=10
    )

fit$coef				#MODEL FIT
sqrt(diag(fit$var$coef))		#STANDARD ERROR


sqrt(diag(fit$vcov))	#ESTIMATED FRAILTY STANDARD DEVIATION

C%*%fit$coef/sqrt(t(C)%*%fit$var$coef%*%C)

