library(SAVE)

#############
# load data
#############

data(spotweldfield,package='SAVE')
data(spotweldmodel,package='SAVE')

##############
# create the SAVE object which describes the problem and
# compute the corresponding mle estimates
##############

sw <- SAVE(response.name="diameter", controllable.names=c("current", "load", "thickness"), 
			calibration.names=c("tuning"), field.data=spotweldfield, model.data=spotweldmodel,
			mean.formula=~1, bestguess=list(tuning=4.0))

# summary of the results

summary(sw)

##############
# obtain the posterior distribution of the unknown parameters 
##############
set.seed(0)
sw <- bayesfit(object=sw, prior=c(uniform("tuning", upper=8, lower=0.8)), n.iter=20000, n.burnin=100, n.thin=2)

# summary of the results
summary(sw)
# traceplots
plot(sw, option="trace", col=2)
# posterior of the calibration parameter
plot(sw, option="calibration", col=4, lty=3)
# posterior of the measurement error precision and
# bias function precision
plot(sw, option="precision", col=2, lty=4)

##############
# validate the computer model at chosen set of controllable
# inputs
###############

aload <- c(4.0,5.3)
curr <- seq(from=20,to=30,length=20)
g <- c(1,2)

xnew <- expand.grid(current = curr, load = load, thickness = g)

set.seed(0)
valsw <- validate(object=sw,newdesign=xnew,calibration.value='mean',n.burnin=100)

# summary of results
summary(valsw)
# plot results
plot(valsw)

# customized plot

av <- (valsw@validate)[,"pure.model"]
tau <- (valsw@validate)[,"tau.pm"] 

par(oma = c(1, 1, 2, 0), mfrow = c(2, 2))
for (i in 1:2) {
  for (j in 1:2) {
    v <- ((i - 1) * 40 + (j - 1) * 20 + 1):((i - 1) * 40 + j * 20)
    plot(curr, av[v], type = "l", ylim = c(3, 9), xlab = "current", ylab = "weld diameter")
    lines(curr, av[v] + tau[v], lty = 3)
    lines(curr, av[v] - tau[v], lty = 3)
    text(22, 9, paste("thickness=", g[i], ", load=", load[j], sep = ""), cex = 0.8, 
         pos = 1)
    # field data that correspond to this situation
    v <- ((i - 1) * 60 + (j - 1) * 30 + 1):((i - 1) * 60 + j * 30)
    data <- spotweldfield$diameter[v]
    inputs <- spotweldfield$current[v]
    points(inputs, data)
  }
}
mtext("Pure-model predictions", side = 3, outer = T, cex = 1.2)
par(mfrow = c(1, 1))

##########
# emulate the output of the model using predictcode
##########

# construct design at which to emulate the model
u <- 3.2
aload <- c(4.0,5.3)
curr <- seq(from=20,to=30,length=20)
g <- c(1,2)

xnewpure <- expand.grid(curr,aload,g)
xnewpure <- cbind(xnewpure,rep(u,dim(xnewpure)[1]))
names(xnewpure) <- c("current","load","thickness","tuning")
xnewpure <- as.data.frame(xnewpure)

set.seed(0)
pcsw<- predictcode(object=sw, newdesign=xnewpure, n.iter=20000, tol=1.E-12)

# Plot results
paths <- pcsw@samples
# compare estimates using simulation and the explicit formulae
avpure <- apply(paths,2,mean)
meanpure <- pcsw@modelmean

qts <- apply(paths,2,quantile,probs=c(0.025,0.975))
qtspure <- 1.96*sqrt(diag(pcsw@covmat))

par(oma=c(1,1,2,0),mfrow=c(2,2))
for(i in 1:2){
  for(j in 1:2){
    v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
    
    # explicit formulae
    plot(curr,meanpure[v],type="l",ylim=c(3,9),
    xlab="current",ylab="weld diameter",col=2)
    lines(curr,meanpure[v]+qtspure[v],lty=1,col=2)
    lines(curr,meanpure[v]-qtspure[v],lty=1,col=2)
    text(22,9,paste("gauge= ",g[i],", 
    aload=",aload[j],sep=""),cex=0.8,pos=1)
    # simulation-based
    lines(curr,avpure[v],type="l",col=1,lty=3)
    lines(curr,qts[1,v],lty=3,col=1)
    lines(curr,qts[2,v],lty=3,col=1)
    
  }
}
mtext("Emulate output of code at u = 3.20 - Spotweld Data",side=3,outer=T,cex=1.2)

#########
# bias-corrected prediction at a set of inputs
# using predictreality
##########

aload <- c(4.0,5.3)
curr <- seq(from=20,to=30,length=20)
g <- c(1,2)

xnew <- as.data.frame(expand.grid(curr,aload,g))
names(xnew)<-c("current","load","thickness")

# Obtain samples
set.seed(0)
prsw <- predictreality(object=sw, newdesign=xnew, tol=1.E-12)

# Plot results
# reality = model + bias
real <- prsw@modelpred+prsw@biaspred

# estimate
av <- apply(real,2,mean)

# tolerance bounds
tmpdata <- matrix(av,ncol=dim(real)[2],nrow=dim(real)[1],
                  byrow=T)
tmpdata <- real - tmpdata
tmpdata <- apply(tmpdata,2,abs)
tau.real <- apply(tmpdata,2,quantile,0.90)

# plot

par(oma=c(1,1,2,0),mfrow=c(2,2))
for(i in 1:2){
  for(j in 1:2){
    v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
    plot(curr,av[v],type="l",ylim=c(3,9),
    xlab="current",ylab="weld diameter")
    lines(curr,av[v]+tau.real[v],lty=3)
    lines(curr,av[v]-tau.real[v],lty=3)
    text(22,9,paste("gauge= ",g[i],", 
    aload=",aload[j],sep=""),cex=0.8,pos=1)
    # field data that correspond to this situation
    v <- ((i-1)*60+(j-1)*30+1):((i-1)*60+j*30)
    data <- spotweldfield$diameter[v]
    inputs <- spotweldfield$current[v]
    points(inputs,data)
    #dev.off()
  }
}
mtext("Bias-corrected predictions - Spotweld Data",side=3,
      outer=T,cex=1.2)

# tolerance bounds for the pure model predictions
# at u = 3.2

# Note: avpure was computed above

tmpdata <- matrix(avpure,ncol=dim(real)[2],nrow=dim(real)[1], byrow=T)
tmpdata <- real - tmpdata
tmpdata <- apply(tmpdata,2,abs)
tau.pure <- apply(tmpdata,2,quantile,0.90)

# plot

par(oma=c(1,1,2,0),mfrow=c(2,2))
for(i in 1:2){
	for(j in 1:2){
		v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
		plot(curr,avpure[v],type="l",ylim=c(3,9),
			 xlab="current",ylab="weld diameter")
        lines(curr,avpure[v]+tau.pure[v],lty=3)
        lines(curr,avpure[v]-tau.pure[v],lty=3)
		text(22,9,paste("gauge= ",g[i],", aload=",aload[j],sep=""),cex=0.8,pos=1)
# field data that correspond to this situation
		v <- ((i-1)*60+(j-1)*30+1):((i-1)*60+j*30)
		data <- spotweldfield$diameter[v]
		inputs <- spotweldfield$current[v]
		points(inputs,data)
#dev.off()
	}
}
mtext("Pure-model predictions - Spotweld Data",side=3,
outer=T,cex=1.2)

# plots using the output of validate()
# Note: valsw was computed above

avpure <- valsw@validate[,"pure.model"]
tau.pure <- valsw@validate[,"tau.pm"]

avbc <- valsw@validate[,"bias.corrected"]
taubc <- valsw@validate[,"tau.bc"]

avbias <- valsw@validate[,"bias"]
biasL <- valsw@validate[,"bias.Lower"]
biasU <- valsw@validate[,"bias.Upper"]

par(mfrow=c(3,4))
for(i in 1:2){
	for(j in 1:2){
		v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
		plot(curr,avpure[v],type="l",ylim=c(3.5,8.5),
			 xlab="current",ylab="weld diameter")
        lines(curr,avpure[v]+tau.pure[v],lty=3)
        lines(curr,avpure[v]-tau.pure[v],lty=3)
        title(main=paste("thickness=",g[i],", load=",aload[j],sep=""),cex=0.8)
		
		# field data that correspond to this situation
		v <- ((i-1)*60+(j-1)*30+1):((i-1)*60+j*30)
		data <- spotweldfield$diameter[v]
		inputs <- spotweldfield$current[v]
		points(inputs,data)
#dev.off()
	}
}

for(i in 1:2){
	for(j in 1:2){
		v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
		plot(curr,avbias[v],type="l",ylim=c(-2,2),
			 xlab="current",ylab="weld diameter")
        lines(curr,biasL[v],lty=3)
        lines(curr,biasU[v],lty=3)
	}
}

for(i in 1:2){
	for(j in 1:2){
		v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
		plot(curr,avbc[v],type="l",ylim=c(3.5,8.5),
			 xlab="current",ylab="weld diameter")
        lines(curr,avbc[v]+taubc[v],lty=3)
        lines(curr,avbc[v]-taubc[v],lty=3)
		
# field data that correspond to this situation
		v <- ((i-1)*60+(j-1)*30+1):((i-1)*60+j*30)
		data <- spotweldfield$diameter[v]
		inputs <- spotweldfield$current[v]
		points(inputs,data)
	}
}

#####
# another application - derivative wrt current
#####

# bias-corrected
aload <- c(4.0)
curr <- seq(from=20,to=30,length=80)
g <- c(1)

xnew <- expand.grid(curr,aload,g)
names(xnew)<-c("current","load","thickness")

# Obtain samples
set.seed(0)
prdersw <- predictreality(object=sw, newdesign=xnew, tol=1.E-12)

delta <- diff(curr)[1]
model <- prdersw@modelpred
dmodel <- diff(t(model))/delta
bias <- prdersw@biaspred
dbias <- diff(t(bias))/delta
dreal <- dmodel + dbias

# bias-corrected prediction
dav <- apply(dreal, 1, mean)

# tolerance bounds
tmpdata <- matrix(dav, ncol = dim(dreal)[2], nrow = dim(dreal)[1], byrow = F)
tmpdata <- dreal - tmpdata
tmpdata <- abs(tmpdata)
tau.real <- apply(tmpdata, 1, quantile, 0.9)

par(mfrow=c(1,1))
plot(curr[-1],dav,ty="l",ylim=c(min(dav-tau.real),max(dav+tau.real)),
     xlab="current",ylab="derivative")
lines(curr[-1],dav+tau.real,lty=2)
lines(curr[-1],dav-tau.real,lty=2)
title(main='Bias-corrected prediction of the derivative')

# pure-model prediction

# construct design at which to emulate the model
u <- 3.2
aload <- 4.0
curr <- seq(from=20,to=30,length=80)
g <- 1
xnewpure <- expand.grid(curr,aload,g,u)
names(xnewpure) <- c("current","load","thickness","tuning")
xnewpure <- as.data.frame(xnewpure)

set.seed(0)
pcdersw <- predictcode(object=sw, newdesign=xnewpure, n.iter=20000, tol=1.E-12)

samples <- pcdersw@samples
dersamples <- diff(t(samples))/diff(curr)[1]
dpure <- apply(dersamples,1,mean)
dpureup <- apply(dersamples,1,quantile,0.975)
dpurelow <- apply(dersamples,1,quantile,0.025)


plot(curr[-1],dpure,ty="l",ylim=c(min(dpurelow),max(dpureup)),
     xlab="current",ylab="derivative")
lines(curr[-1],dpureup,col=2)
lines(curr[-1],dpurelow,col=2)
title(main='Pure-mdoel prediction of the derivative')

# plot
#pdf(file="deriv.pdf")
plot(curr[-1],dav,ty="l",ylim=c(min(dav-tau.real), 
				max(dav+tau.real)),xlab="current",ylab="derivative")
lines(curr[-1],dav+tau.real,lty=2)
lines(curr[-1],dav-tau.real,lty=2)
lines(curr[-1],dpure,lty=4,lwd=2)
title(main='Bias-corrected and pure-morel prediction of the derivative')
#dev.off()
