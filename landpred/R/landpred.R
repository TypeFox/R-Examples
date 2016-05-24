
###################################################
##calculate probability given no info ###########  
###################################################

Prob.Null <- function(t0, tau, data, weight=NULL, newdata = NULL) {
	Xi.long = data[,1]; Di.long = data[,2]; Xi.short = data[,3];  Di.short = data[,4];  Zi = data[,5]
	if(sum(Xi.long > t0) == 0) {stop("No long term events past the landmark time.")}
	if(is.null(weight))	{W2i <- Wi.FUN(data = cbind	(Xi.long,Di.long),t0=t0,tau=tau)}
	else{W2i=weight}
	Prob = (sum(1*W2i*(Xi.long <= t0 + tau)*(Xi.long > t0)))/(sum	(1*W2i*(Xi.long > t0)))
	data.column = vector(length = dim(data)[1])
	data.column[data[, 1] <= t0] = NA
	data.column[data[, 1] > t0] = Prob
	data = cbind(data, data.column)
    if(!is.null(newdata)) {
    	newdata.column = vector(length = dim(newdata)[1])
	    newdata.column[newdata[, 1] <= t0] = NA
	    newdata.column[newdata[, 1] > t0] = Prob
		newdata = cbind(newdata, newdata.column)   
    }
    return(list("Prob" = Prob, 
            "data" = data, "newdata" = newdata))

}

###################################################
##calculate IPCW##################################  
###################################################
Wi.FUN <- function(data, t0, tau, weight.given=NULL)	{
    Xi <- data[,1]; Di <- data[,2]; wihat <- rep(0, length(Xi))
    tmpind1 <- (Xi > t0)&(Xi <= (t0+tau)); tt0 <- c(Xi[tmpind1], t0 + tau); Ghat.tt0 <- Ghat.FUN(tt0,data, weight.given=weight.given)
    wihat[Xi > (t0+tau)] <- 1/Ghat.tt0[length(tt0)]
    wihat[tmpind1] <- Di[tmpind1]/Ghat.tt0[-length(tt0)]
    wihat
}


Ghat.FUN <- function(tt, data,type='fl', weight.given)	{
    tmpind <- rank(tt); Xi <- data[,1]; Di <- data[,2]
    summary(survfit(Surv(Xi,1-Di)~1,se.fit=F,type=type, weights=weight.given),sort(tt))$surv[tmpind]
}

###################################################
##calculate probability given covariate ###########  
###################################################

Prob.Covariate <- function(t0, tau, data, weight=NULL, short=TRUE, newdata = NULL) {
 	Xi.long = data[,1]
 	Di.long = data[,2] 
 	if(sum(Xi.long > t0) == 0) {stop("No long term events past the landmark time.")}
	if(short) { Zi = data[,5]}
	if(!short) { Zi = data[,3]}
    if (is.null(weight)) {
        W2i <- Wi.FUN(data = cbind(Xi.long, Di.long), t0 = t0, tau = tau)
    }
    else {
        W2i = weight
    }
    zi.cat = unique(Zi)
    covariate.results = matrix(nrow = length(zi.cat), ncol = 2)
    data.column = vector(length = dim(data)[1])
	data.column[data[, 1] <= t0] = NA
	if(!is.null(newdata)) {
		newdata.column = vector(length = dim(newdata)[1])
		newdata.column[newdata[, 1] <= t0] = NA
		newdata.Zi.column = dim(newdata)[2]
	}
    for (j in 1:length(zi.cat)) {
    	covariate.results[j, 1] = zi.cat[j]
        c = zi.cat[j]
        if(sum(Zi==c) < 10) {print(paste("Warning: Very few individuals with covariate value = ",c))}
        covariate.results[j, 2] = (sum(1 * (Zi == c) * W2i * (Xi.long <= t0 + tau) * (Xi.long > t0)))/(sum(1 * (Zi == c) * W2i * (Xi.long > t0)))
        data.column[data[, 1] > t0 & Zi == c] = covariate.results[j, 2]
        if(!is.null(newdata)) {
      		newdata.column[newdata[, 1] > t0 & newdata[,newdata.Zi.column] == c] = covariate.results[j, 2]
		}
        }
    data =cbind(data, data.column)
    if(!is.null(newdata)) {newdata = cbind(newdata, newdata.column)
    	newdata=as.data.frame(newdata)
		if(dim(newdata)[2] == 6) {names(newdata) = c("XL", "DL", "XS", "DS", "Z", "Probability") }
		if(dim(newdata)[2] == 4) {names(newdata) = c("XL", "DL", "Z", "Probability")}
	}
    data=as.data.frame(data)
	if(short) {names(data) = c("XL", "DL", "XS", "DS", "Z", "Probability") }
	if(!short) {names(data) = c("XL", "DL", "Z", "Probability")}
    return(list("Prob" = covariate.results, 
            "data" = data, "newdata" = newdata))
    }

###################################################
##calculate probability given covariate and TS ####  
###################################################

Prob.Covariate.ShortEvent <- function(t0, tau, data, weight=NULL, bandwidth = NULL, newdata=NULL) {
	data[,3] = log(data[,3])
	Xi.long = data[,1]
	Di.long = data[,2]
	Xi.short = data[,3]
	Di.short = data[,4]
	Zi = data[,5]
	if(sum(Xi.long > t0) == 0) {stop("No long term events past the landmark time.")}
	if(is.null(weight))	{W2i <- Wi.FUN(data = cbind(Xi.long,Di.long),t0=t0,tau=tau)}
	else{W2i=weight}
	if(is.null(bandwidth)) {
		h = bw.nrd(Xi.short[Xi.short<log(t0) & Xi.long > t0])
		bandwidth = h*sum(Xi.short<log(t0) & Xi.long > t0)^(-.10)
	}
	zi.cat = unique(Zi)
	data.column = vector(length = dim(data)[1])
	data.column[data[,1] <= t0] = NA
	for(j in 1:length(zi.cat)) {
		c = zi.cat[j]
		if(sum(Zi==c) < 10) {print(paste("Warning: Very few individuals with covariate value = ",c))}
		if(sum(data[,1] > t0 & data[,3]<log(t0) & data[,5] == c) < 50) {print("Warning: Smoothing over very few short term events")}
		if(sum(data[,1] > t0 & data[,3]>=log(t0) & data[,5] == c) < 10) {print(paste("Warning: Very few individuals with short term event past t0 and with covariate value = ",c))}
		short.ind = (data[,1] > t0 & data[,3]<log(t0) & data[,5] == c)
		Pr.2.t = Prob2.k.t(t=Xi.short[short.ind], t0=t0, tau=tau,data.use=data,bandwidth=bandwidth,covariate.value=c)
		data.column[short.ind] = Pr.2.t
		Pr.2 = Prob2(t0=t0, tau=tau,data=data,covariate.value = c)
		data.column[data[,1] > t0 & data[,3] >= log(t0) & data[,5] == c] = Pr.2
		}
	if(!is.null(newdata)) {
		if(is.null(names(newdata))) {
			newdata = as.data.frame(newdata)
		} 
		newdata.column = vector(length = dim(newdata)[1])
		newdata.column[newdata[,1] <= t0] = NA
		for(j in 1:length(zi.cat)) {
			c = zi.cat[j]
			short.ind = (newdata[,1] > t0 & newdata[,3]<t0 & newdata[,5] == c)
			Pr.2.t = Prob2.k.t(log(newdata[short.ind,3]), t0=t0, tau=tau,data.use=data,bandwidth=bandwidth,covariate.value=c)
			newdata.column[short.ind] = Pr.2.t
			Pr.2 = Prob2(t0=t0, tau=tau,data=data,covariate.value = c)
			newdata.column[newdata[,1] > t0 & newdata[,3] >= t0 & newdata[,5] == c] = Pr.2
			}
		newdata = cbind(newdata, newdata.column)
		names(newdata) = c("XL", "DL", "XS", "DS", "Z", "Probability")
	}
	data = cbind(data[,c(1:2)], exp(data[,3]), data[,c(4:5)], data.column)
	data=as.data.frame(data)
	names(data) = c("XL", "DL", "XS", "DS", "Z", "Probability")
	return(list("data" = data, "newdata" = newdata))			
}

prob2.single= function(K,W2i,Xi.long,tau,Di.short,Xi.short, Zi, t0,covariate.value) {
	P.2 = (sum(1*(Zi==covariate.value)*W2i*(Xi.long <= t0 + tau)*(Xi.long > t0)*(Di.short == 1)*(Xi.short<log(t0))*(Xi.short<log(t0))*K))/(sum(1*(Zi==covariate.value)*W2i*(Xi.long > t0)*(Di.short == 1)*(Xi.short<log(t0))*(Xi.short<log(t0))*K))
 		return(P.2)		
 		}

Prob2.k.t <- function(t,t0, tau, data.use,bandwidth, covariate.value, weight=NULL) {
	Xi.long = data.use[,1]
	Di.long = data.use[,2]
	Xi.short = data.use[,3]  
	Di.short = data.use[,4]  
	Zi = data.use[,5]
	if(is.null(weight))	{W2i <- Wi.FUN(data = cbind	(Xi.long,Di.long),t0=t0,tau=tau)}
	else{W2i=weight}
	K = Kern.FUN(Xi.short,t,bandwidth)
	P.2.return = apply(K, 1, prob2.single,W2i=W2i,Xi.long=Xi.long,tau=tau,Di.short =Di.short,Xi.short=Xi.short,Zi=Zi, t0=t0,covariate.value=covariate.value) 
	P.2.return[t>=log(t0)] = NA
	return(P.2.return)
}


Prob2 <- function(t0, tau, data, covariate.value, weight=NULL) {
	Xi.long = data[,1]; Di.long = data[,2]; Xi.short = data[,3];  Di.short = data[,4];  Zi = data[,5]
	if(is.null(weight))	{W2i <- Wi.FUN(data = cbind	(Xi.long,Di.long),t0=t0,tau=tau)}
	else{W2i=weight}	
	P.1 = (sum(1*(Zi==covariate.value)*W2i*(Xi.long <= t0 + tau)*(Xi.long > t0)*(Xi.short > log(t0))))/(sum	(1*(Zi==covariate.value)*W2i*(Xi.long > 	t0)*(Xi.short > log(t0))))
	return(P.1)
}

Kern.FUN <- function(zz,zi,bw) ## returns an (n x nz) matrix 
{
	out = (VTM(zz,length(zi))- zi)/bw
   	norm.k = dnorm(out)/bw
   	norm.k
}

VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }

mse.BW <- function(data, t0,tau,h, folds = 3,reps=2)
{	BW.vec = vector(length=folds)
	Xi.long = data[,1]
	Di.long = data[,2]
	Xi.short = data[,3]
	Di.short = data[,4]
	Zi = data[,5]
	W2i <- Wi.FUN(data = cbind	(Xi.long,Di.long),t0=t0,tau=tau)
	num.observ = dim(data)[1]
	nv = floor(num.observ/folds)
	replic = matrix(nrow = reps, ncol =folds)
	for(j in 1:reps) {
	for(k in 1:folds) {
		ind.val = sample(1:num.observ, nv)
		ind.tra = setdiff(1:num.observ,ind.val)
		d = Prob.Covariate.ShortEvent(t0, tau, data = data[ind.tra,], newdata = data[ind.val,], bandwidth = h, weight = W2i[ind.tra])
		prob.val = d$newdata
		BW = sum((((Xi.long[ind.val]<=t0+tau) - prob.val$Probability)^2)*(Xi.long[ind.val]>t0)*(Xi.short[ind.val]<t0)*(W2i[ind.val])*(Di.short[ind.val]==1),na.rm=T)
		BW.vec[k] = BW
	}
	replic[j,] = BW.vec
	}
return(mean(replic))
}

optimize.mse.BW = function(data, t0,tau,h.grid=seq(.01,2,length=50), folds=3, reps=2){
	opt = optimize(f = mse.BW, data=data,t0=t0,tau=tau, lower =min(h.grid), upper = max(h.grid), folds=folds, reps=reps)
	h=opt$minimum
	print("Consider undersmoothing to obtain optimal order.")
 return(list("h" = h))
}


BS.landmark <- function(t0, tau, data, short = TRUE, weight=NULL) {
	Xi.long = data[,1]
	Di.long = data[,2]
	if(short) {
		Xi.short = data[,3]
		Di.short = data[,4]
		Zi = data[,5]
		Prob.est = data[,6]
		}
	if(!short) {
		Zi = data[,3]
		Prob.est = data[,4]
		}
	if(sum(Xi.long > t0) == 0) {stop("No long term events past the landmark time.")}
	if(is.null(weight))	{W2i <- Wi.FUN(data = cbind	(Xi.long,Di.long),t0=t0,tau=tau)}
	else{W2i=weight}
	brier.score= sum( (W2i[Xi.long>t0])*((1*(Xi.long[Xi.long>t0] <= t0+tau) - Prob.est[Xi.long>t0])^2)  )/sum((1*(Xi.long>t0))*(W2i))
    return(list("Brier.score" = brier.score))

}

AUC.landmark <- function(t0, tau, data, short = TRUE, weight=NULL) {
	Xi.long = data[,1]
	Di.long = data[,2]
	if(short) {
		Xi.short = data[,3]
		Di.short = data[,4]
		Zi = data[,5]
		Prob.est = data[,6]
		}
	if(!short) {
		Zi = data[,3]
		Prob.est = data[,4]
		}
	if(sum(Xi.long > t0) == 0) {stop("No long term events past the landmark time.")}
	if(is.null(weight))	{W2i <- Wi.FUN(data = cbind	(Xi.long,Di.long),t0=t0,tau=tau)}
	else{W2i=weight}
	Si <- Prob.est; ss <- unique(sort(Si))
	AUC.est <- sum((helper.si(Si, "<=", Si, W2i*(Xi.long <= t0 + tau)*(Xi.long > t0))   + helper.si(Si, "<", Si, W2i*(Xi.long <= t0 + tau)*(Xi.long 	> t0)))*W2i*(Xi.long > t0 + tau)/2)/(sum(W2i*(Xi.long <= t0+tau)*(Xi.long > t0))*sum(W2i*(Xi.long > t0+tau)))  
    return(list("AUC.est" = AUC.est))

}


helper.si <- function(yy,FUN,Yi,Vi=NULL)  
  {  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
    if(is.null(Vi)){return(pos)}else{
      Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
      out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
      out[pos!=0,] <- Vi[pos,]
      if(is.null(dim(Vi))) out <- c(out)
      return(out) ## n.y x p
    }
  } 
  
cumsum2 <- function(mydat) 
  {
    if(is.null(dim(mydat))) return(cumsum(mydat))
    else{
      out <- matrix(cumsum(mydat), nrow=nrow(mydat))
      out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
      return(out)
    }
  }
