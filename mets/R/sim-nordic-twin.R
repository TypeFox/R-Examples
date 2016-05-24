
F1addfg<-function(t,lam0=0.13,beta=c(-0.5),x=0) # FG
{ ## {{{
baset <- lam0*pnorm((t-.70)/0.15)
return( 1 - exp(-baset*exp(c(x * beta)))) 
} ## }}}

##' @export 
corsim.prostate <- function(n,theta=1,thetaslope=0,censS=c(0,1),pcens=0.5,test=0,mt=1,same.cens=TRUE,country=TRUE,
   delayed=FALSE,ptrunc=0.5,lam0=0.13,truncS=c(0,1)) 
{ ## {{{
###n <- 10; theta <- 1; thetaslope <- 0; mt <- 1
if (country==TRUE) xl <- sample(1:4,n,replace=TRUE) else xl <- rep(1,n)
x <- (xl==1)
tt<-seq(0,1,length=100)
###
###n=100;theta=1;lam0=0.5;beta=0.3;crate=2
thetat <- exp(log(theta))
F11x<-F1addfg(1,x=x,lam0=lam0)
F12x<-F1addfg(1,x=x,lam0=lam0)
###
thetaslut <- exp(log(theta)+thetaslope*(1-1/2))
p11 <- thetaslut*F11x*F12x/((1-F11x)+thetaslut*F11x)
p12 <- F11x-p11
p21 <- F12x-p11
p22 <- 1-F12x-F11x+p11
###print(apply(cbind(p11,p12,p21),1,sum))
###print(cbind(p11,p12,p21,p22))
if (test==1) { ## {{{
for (i in 1:2) {
print(x[i,]); 
F11xt<-F1addfg(tt,x=x[i,])
F12xt<-F1addfg(tt,x=x[i,])
p11t <- thetat* F11xt*F12xt/((1-F11xt)+thetat*F11xt)
cortt <- ((p11t)/(F12xt-p11t))/(F11xt/(1-F11xt))
###plot(tt,log(cortt))
if (i==1) { 
plot(tt,p11t,type="l",ylim=c(0,0.1),xlim=c(0,mt))
###lines(tt,F11x[i]-p11t,col=2)
###lines(tt,F12x[i]-p11t,col=2)
} else lines(tt,p11t,col=2);
###if (sum(diff(p11t<0))>0) stop("dec\n"); 
###p11 <- max(p11t)
###p12 <- F11x[i]-p11
###p21 <- F12x[i]-p11
###p22 <- 1- F12x[i]-F11x[i]+p11
###pnn <- 1- F12x[i]-F11x[i]+p11
}
} ## }}}
###apply(cbind(p11,p12,p21,p22),1,sum)
###
print(table(F11x))
print(table(p11))

types <- rep(0,n)
causes <- matrix(0,n,2)
stime<-matrix(1+1,n,2); 
for (i in 1:n)
{ ## {{{ 
ptype <- runif(1)
if (ptype<=p11[i]) {
	types[i] <- 1
	myhazx<-F1addfg(tt,x=x[i,])/F12x[i]
###	if (abs(max(myhazx)-1)> 0.001) stop("not dist\n"); 
	stime[i,2]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
	f1<- F1addfg(tt,x=x[i,])
       	myhazx<- (F12x[i]/p11[i]) * (thetat*f1/((1-f1)+thetat*f1))
###	if (abs(max(myhazx)-1)> 0.001) stop("not dist\n"); 
	stime[i,1]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
	causes[i,] <- c(1,1)
}
if ((ptype>p11[i]) & (ptype<=p12[i]+p11[i])) {
	types[i] <- 2
	f1 <- F1addfg(tt,x=x[i,])
	myhazx<- ( f1 - thetat*F12x[i]*f1/((1-f1)+thetat*f1))/p12[i]; 
	myhazx <- f1/F11x[i]
###	if (abs(max(myhazx)-1)> 0.001) stop("not dist 2 \n"); 
	stime[i,1]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
	causes[i,] <- c(1,2)
	stime[i,2] <- runif(1)*1
}
if ((ptype>p11[i]+p12[i]) && (ptype<=p21[i]+p12[i]+p11[i])) {
	types[i] <- 3
        f2 <- F1addfg(tt,x=x[i,])
	myhazx <-  (f2 - (thetat*F11x[i]*f2/((1-F11x[i])+thetat*F11x[i])))/p21[i]; 
	myhazx <- f2/F12x[i]
###	if (abs(max(myhazx)-1)> 0.001) stop("not dist3 \n"); 
	stime[i,2]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
	causes[i,] <- c(2,1)
	stime[i,1] <- runif(1)*1
}
if (ptype>p11[i]+p12[i]+p21[i] ) {
	types[i] <- 4
	causes[i,] <- c(2,2)
	stime[i,1:2] <- runif(2)*1
}

} ## }}} 
stime <- c(stime)
cause <- c(causes)

###print(summary(stime))
###print(sum(stime==mt))

###same.cens=TRUE
if (same.cens==TRUE) {
	ctime <- rep(rbinom(n,1,pcens)*runif(n,censS),each=2)
        ctime[ctime==0] <- 1;
}
else {
	ctime<- rbinom(2*n,1,pcens)*runif(2*n,censS)
        ctime[ctime==0] <- 1;
}

cens <- (ctime< stime)
time <- ifelse(cens,ctime,stime)
cause <- ifelse(cens,0,cause)
id <- rep(1:n,rep(2,n))

country <- c()
country[xl==1] <- "SWE"
country[xl==2] <- "DK"
country[xl==3] <- "FIN"
country[xl==4] <- "NOR"


if (delayed) {
if (same.cens==TRUE) {
    etime <- rep(rbinom(n,1,ptrunc)*(runif(n,truncS)),each=2)
} else  etime<- rbinom(2*n,1,ptrunc)*(runif(2*n,truncS))
} else etime <- rep(0,2*n)


data<-data.frame(time=mt*time,cause=cause,xl=rep(xl,each=2),
		 country=rep(country,each=2),id=id,cens=cens,stime=mt*stime,type=rep(types,each=2),
		 f1=rep(F11x,each=2),p11=rep(p11,each=2),p12=rep(p12,each=2),p21=rep(p21,each=2),
		 p22=rep(p22,each=2),entry=mt*etime,truncated=(time<etime))
return(data)
} ## }}}

##' @export 
simnordic <- function(n,cordz=2,cormz=3,
		      censSmz=c(0,1),censSdz=c(0,1),
		      pcensmz=0.5,pcensdz=0.5,
		      ptrunc=0.5,
		      truncSmz=c(0,0.5),truncSdz=c(0,0.5),
		      country=TRUE,same.cens=TRUE,
                      delayed=FALSE,only.delayed=FALSE,lam0=0.13) 
{ ## {{{
outdz <- corsim.prostate(n,theta=cordz,censS=censSdz,pcens=pcensdz,mt=1,same.cens=same.cens,test=0,country=country,
	 delayed=delayed,ptrunc=ptrunc,lam0=lam0,truncS=truncSdz) 
outmz <- corsim.prostate(n,theta=cormz,censS=censSmz,pcens=pcensmz,mt=1,same.cens=same.cens,test=0,country=country,
	 delayed=delayed,ptrunc=ptrunc,lam0=lam0,truncS=truncSmz)
outdz$zyg <- "DZ" 
outmz$zyg <-  "MZ"
outmz$id <- outmz$id+nrow(outdz)
###
out <- rbind(outdz,outmz)
out$time <- out$time*100
out$entry <- out$entry*100
###table(out$type,out$country)
###table(out$type,out$cause)
if (country==TRUE) out$country <- relevel(factor(out$country),ref="SWE")
###table(out$country)
outk <- out[,c("country","cause","id","time","zyg","type","entry","truncated")]

if (only.delayed) outk <- outk[!out$truncated,]

return(outk)
} ## }}}

## {{{ simulation for gamma distributed cif model 

lap<-function(theta,t) { return( (1+t/theta)^(-theta)) }
ilap<-function(theta,t) {
	itheta<-1/theta; return((t^(-itheta)-1)/(itheta)) }

F1clust<-function(t,rtheta=1,theta=1,lam0=0.5,x=0) {
	return(1-exp(-rtheta*ilap(theta,1-F1addfg(t,lam0=lam0,x=x))))
}

F1<-function(t,lam0=0.5,beta=0.3,x=0) # additive version
{ return( 1 - exp(-t*lam0-t*x*beta)) }

F2<-function(t,lam0=0.5,beta=0.3,x=0) # additive version
{ return( 1 - exp(-t*lam0-t*x*beta)) }

sim.F1F2<-function(n,theta=1,lam0=0.1,beta=0.3,lam02=0.1,beta2=0,crate=3,cstart=0) 
{ ## {{{
	x<-rbinom(n,1,0.5); 
	tt<-seq(0,3,length=100)
	F1x<-F1(3,x=x,beta=beta,lam0=lam0)
	F2x<-F2(3,x=x,beta=beta2,lam0=lam02)

	okX <- (F1x+F2x<=1)*1
	n <- sum(okX)
	F1x <- F1x[okX]
	F2x <- F1x[okX]

	### death or alive
	cause12<-rbinom(n,1,F1x+F2x)
	### if death 1 or 2
	cause1e2<-rbinom(n,1,F1x/(F1x+F2x))
	cause <- rep(0,n)
	stime<-rep(3,n); 

	for (i in 1:n)
	{
		if (cause12[i]==1) {
			if (cause1e2[i]==1) {
				cause[i] <- 1
				myhazx<-F1(tt,x=x[i],beta=beta,lam0=lam0)/F1x[i]
				stime[i]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
			} else  {
				cause[i] <- 2
				myhazx<-F2(tt,x=x[i],beta=beta2,lam0=lam02)/F2x[i]
				stime[i]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
			}
		}
	}

	### censoring on top
	ctime<-(cstart+runif(n)*crate)
	time<-apply(cbind(ctime,stime),1,min)
	status<-(stime<ctime); 
	if (sum(status==0)>0) cause[status==0]<-0; 
	data<-data.frame(time=time,ctime=ctime,status=status,X=x,cause=cause)
	return(data)
} ## }}}

sim.F1<-function(n,lam0=0.5,beta=0.3,Cint=c(0,1)) 
{ ## {{{
	x<-runif(n);
	tt<-seq(0,1,length=100)
	F11x<-F1(1,x=x,beta=beta,lam0=lam0)
	cause1<-rbinom(n,1,F11x)
	###
	stime<-rep(2,n); 
	for (i in 1:n) {
		if (cause1[i]==1) {
			myhazx<-F1(tt,x=x[i],beta=beta,lam0=lam0)/F11x[i]
			stime[i]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
		} 
	}
	ctime<-runif(n,Cint)
	time<-pmin(ctime,stime)
	status<-(stime<ctime); 
	cause1[status==0]<-0; 
	data<-data.frame(time=time,status=status,X=x,cause=cause1)
	return(data)
} ## }}}

sim.F1clust<-function(n,theta=1,lam0=0.5,Clims=c(0,1),
		      beta=-0.5,same.cens=FALSE,fix.cens=FALSE)
{ ## {{{ 
	k<-n/2; tt<-seq(0,1,length=100)
	theta <- 1/theta
	rtheta<-rgamma(k,theta,scale=1/theta)
	stime<-c();cause1<-c();id<-c();vtheta<-c(); X<-c()

	for (i in 1:k)
	{  ## {{{  
		x <- rbinom(1,1,0.25)
		x <- rep(x,2)
		X<-c(X,x); 
		F11x<-F1clust(1,rtheta=rtheta[i],theta=theta,x=x,lam0=lam0) 
		cause<-rbinom(2,1,F11x); 
		###cause1<-c(cause1,cause); 
		id<-c(id,rep(i,2)); vtheta<-c(vtheta,rep(rtheta[i],2))

		for (j in 1:2) {
			if (cause[j]==1) {
				myhazx<-F1clust(tt,x=x[j],rtheta=rtheta[i],theta=theta,lam0=lam0)/F11x[j]
				stime<-c(stime,Cpred(cbind(myhazx,tt),runif(1))[1,2]+ runif(1,0,0.001))
				cause1 <- c(cause1,1);
			} else { stime<-c(stime,runif(1)); cause1 <- c(cause1,2);}
		}
	} ## }}} 

	if (same.cens) ctime <- rep(runif(k,min=Clims[1],max=Clims[2]),each=2) else ctime<-runif(n,min=Clims[1],max=Clims[2]);
	if (fix.cens) ctime <- rep(Clims[2],n)

	time<-apply(cbind(ctime,stime),1,min)
	status<-(stime<ctime); 
	cause1[status==0]<-0; 
	data<-data.frame(time=time,X=X,cause=cause1,stime=stime,ctime=ctime,id=id,theta=vtheta)
	return(data)
} ## }}}
## }}}

##' @export 
corsim.prostate.random <- function(n,theta=1,censS=c(0,1),
   pcens=0.5,test=0,mt=1,same.cens=TRUE,country=TRUE,
   delayed=FALSE,ptrunc=0.5,lam0=0.13,truncS=c(0.5,1)) 
{ ## {{{
### n <- 1000; theta <- 2; 
k<-n/2; 
theta <- 1/theta
tt<-seq(0,1,length=100)
rtheta<-rgamma(k,theta,scale=1/theta)
stime<-c();cause1<-c();id<-c();vtheta<-c(); X<-c()

for (i in 1:k)
{  ## {{{  
if (country==TRUE) x <- rbinom(1,1,0.25) else x <- 0
x <- rep(x,2)
X <- c(X,x); 
F11x <- F1clust(1,rtheta=rtheta[i],theta=theta,x=x,lam0=lam0) 
cause <- rbinom(2,1,F11x); 
###cause1<-c(cause1,cause); 
id<-c(id,rep(i,2)); vtheta<-c(vtheta,rep(rtheta[i],2))

for (j in 1:2) {
if (cause[j]==1) {
myhazx<-F1clust(tt,x=x[j],rtheta=rtheta[i],theta=theta,lam0=lam0)/F11x[j]
stime<-c(stime,Cpred(cbind(myhazx,tt),runif(1))[1,2]+ runif(1,0,0.001))
cause1 <- c(cause1,1);
} else { stime<-c(stime,runif(1)); cause1 <- c(cause1,2);}
}
} ## }}} 

###same.cens=TRUE
if (same.cens==TRUE) {
    ctime <- rep(rbinom(n,1,pcens)*runif(n,min=censS[1],max=censS[2]),each=2)
    ctime[ctime==0] <- 1;
}
else {
     ctime<- rbinom(2*n,1,pcens)*runif(2*n,min=censS[1],max=censS[2])
     ctime[ctime==0] <- 1;
}

cens <- (ctime < stime)
time <- ifelse(cens,ctime,stime)
cause <- ifelse(cens,0,cause1)
id <- rep(1:n,each=2)

if (country==TRUE) {
country <- c()
country[X==1] <- "DK"
no <- sum(X==0)
country[X==0] <- rep(sample(c("SWE","FIN","NOR"),no/2,replace=TRUE),each=2)
country <- relevel(factor(country),ref="DK")
} else  country <- rep("same",2*n)

if (delayed) {
if (same.cens==TRUE) etime <- rep(rbinom(n,1,ptrunc)*runif(n,min=truncS[1],max=truncS[2]),each=2)
 else  etime<- rbinom(2*n,1,ptrunc)*(runif(2*n,min=truncS[1],max=truncS[2]))
} else etime <- rep(0,2*n)

data<-data.frame(time=time,cause=cause,x=X,
	 country=country,id=id,cens=cens,
	 stime=stime, entry=etime,truncated=(time<etime))
return(data)
} ## }}}

##' @export 
simnordic.random <- function(n,cordz=2,cormz=3,
     censSmz=c(0,1),censSdz=c(0,1),pcensmz=0.5,pcensdz=0.5,
     ptrunc=0.5,truncSmz=c(0.5,1),truncSdz=c(0.5,1),
     country=TRUE,same.cens=TRUE,
     delayed=FALSE,only.delayed=FALSE,lam0=0.13) 
{ ## {{{
outdz <- corsim.prostate.random(n,theta=cordz,censS=censSdz,
 pcens=pcensdz,mt=1,same.cens=same.cens,test=0,country=country,
 delayed=delayed,ptrunc=ptrunc,lam0=lam0,truncS=truncSdz) 
outmz <- corsim.prostate.random(n,theta=cormz,censS=censSmz,
 pcens=pcensmz,mt=1,same.cens=same.cens,test=0,country=country,
 delayed=delayed,ptrunc=ptrunc,lam0=lam0,truncS=truncSmz)
outdz$zyg <- "DZ" 
outmz$zyg <-  "MZ"
outmz$id <- outmz$id+nrow(outdz)
###
out <- rbind(outdz,outmz)
out$time <- out$time*100
out$entry <- out$entry*100
out$zyg <- relevel(factor(out$zyg),ref="MZ")
if (country==TRUE) out$country <- relevel(factor(out$country),ref="DK")
###
if (only.delayed) out <- out[!out$truncated,]

return(out)
} ## }}}

