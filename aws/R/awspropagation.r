awstestprop <- function(dy,hmax,theta=1,family="Gaussian",
                 lkern="Triangle",aws=TRUE,memory=FALSE,shape=2,
                 homogeneous=TRUE,varadapt=FALSE,ladjust=1,spmin=0.25,seed=1,
                 minlevel=1e-6,maxz=25,diffz=.5,maxni=FALSE,verbose=FALSE){
if(length(dy)>3) {
   cat("maximum array dimension is 3\n contents of first argument will be interpreted as array of
   parameters\n")
   nnn <- length(dy)
} else {
   nnn <- prod(dy)
}
if(minlevel < 5/nnn) {
minlevel <- 5/nnn
cat("minlevel reset to",minlevel,"due to insufficient size of test sample\n")
}
set.seed(seed)
par(mfrow=c(1,1),mar=c(3,3,3,3),mgp=c(2,1,0))
if(length(dy)<=3){
   y <- array(switch(family,"Gaussian"=rnorm(nnn),
                         "Poisson"=rpois(nnn,theta),
                         "Exponential"=rexp(nnn,1),
                         "Bernoulli"=rbinom(nnn,1,theta),
                         "Volatility"=rnorm(nnn),
                         "Variance"=rchisq(nnn,shape)/shape,
                         "NCchi"=sqrt(rchisq(nnn,shape,theta^2))),dy)
} else {
   ddy <- if(!is.null(dim(dy))) dim(dy) else length(dy)
   y <- array(switch(family,"Gaussian"=rnorm(nnn,dy-mean(dy)),
                         "Poisson"=rpois(nnn,dy-mean(dy)+theta),
                         "Exponential"=rexp(nnn,dy-mean(dy)+1),
                         "Bernoulli"=rbinom(nnn,1,dy-mean(dy)+theta),
                         "Volatility"=rnorm(nnn,dy-mean(dy)),
                         "Variance"=rchisq(nnn,dy-mean(dy)+shape)/(dy-mean(dy)+shape),
                         "NCchi"=sqrt(rchisq(nnn,shape,(dy-mean(dy)+theta)^2))),ddy)
   dy <- ddy
}
cat("minlevel ",minlevel,"due to insufficient size of test sample\n")
if(family=="NCchi"){
varstats <- sofmchi(shape/2) # precompute table of mean, sd and var for 
#
#   NCchi for noncentral chi with shape=degrees of freedom and theta =NCP
#
}
z <- seq(0,maxz,diffz)
nz <- length(z)
elevel <- trunc(log(1e-6,10))
levels <- as.vector(outer(c(.5,.2,.1),10^(-0:elevel)))
levels <- levels[levels>=minlevel]
wghts <- switch(length(dy),c(0,0),c(1,0),c(1,1))
cpar<-setawsdefaults(dy,mean(y),family,lkern,"Uniform",aws,memory,ladjust,hmax,shape,wghts)
lambda <- cpar$lambda
hmax <- cpar$hmax
shape <- cpar$shape
d <- cpar$d
n<-length(y)
if(!homogeneous&family=="Gaussian"){
sigma2 <- array(rchisq(prod(dy),shape)/shape,dy)
} else sigma2 <- 1
zfamily <- awsfamily(family,y,sigma2,shape,0,lambda,cpar)
cpar <- zfamily$cpar
lambda <- zfamily$lambda
sigma2 <- zfamily$sigma2
h0 <- zfamily$h0
y <- zfamily$y
lkern <- cpar$lkern
rm(zfamily)
n1 <- switch(d,n,dy[1],dy[1])
n2 <- switch(d,1,dy[2],dy[2])
n3 <- switch(d,1,1,dy[3])
maxvol <- cpar$maxvol
# k <- cpar$k
k <- 1  
kstar <- cpar$kstar
h <- numeric(kstar)
if(k>1) h[1:(k-1)] <- 1+(0:(k-2))*.001
exceedence  <- exceedencena  <- matrix(0,nz,kstar) # this is used to store exceedence probabilities for adaptive and nonadaptive estimates
zobj<-zobj0<-list(ai=y, bi=rep(1,n))
bi <- rep(1,n)
yhat <- y/shape
hhom <- rep(1,n)
lambda0<-1e50
total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
#
#  get initial conditions for a comparison
#
if(family=="Bernoulli") y0 <- (10*y+1)/12
if(family=="Poisson") y0 <- y+.1
#
#  this corresponds to the regularization used to avoid Inf distances
#
kldistnorm1 <- function(th1,y,df){
L <- df/2
m1 <- sqrt(pi/2)*gamma(L+1/2)/gamma(1.5)/gamma(L)*hyperg_1F1(-0.5,L, -th1^2/2, give=FALSE, strict=TRUE)
(m1-y)^2/2/(2*L+th1^2-m1^2)
}
KLdist0 <- switch(family,"Gaussian"=y^2/2,
                         "Poisson"=(theta-y0+y0*(log(y0)-log(theta))),
                         "Exponential"=(log(y)-1+1/y),
                         "Bernoulli"=(y0*log(y0/theta)+
                                     (1-y0)*log((1-y0)/(1-theta))),
                         "Volatility"=(log(y)-1+1/y)/2,
                         "Variance"=shape/2*(log(y)-1+1/y),
                         "NCchi"=kldistnorm1(theta,y,shape)
                         )
exceedence0 <- .Fortran("exceed",
                           as.double(KLdist0),
                           as.integer(length(KLdist0)),
                           as.double(z),
                           as.integer(nz),
                           exprob=double(nz),
                           PACKAGE="aws",DUP=FALSE)$exprob
#
#  now iterate
#
t0 <- Sys.time()
cat("using lambda=",lambda,"\n")
while (k<=kstar){
      t1 <- Sys.time()
      hakt0 <- gethani(1,1.25*hmax,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,1.25*hmax,lkern,1.25^k,wghts,1e-4)
      cat("step",k,"hakt",hakt,"\n")
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hakt <- hakt*0.42445*4
    }
h[k] <- hakt
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
#
#   get nonadaptive estimate
#
if(!homogeneous&family=="Gaussian"){
zobj0 <- .Fortran("chaws1",as.double(y),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       bi=double(n),
                       double(n),
                       double(n),
                       double(n),#vred
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","ai","hakt")]
} else {
zobj0 <- .Fortran("caws1",as.double(y),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       bi=double(n),
                       double(n),
                       double(n),
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","ai","hakt")]
}
if(family%in%c("Bernoulli","Poisson")) zobj0<-regularize(zobj0,family)
yhat0 <- zobj0$ai/zobj0$bi
dim(yhat0) <- dy
bi <- zobj0$bi
ih <- as.integer(hakt)
ind1 <- (ih+1):(dy[1]-ih)
if(length(dy)>1) ind2 <- (ih+1):(dy[2]-ih)
if(length(dy)>2) ind3 <- (ih+1):(dy[3]-ih)
yhat0 <- switch(length(dy),yhat0[ind1],yhat0[ind1,ind2],yhat0[ind1,ind2,ind3])
ni <- max(zobj0$bi)
KLdist0 <- switch(family,"Gaussian"=yhat0^2/2,
                         "Poisson"=(theta-yhat0+yhat0*(log(yhat0)-log(theta))),
                         "Exponential"=(log(yhat0)-1+1/yhat0),
                         "Bernoulli"=(yhat0*log(yhat0/theta)+
                                     (1-yhat0)*log((1-yhat0)/(1-theta))),
                         "Volatility"=(log(yhat0)-1+1/yhat0)/2,
                         "Variance"=shape/2*(log(yhat0)-1+1/yhat0),
                         "NCchi"=kldistnorm1(theta,yhat0,shape))
exceedencena[,k] <- .Fortran("exceed",
                           as.double(KLdist0),
                           as.integer(length(KLdist0)),
                           as.double(z/ni),
                           as.integer(nz),
                           exprob=double(nz),
                           PACKAGE="aws",DUP=FALSE)$exprob
#
#   get adaptive estimate
#
if(!homogeneous&family=="Gaussian"){
zobj <- .Fortran("chaws",as.double(y),
                       as.logical(rep(FALSE,n)),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(yhat),
                       bi=as.double(bi),
                       bi2=double(n),
                       double(n),
                       double(n),#vred
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi2","ai","hakt")]
} else {
if(cpar$mcode!=6){
   zobj <- .Fortran("caws",as.double(y),
                       as.logical(rep(FALSE,n)),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(yhat),
                       bi=as.double(bi),
                       bi2=double(n),
                       double(n),
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi2","ai","hakt")]
} else {
   zobj <- .Fortran("caws6",as.double(y),
                       as.logical(rep(FALSE,n)),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(yhat),
                       as.double(fncchiv(yhat,varstats)/2),
                       bi=as.double(bi),
                       bi2=double(n),
                       double(n),
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi2","ai","hakt")]
}
}
if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
dim(zobj$ai)<-dy
yhat <-zobj$ai/zobj$bi
dim(yhat)<-dy
if(varadapt) bi <- bi^2/zobj$bi2
if(maxni) bi <- pmax(bi,zobj$bi) else bi <- zobj$bi
lambda0 <- lambda
yhat0 <- switch(length(dy),yhat[ind1],yhat[ind1,ind2],yhat[ind1,ind2,ind3])
KLdist1 <- switch(family,"Gaussian"=yhat0^2/2,
                         "Poisson"=(theta-yhat0+yhat0*(log(yhat0)-log(theta))),
                         "Exponential"=(log(yhat0)-1+1/yhat0),
                         "Bernoulli"=(yhat0*log(yhat0/theta)+
                                     (1-yhat0)*log((1-yhat0)/(1-theta))),
                         "Volatility"=(log(yhat0)-1+1/yhat0)/2,
                         "Variance"=shape/2*(log(yhat0)-1+1/yhat0),
                         "NCchi"=kldistnorm1(theta,yhat0,shape))
exceedence[,k] <- .Fortran("exceed",
                           as.double(KLdist1),
                           as.integer(length(KLdist1)),
                           as.double(z/ni),
                           as.integer(nz),
                           exprob=double(nz),
                           PACKAGE="aws",DUP=FALSE)$exprob
                           
contour(z,0:k,cbind(exceedence0,exceedence[,1:k]),levels=levels,ylab="step",xlab="z",
       main=paste(family,length(dy),"-dim. ladj=",ladjust," Exceed. Prob."))
contour(z,0:k,cbind(exceedence0,exceedencena[,1:k]),levels=levels,ylab="step",xlab="z",
       add=TRUE,col=2,lty=3)
       yaxp <- par("yaxp")
       at <- unique(as.integer(seq(yaxp[1],yaxp[2],length=yaxp[3]+1)))
       at <- at[at>0&at<=k]
       axis(4,at=at,labels=as.character(signif(h[at],3)))
       mtext("bandwidth",4,1.8)
if (max(total) >0) {
      cat("Progress:",signif(total[k],2)*100,"% . ",sep="")
     }
t2 <- Sys.time()
tpar <- if(family%in%c("Bernoulli","Poisson"))  paste("theta=",theta) else  ""
cat(family,"(dim:",length(dy),tpar,") ni=",ni," Time: Step",format(signif(difftime(t2,t1),3)),"Total",format(signif(difftime(t2,t0),3)),"\n")
k <- k+1
gc()
}
if(family%in%c("Bernoulli","Poisson")) y <- y0

list(h=h,z=z,prob=exceedence,probna=exceedencena,y=if(verbose) y else NULL , theta= if(verbose) yhat else NULL, levels=levels, family=family)
}

awsweights <- function(awsobj,spmin=0.25){
if(awsobj@degree!=0 || awsobj@varmodel!="Constant"||
any(awsobj@scorr!=0)) stop("Not implemented")
dy <- awsobj@dy
n1 <- dy[1]
ldy <- length(dy)
if(is.null(ldy)) ldy <- 1
if(ldy>1) n2 <- dy[2] else n2 <- 1
if(ldy==3) n3 <- dy[3] else n3 <- 1
hakt <- awsobj@hmax
n <- n1*n2*n3
lambda0 <- awsobj@lambda
yhat <- awsobj@theta
bi <- awsobj@ni
mcode <- switch(awsobj@family,
                Gaussian=1,
		Bernoulli=2,
		Poisson=3,
		Exponential=4,
		Volatility=4,
		Variance=5)
lkern <- awsobj@lkern
hakt <- rep(hakt,ldy)
dlw<-(2*trunc(hakt)+1)
zobj <- .Fortran("cawsw",
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(yhat),
                       bi=as.double(bi),
                       as.integer(mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       wghts=double(n*n),
                       PACKAGE="aws",DUP=TRUE)$wghts
array(zobj,c(dy,dy))
}