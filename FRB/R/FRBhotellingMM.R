FRBhotellingMM <-function(X,...) UseMethod("FRBhotellingMM")

FRBhotellingMM.formula <-function(formula, data=NULL, ...)
{
# --------------------------------------------------------------------

# Returns response of formula in nice way

model.multiregresp<-function (data, type = "any") 
{
    if (attr(attr(data, "terms"), "response")) {
        if (is.list(data) | is.data.frame(data)) {
  		v <- data[[1L]]
		if (is.data.frame(data) && is.vector(v)) v <- data[,1L,drop=FALSE]
            if (type == "numeric" && is.factor(v)) {
                warning("using type=\"numeric\" with a factor response will be ignored")
            }
            else if (type == "numeric" | type == "double") 
                storage.mode(v) <- "double"
            else if (type != "any") 
                stop("invalid response type")
            if (is.matrix(v) && ncol(v) == 1L){ 
                if (is.data.frame(data)) {v=data[,1L,drop=FALSE]}
	          else {dim(v) <- NULL}}
            rows <- attr(data, "row.names")
            if (nrows <- length(rows)) {
                if (length(v) == nrows) 
                  names(v) <- rows
                else if (length(dd <- dim(v)) == 2L) 
                  if (dd[1L] == nrows && !length((dn <- dimnames(v))[[1L]])) 
                    dimnames(v) <- list(rows, dn[[2L]])
            }
            return(v)
        }
        else stop("invalid 'data' argument")
    }
    else return(NULL)
}

    mt <- terms(formula, data = data)
    if (attr(mt, "response") == 0L) stop("response is missing in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    miss <- attr(mf,"na.action")
#print((model.response(mf))
#    X <- as.matrix(model.response(mf))
     X <- model.multiregresp(mf)
#print(X)
    Terms <- attr(mf, "terms")
    Y <- model.matrix(Terms, mf)
    Yint <- match("(Intercept)", colnames(Y), nomatch = 0L)
    if (Yint > 0 && ncol(Y)==1)
	{
	x<-X
	y<-NULL
	}
    else
	{
	if (Yint > 0) Y <- Y[, -Yint, drop = FALSE]
	if (ncol(Y) > 1L) stop("too many predictors in formula")
	Y<- as.factor(Y)
      if (nlevels(Y) != 2) stop("predictor must have two levels")
    	x <- X[which(Y==levels(Y)[1]),,drop=FALSE]
	y <- X[which(Y==levels(Y)[2]),,drop=FALSE]
     }
#print(x)
	res <- FRBhotellingMM.default(x, y, ...)
    	res$terms <- Terms
    	cl <- match.call()
    	cl[[1L]] <- as.name("FRBhotellingMM")
    	res$call <- cl
    	if (!is.null(miss)) res$na.action <- miss
	return(res)
}                                                 



FRBhotellingMM.default <-function(X,Y=NULL,mu0=0,R=999,conf=0.95,method=c("HeFung","pool"),control=MMcontrol(...),na.action=na.omit, ...)
{
# performs robust Hotelling test based on multivariate MM estimates 
# with fast and robust bootstrap
#
# calls: MMest_loccov(), MMboot_loccov(), MMest_twosample(), MMboot_twosample()
#
# Input
# Xdata: (n x p) data set
# mu0: mean under null hypothesis in case of one sample
# R: number of bootstrap samples
# method: "pool" uses the pooled covariance matrix  
#          "HeFung" uses the covariance matrix of He and Fung
#          to estimate the common covariance matrix 
# Output
# p-value from the robust Hotelling test
#-------------------------------------------------------------------------
vecop <- function(mat) {
# performs vec-operation (stacks colums of a matrix into column-vector)

nr <- nrow(mat)
nc <- ncol(mat)

vecmat <- rep(0,nr*nc)
for (col in 1:nc) {
    startindex <- (col-1)*nr+1
    vecmat[startindex:(startindex+nr-1)] <- mat[,col]
}
return(vecmat)
}

# --------------------------------------------------------------------

reconvec <- function(vec,ncol) {
# reconstructs vecop'd matrix

lcol <- length(vec)/ncol
rec <- matrix(0,lcol,ncol)
for (i in 1:ncol)
    rec[,i] <- vec[((i-1)*lcol+1):(i*lcol)]

return(rec)
}

#--------------------------------------------------------------------------
hottest1MM<-function(Xdata,mu0,R,conf,control=MMcontrol())
{
#Robust one-sample Hotelling test based on MM-estimator

    mu0<- as.vector(mu0)
    n <- nrow(Xdata)
    p <- ncol(Xdata)
    dimens <- p+p*p
    

    MMests <- MMest_loccov(Xdata,control=control) 
    MMmu <- MMests$Mu
    MMSigma <- MMests$Sigma
    MMGamma <- MMests$Gamma
    Smu <- MMests$SMu
    SSigma <- MMests$SSigma
    teststatMM<-n*(MMmu -mu0) %*%solve(MMSigma)%*%t(MMmu-mu0)
    testvectMM<-c()
    bootresMM<-MMboot_loccov(Xdata,R=R,ests=MMests)
    bmatrixSenmm <- bootresMM$centered

for (r in 1:R) {
        estimates<-bmatrixSenmm[,r]
        estmubminmuhat<-estimates[1:p]
        estmubminmuhat<-as.matrix(estmubminmuhat)	
        estgammabmingammahat<-estimates[(p+1):dimens]
        estgammab<-estgammabmingammahat+vecop(MMGamma)
        estgammab<-reconvec(estgammab,p)
        estsigmabminsigmahat<-estimates[(dimens+1):(dimens+p*p)]
        estsigmab<-estsigmabminsigmahat+vecop(SSigma)
        estsigmab<-reconvec(estsigmab,p)
        if ((det(estsigmab) >0) && (det(estgammab) >0)) {
            estsigmammb<-det(estsigmab)^(1/p)*estgammab
            testvectMMwaarde<-n*t(estmubminmuhat)%*%solve(estsigmammb)%*%estmubminmuhat
            testvectMM<-c(testvectMM,testvectMMwaarde)
        }    
}    
pvalue<-mean(testvectMM>=as.numeric(teststatMM))
Rok <- length(testvectMM)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvectMM[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
	   conf.int[1,i] <- MMmu[i]-sqrt(1/n*quantval*MMSigma[i,i])
     conf.int[2,i] <- MMmu[i]+sqrt(1/n*quantval*MMSigma[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata)[[2]]) 
rownames(MMmu) <- ("   Estimate")


return(list(teststat=teststatMM,testvect=testvectMM,pvalue=pvalue,loc=MMmu,cov=MMSigma,confint=conf.int,w=MMests$w,outFlag=MMests$outFlag,ROK=Rok))

}

#-----------------------------------------------------------------------------

hottest2MMpool <-function(Xdata1,Xdata2,R,conf,control=MMcontrol())
{
#Robust two-sample Hotelling test based on MM-estimator and pooled covariance matrix

n1 <- nrow(Xdata1)
p <- ncol(Xdata1)
n2 <- nrow(Xdata2)
n <- n1*n2/(n1+n2)
dimens<-p+p*p
MMests1 <- MMest_loccov(Xdata1,control=control)
MMests2 <- MMest_loccov(Xdata2,control=control)
w <- c(MMests1$w, MMests2$w)
outFlag <- c(MMests1$outFlag, MMests2$outFlag)
MMmu1 <- MMests1$Mu
MMSigma1 <- MMests1$Sigma
MMGamma1 <- MMests1$Gamma
Smu1<- MMests1$SMu
SSigma1 <- MMests1$SSigma
MMmu2 <- MMests2$Mu
MMSigma2 <- MMests2$Sigma
Smu2<-MMests2$SMu
SSigma2<-MMests2$SSigma
MMGamma2 <- MMests2$Gamma

MMSigmap<-((n1-1)*MMSigma1+(n2-1)*MMSigma2)/(n1+n2-2)
SSigmap<-((n1-1)*SSigma1+(n2-1)*SSigma2)/(n1+n2-2)
teststatMM<-((n1*n2)/(n1+n2))*(MMmu1-MMmu2)%*%solve(MMSigmap)%*%t(MMmu1-MMmu2)


bootresMM1<-MMboot_loccov(Xdata1,R=R,ests=MMests1)
bootresMM2<-MMboot_loccov(Xdata2,R=R,ests=MMests2)
bmatrixSenmm1 <- bootresMM1$centered
bmatrixSenmm2 <- bootresMM2$centered

testvectMM<-c()

    
for (r in 1:R) {
    estimates1<-bmatrixSenmm1[,r]
    estmubminmuhat1<-estimates1[1:p]
    estmubminmuhat1 <- as.matrix(estmubminmuhat1)
    estgammabmingammahat1<-estimates1[(p+1):dimens]
    estgammab1<-estgammabmingammahat1+vecop(MMGamma1)
    estgammab1<-reconvec(estgammab1,p)
    estsigmabminsigmahat1<-estimates1[(dimens+1):(dimens+p*p)]
    estsigmab1<-estsigmabminsigmahat1+vecop(SSigma1)
    estsigmab1<-reconvec(estsigmab1,p)

    estimates2<-bmatrixSenmm2[,r]
    estmubminmuhat2<-estimates2[1:p]
    estmubminmuhat2 <- as.matrix(estmubminmuhat2)
    estgammabmingammahat2<-estimates2[(p+1):dimens]
    estgammab2<-estgammabmingammahat2+vecop(MMGamma2)
    estgammab2<-reconvec(estgammab2,p)
    estsigmabminsigmahat2<-estimates2[(dimens+1):(dimens+p*p)]
    estsigmab2<-estsigmabminsigmahat2+vecop(SSigma2)
    estsigmab2<-reconvec(estsigmab2,p)
    if (is.finite(det(estsigmab1)) && is.finite(det(estsigmab2)) && is.finite(det(estgammab1)) && is.finite(det(estgammab2)) &&
    (det(estsigmab1) >0) && (det(estsigmab2)>0) && (det(estgammab1) >0) && (det(estgammab2)>0)) {
        estsigmammb1<-det(estsigmab1)^(1/p)*estgammab1
        estsigmammb2<-det(estsigmab2)^(1/p)*estgammab2
        estsigmaMMP<-((n1-1)*estsigmammb1+(n2-1)*estsigmammb2)/(n1+n2-2)
        testvectMMwaarde<-((n1*n2)/(n1+n2))*t(estmubminmuhat1-estmubminmuhat2)%*%solve(estsigmaMMP)%*%(estmubminmuhat1-estmubminmuhat2)

        testvectMM=c(testvectMM,testvectMMwaarde)
    }
}

pvalue<-mean(testvectMM>=as.numeric(teststatMM))
Rok<-length(testvectMM)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvectMM[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
	   conf.int[1,i] <- MMmu1[i]-MMmu2[i]-sqrt(1/n*quantval*MMSigmap[i,i])
	   conf.int[2,i] <- MMmu1[i]-MMmu2[i]+sqrt(1/n*quantval*MMSigmap[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata1)[[2]]) 
rownames(MMmu1) <- rownames(MMmu2) <- ("   Estimate")


return(list(teststat=teststatMM,testvect=testvectMM,pvalue=pvalue,loc1=MMmu1,loc2=MMmu2,cov=MMSigmap,confint=conf.int,w=w,outFlag=outFlag, ROK=Rok))

}

#----------------------------------------------------------------------------

hottest2MMHe<-function(Xdata1,Xdata2,R,conf,control=MMcontrol()){
#Robust two-sample Hotelling test based on MM-estimator of He and Fung 

n1 <- nrow(Xdata1)
n2 <- nrow(Xdata2)
n <- n1*n2/(n1+n2)
p <- ncol(Xdata1)
dimens <- 2*p+p*p
Xdata <- rbind(Xdata1,Xdata2)
groups <- c(rep(1,n1),rep(2,n2))
estsMM <- MMest_twosample(Xdata,groups,control=control)

MMmu1 <- estsMM$Mu1
Smu1 <- estsMM$SMu1
MMmu2 <- estsMM$Mu2
Smu2 <- estsMM$SMu2
SSigmap <- estsMM$SSigma
MMGammap <- estsMM$Gamma
MMSigmap <- estsMM$Sigma

teststatMM <-((n1*n2)/(n1+n2))*(MMmu1-MMmu2)%*%solve(MMSigmap)%*%t(MMmu1-MMmu2)

bootresMM <- MMboot_twosample(Xdata, groups=groups, R, ests=estsMM)
bmatrixSenmm <- bootresMM$centered
testvectMM <- c()

for (r in 1:R) {
    estimates <- bmatrixSenmm[,r]
    estmubminmuhat1 <- estimates[1:p]
    estmubminmuhat2 <- estimates[(p+1):(2*p)]
    estmubminmuhat1 <- as.matrix(estmubminmuhat1)
    estmubminmuhat2 <- as.matrix(estmubminmuhat2)
    estgammabmingammahat <- estimates[(2*p+1):(2*p+p*p)]
    estgammaMMP <- estgammabmingammahat+vecop(MMGammap)
    estgammaMMP <- reconvec(estgammaMMP,p)
    estsigmabminsigmahatSP <- estimates[(2*p+p*p+1):(2*p+2*p*p)]
    estsigmabSP <- estsigmabminsigmahatSP+vecop(SSigmap)
    estsigmabSP <- reconvec(estsigmabSP,p)
    if (is.finite(det(estsigmabSP)) && is.finite(det(estgammaMMP)) && (det(estsigmabSP) >0) && (det(estgammaMMP) >0)) {
        estsigmaMMP <- det(estsigmabSP)^(1/p)*estgammaMMP
        testvectMMwaarde<-((n1*n2)/(n1+n2))*t(estmubminmuhat1-estmubminmuhat2)%*%solve(estsigmaMMP)%*%(estmubminmuhat1-estmubminmuhat2)
        testvectMM=c(testvectMM,testvectMMwaarde)
    }
}

pvalue<-mean(testvectMM>=as.numeric(teststatMM))
Rok<-length(testvectMM)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvectMM[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
    conf.int[1,i] <- MMmu1[i]-MMmu2[i]-sqrt(1/n*quantval*MMSigmap[i,i])
		conf.int[2,i] <- MMmu1[i]-MMmu2[i]+sqrt(1/n*quantval*MMSigmap[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata)[[2]]) 
rownames(MMmu1) <- rownames(MMmu2) <- ("   Estimate")

return(list(teststat=teststatMM,testvect=testvectMM,pvalue=pvalue,loc1=MMmu1,loc2=MMmu2,cov=MMSigmap,confint=conf.int,w=estsMM$w,outFlag=estsMM$outFlag,ROK=Rok))
 
}

#----------------------------------------------------------------------------
#-                         main function                                     -
#----------------------------------------------------------------------------

method <- match.arg(method)

Xdatam <- as.matrix(X)
xnam=colnames(Xdatam)
Xdatam=na.action(Xdatam)

p <- ncol(Xdatam)
n <- nrow(Xdatam)

if (p < 1L) stop("at least one variable needed")
if (n < p) stop("For Hotelling tests the number of observations cannot be smaller than the 
number of variables")

if (is.null(xnam))
    colnames(Xdatam) <- paste("V",1:p,sep="")

if (is.null(Y)==TRUE)
 {nrsamples <- 1}
else
 {Ydatam <- as.matrix(Y)
  Ydatam=na.action(Ydatam)
  q <- ncol(Ydatam)
  if (p != q) stop("For two-sample Hotelling test both samples must contain the same variables")
  if (is.null(colnames(Ydatam)))
    colnames(Ydatam) <- paste("V",1:p,sep="")
  nrsamples <- 2}

if (length(mu0)==1)
      {mu0=rep(mu0,p)}


if (nrsamples==1)
   {
   meth=paste("One sample Hotelling test based on multivariate MM-estimates (bdp = ", control$bdp,", eff = ", control$eff, ")", sep="")
   res=hottest1MM(Xdatam,mu0=mu0,R=R,conf=conf,control=control)
   names(res$teststat)="T^2_R"
   rownames(res$loc)="MM-loc vector"
   res$alt=paste("true mean vector is not equal to",
   paste("(", paste(round(mu0, digits=3), collapse = ","), ")", sep = ""), "\n")
   z <- list(p.value=res$pvalue,statistic=res$teststat,teststat.boot=res$testvect,
   estimate=round(res$loc,digits=3),Sigma=res$cov,alternative=res$alt,
   CI=res$confint,Mu0=mu0,conf=conf,data.name=deparse(substitute(X)),
   method=meth,X=Xdatam,w=res$w,outFlag=res$outFlag,ROK=res$ROK)
}
else 
    {
    if (method == "pool") {
      meth=paste("Two sample Hotelling test based on multivariate MM-estimates (bdp = ", control$bdp,", eff = ", control$eff, ")\n",
          "(common covariance estimated by pooled covariance matrix)", sep="")
      res=hottest2MMpool(Xdatam,Ydatam,R=R,conf=conf,control=control)
    }
    else { 
      meth=paste("Two sample Hotelling test based on multivariate MM-estimates (bdp = ", control$bdp,", eff = ", control$eff, ")\n",
          "(common covariance estimated by He and Fung method)", sep="")
      res=hottest2MMHe(Xdatam,Ydatam,R=R,conf=conf,control=control)  
    }
   names(res$teststat)="T^2_R"
   res$loc=rbind(res$loc1,res$loc2)
   rownames(res$loc)=c("MM-loc x-vector","MM-loc y-vector")
   res$alt=paste("true difference in mean vectors is not equal to",
   paste("(", paste(round(mu0, digits=3), collapse = ","), ")", sep = ""), "\n")    
   z <- list(p.value=res$pvalue,statistic=res$teststat,teststat.boot=res$testvect,
   Mu1=res$loc1,Mu2=res$loc2,estimate=round(res$loc,digits=3),Sigma=res$cov,
   alternative=res$alt,CI=res$confint,conf=conf,
   data.name=paste(deparse(substitute(x)),"and",deparse(substitute(y))),
   method=meth,X=Xdatam,Y=Ydatam,w=res$w,outFlag=res$outFlag,ROK=res$ROK)
}

#class(z) <- "FRBhot"
class(z)=c("htest","FRBhot")
return(z)

}  
#-----------------------------------------------------------------------------------------









