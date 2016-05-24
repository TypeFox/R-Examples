FRBhotellingS <-function(X,...) UseMethod("FRBhotellingS")

FRBhotellingS.formula <-function(formula, data=NULL, ...)
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
#    X <- as.matrix(model.response(mf))
     X <- model.multiregresp(mf)
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
	res <- FRBhotellingS.default(x, y, ...)
    	res$terms <- Terms
    	cl <- match.call()
    	cl[[1L]] <- as.name("FRBhotellingS")
    	res$call <- cl
    	if (!is.null(miss)) res$na.action <- miss
	return(res)
}                                                 


FRBhotellingS.default <-function(X,Y=NULL,mu0=0,R=999,bdp=0.5,conf=0.95,method=c("HeFung","pool"),control=Scontrol(...),na.action=na.omit, ...)
{
# performs robust Hotelling test based on multivariate S estimates 
# with fast and robust bootstrap
#
# calls: Sest_loccov(), Sboot_loccov(), Sest_twosample(), Sboot_twosample()
#
# Input
# Xdata: (n x p) data set
# groups:  vector of 1's and 2's, indicating group numbers
# mu0: mean under null hypothesis in case of one sample
# bdp: breakdown point of the robust estimator
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
hottest1S <-function(Xdata,mu0=mu0,bdp,R,conf,control=Scontrol())
{
#Robust one-sample Hotelling test based on S-estimator 
 

mu0<- as.vector(mu0)
n <- nrow(Xdata)
p <- ncol(Xdata)
dimens<-p+p*p

Sests <- Sest_loccov(Xdata,bdp=bdp,control=control)
Smu <- Sests$Mu
SSigma <- Sests$Sigma
teststat<-n*(Smu -mu0) %*%solve(SSigma)%*%t(Smu-mu0)
testvect<-c()
bootresS<-Sboot_loccov(Xdata,R=R,ests=Sests)
bmatrixs <- bootresS$centered
    
for (r in 1:R){
    estimates<-bmatrixs[,r]
    estmubminmuhat<-estimates[1:p]
    estmubminmuhat <- as.matrix(estmubminmuhat)
    estsigmabminsigmahat<-estimates[(p+1):dimens]
    estsigmab<-estsigmabminsigmahat+vecop(SSigma)
    estsigmab<-reconvec(estsigmab,p)
    if (det(estsigmab) >0)
        {testvectwaarde<-n*t(estmubminmuhat)%*%solve(estsigmab)%*%estmubminmuhat
        testvect<-c(testvect,testvectwaarde)
    }  
} 
pvalue<-mean(testvect>=as.numeric(teststat))
Rok<-length(testvect)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvect[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
 		conf.int[1,i] <- Smu[i]-sqrt(1/n*quantval*SSigma[i,i])
		conf.int[2,i] <- Smu[i]+sqrt(1/n*quantval*SSigma[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata)[[2]]) 
rownames(Smu) <- ("   Estimate")

return(list(teststat=teststat,testvect=testvect,pvalue=pvalue,loc=Smu,cov=SSigma,confint=conf.int,w=Sests$w,outFlag=Sests$outFlag,ROK=Rok))

}
#-------------------------------------------------------------------------------

hottest2Spool <- function(Xdata1,Xdata2,bdp,R,conf,control=Scontrol())
{
#Robust two-sample Hotelling test based on S-estimator and pooled covariance matrix


n1 <- nrow(Xdata1)
n2 <- nrow(Xdata2)
n <- n1*n2/(n1+n2)
p <- ncol(Xdata1)
dimens <- p+p*p
Sests1 <- Sest_loccov(Xdata1,bdp=bdp,control=control)
Sests2 <- Sest_loccov(Xdata2,bdp=bdp,control=control)
w <- c(Sests1$w, Sests2$w)
outFlag <- c(Sests1$outFlag, Sests2$outFlag)
Smu1 <- Sests1$Mu
SSigma1 <- Sests1$Sigma
Smu2 <- Sests2$Mu
SSigma2 <- Sests2$Sigma
SSigmap <- ((n1-1)*SSigma1+(n2-1)*SSigma2)/(n1+n2-2)
teststatS <- ((n1*n2)/(n1+n2))*(Smu1-Smu2)%*%solve(SSigmap)%*%t(Smu1-Smu2)

bootres1 <- Sboot_loccov(Xdata1,R,ests=Sests1)
bootres2 <- Sboot_loccov(Xdata2,R,ests=Sests2)
bmatrixs1 <- bootres1$centered
bmatrixs2 <- bootres2$centered
testvectS <- c()
    
 for (r in 1:R) {
     estimates1<-bmatrixs1[,r]
     estmubminmuhat1<-estimates1[1:p]
     estmubminmuhat1 <- as.matrix(estmubminmuhat1)
     estsigmabminsigmahat1<-estimates1[(p+1):dimens]
     estsigmab1<-estsigmabminsigmahat1+vecop(SSigma1)
     estsigmab1<-reconvec(estsigmab1,p)

     estimates2<-bmatrixs2[,r]
     estmubminmuhat2<-estimates2[1:p]
     estmubminmuhat2 <- as.matrix(estmubminmuhat2)
     estsigmabminsigmahat2<-estimates2[(p+1):dimens]
     estsigmab2<-estsigmabminsigmahat2+vecop(SSigma2)
     estsigmab2<-reconvec(estsigmab2,p)
     if (is.finite(det(estsigmab1)) && is.finite(det(estsigmab2)) &&
(det(estsigmab1) >0) && (det(estsigmab2)>0)) {
         estsigmaSP<-((n1-1)*estsigmab1+(n2-1)*estsigmab2)/(n1+n2-2)
         testvectSwaarde<-((n1*n2)/(n1+n2))*t(estmubminmuhat1-estmubminmuhat2)%*%solve(estsigmaSP)%*%(estmubminmuhat1-estmubminmuhat2)
         testvectS<-c(testvectS,testvectSwaarde)
     }
 }
 
pvalue <- mean(testvectS>=as.numeric(teststatS))
Rok <- length(testvectS)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvectS[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
 		conf.int[1,i] <- Smu1[i]-Smu2[i]-sqrt(1/n*quantval*SSigmap[i,i])
		conf.int[2,i] <- Smu1[i]-Smu2[i]+sqrt(1/n*quantval*SSigmap[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata1)[[2]]) 
rownames(Smu1) <- rownames(Smu2) <- ("   Estimate")

return(list(teststat=teststatS,testvect=testvectS,pvalue=pvalue,loc1=Smu1,loc2=Smu2,cov=SSigmap,confint=conf.int,w=w,outFlag=outFlag,ROK=Rok))

}

#-------------------------------------------------------------------------------

hottest2SHe<-function(Xdata1,Xdata2,bdp,R,conf,control=Scontrol())
{
#Robust two-sample Hotelling test based on S-estimator of He and Fung

n1 <- nrow(Xdata1)
n2 <- nrow(Xdata2)
n <- n1*n2/(n1+n2)
p <- ncol(Xdata1)
dimens <- 2*p+p*p
Xdata <- rbind(Xdata1,Xdata2)
groups <- c(rep(1,n1),rep(2,n2))
Sests <- Sest_twosample(Xdata, groups, bdp, control=control)
Smu1 <- Sests$Mu1
Smu2 <- Sests$Mu2
SSigmap <- Sests$Sigma

teststatS <- ((n1*n2)/(n1+n2))*(Smu1-Smu2)%*%solve(SSigmap)%*%t(Smu1-Smu2)
bootres <- Sboot_twosample(Xdata, groups,R,ests=Sests)
bmatrixs <- bootres$centered

testvectS <- c()

for (r in 1:R) {
    estimates <- bmatrixs[,r]
    estmubminmuhat1 <- estimates[1:p]
    estmubminmuhat2 <- estimates[(p+1):(2*p)]
    estmubminmuhat1 <- as.matrix(estmubminmuhat1)
    estmubminmuhat2 <- as.matrix(estmubminmuhat2)
    estsigmabminsigmahatP <- estimates[(2*p+1):dimens]
    estsigmaSP <- estsigmabminsigmahatP+vecop(SSigmap)
    estsigmaSP <- reconvec(estsigmaSP,p)
    if (is.finite(det(estsigmaSP)) && (det(estsigmaSP)) >0) {
        testvectSwaarde <-((n1*n2)/(n1+n2))*t(estmubminmuhat1-estmubminmuhat2)%*%solve(estsigmaSP)%*%(estmubminmuhat1-estmubminmuhat2)

        testvectS<-c(testvectS,testvectSwaarde)
        }
}

pvalue <- mean(testvectS>=as.numeric(teststatS))
Rok <- length(testvectS)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvectS[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
 		conf.int[1,i] <- Smu1[i]-Smu2[i]-sqrt(1/n*quantval*SSigmap[i,i])
		conf.int[2,i] <- Smu1[i]-Smu2[i]+sqrt(1/n*quantval*SSigmap[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata)[[2]]) 
rownames(Smu1) <- rownames(Smu2) <- ("   Estimate")

return(list(teststat=teststatS,testvect=testvectS,pvalue=pvalue,loc1=Smu1,loc2=Smu2,cov=SSigmap,confint=conf.int,w=Sests$w, outFlag=Sests$outFlag,ROK=Rok))
}


#----------------------------------------------------------------------------
#-                         main function                                     -
#----------------------------------------------------------------------------


method <- match.arg(method)

Xdatam <- as.matrix(X)
Xdatam=na.action(Xdatam)

p <- ncol(Xdatam)
n <- nrow(Xdatam)

if (p < 1L) stop("at least one variable needed")
if (n < p) stop("For Hotelling tests the number of observations cannot be smaller than the 
number of variables")

if (is.null(colnames(Xdatam)))
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

if (nrsamples ==1)
   {
   meth = paste("One sample Hotelling test based on multivariate S-estimates (breakdown point = ", bdp, ")", sep="")
   res=hottest1S(Xdatam,mu0=mu0,bdp=bdp,R=R,conf=conf,control=control)
   names(res$teststat)="T^2_R"
   rownames(res$loc)="S-loc vector"
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
      meth = paste("Two sample Hotelling test based on multivariate S-estimates (breakdown point = ", bdp, ")\n",
          "(common covariance estimated by pooled covariance matrix)", sep="")
      res=hottest2Spool(Xdatam,Ydatam,bdp=bdp,R=R,conf=conf,control=control)
    }
    else {
      meth = paste("Two sample Hotelling test based on multivariate S-estimates (breakdown point = ", bdp, ")\n",
          "(common covariance estimated by He and Fung method)", sep="")
      res=hottest2SHe(Xdatam,Ydatam,bdp=bdp,R=R,conf=conf,control=control)
    }
   names(res$teststat)="T^2_R"
   res$loc=rbind(res$loc1,res$loc2)
   rownames(res$loc)=c("S-loc x-vector","S-loc y-vector")
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










