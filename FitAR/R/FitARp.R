FitARp <-
function(z,p,lag.max="default",MLEQ=FALSE)
{
stopifnot(length(z)>0, length(z)>length(p), length(p)>0)
is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
stopifnot(is.wholenumber(p), p>0)
n<-length(z)
if (lag.max=="default")
    MaxLag <- min(300, ceiling(length(z)/5))
else
    MaxLag <- lag.max
pvec <- sort(p)
pvec<-pvec[pvec>0]
if (length(pvec)==0)
    pvec<-0
PMAX<-max(pvec)
if (length(p)==1 && PMAX!=0)  
    pvec <- 1:p
SubsetQ <- length(pvec)<PMAX
if (PMAX == 0) SubsetQ<-FALSE
mz <- mean(z) 
y <- z 
#get parameter estimates
if (MLEQ){
    ans<-GetFitARpMLE(y,pvec)
    FitMethod<-"MLE"
    MeanMLEQ<-TRUE    
    }
else {
    ans<-GetFitARpLS(y,pvec)
    FitMethod<-"LS"
    MeanMLEQ<-FALSE
    }
phiHat<-ans$phiHat
res<-BackcastResidualsAR(y-mz, phiHat)
fits<-y-res
sigsq<-sum(res^2)/n
racf<-(acf(res, plot=FALSE, lag.max=MaxLag)$acf)[-1]
#covariance matrix via inverse Fisher information matrix
#sd of racf
if (SubsetQ){
    varNames<-paste("phi(",pvec,")",sep="")
    covHat<-solve(InformationMatrixARp(phiHat,pvec))/n
    dimnames(covHat)<-list(varNames,varNames)
    sdRacf<-sqrt(diag(VarianceRacfARp(phiHat,pvec,MaxLag,n)))
}
else {
    if (PMAX>0) {
        varNames<-paste("phi(",1:PMAX,")",sep="")
        covHat<-SiddiquiMatrix(phiHat)/n
        dimnames(covHat)<-list(varNames,varNames)
        sdRacf<-sqrt(diag(VarianceRacfAR(phiHat,MaxLag,n)))
        }
    else {
        varNames<-character(0)
        covHat<-numeric(0)
        sdRacf<-rep(1/sqrt(n),MaxLag)
    }
}
if (SubsetQ) {
    ModelTitle<-deparse(as.numeric(pvec),width.cutoff=180)
    ModelTitle<-paste("ARp",substr(ModelTitle,2,nchar(ModelTitle)),sep="")
    ModelTitle<-gsub(" ", "", ModelTitle)
}
else 
    ModelTitle<-paste("AR(",p,")",sep="")
#
LBQ<-LjungBoxTest(res, lag.max=MaxLag, k=length(pvec))
RacfMatrix<-matrix(c(racf,sdRacf),ncol=2)
dimnames(RacfMatrix)<-list(1:MaxLag, c("ra", "Sd(ra)"))
zetaHat<-ARToPacf(phiHat)
#
if (!MLEQ) { #for LS, report usual LS covariance matrix
    covHat<-ans$covHat
    covHat <- covHat[-1,-1,drop=FALSE] #remove intercept
    dimnames(covHat)<-list(varNames,varNames)
    }
ans<-list(loglikelihood=ans$loglikelihood,phiHat=phiHat,sigsqHat=sigsq,muHat=mz,covHat=covHat,zetaHat=zetaHat,
          RacfMatrix=RacfMatrix,LjungBoxQ=LBQ,res=res,fits=fits,SubsetQ=SubsetQ,pvec=pvec,FitMethod=FitMethod,
          MeanMLE=MeanMLEQ, ARModel="ARp", tsp=tsp(z),
          call=match.call(),DataTitle=attr(z,"title"),ModelTitle=ModelTitle,z=z)
class(ans)<-"FitAR"
ans
}

