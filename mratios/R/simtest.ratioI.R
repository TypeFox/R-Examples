"simtest.ratioI" <-
function(Response, Treatment, alternative="two.sided", Margin.vec=NULL, FWER=0.05, Num.Contrast, Den.Contrast) {

CMat <- Num.Contrast
DMat <- Den.Contrast

n.Treat <- tapply(Response,Treatment,length)
ybar.Treat <- tapply(Response,Treatment,mean)
var.Treat <-  tapply(Response,Treatment,var)
d.freedom <- sum(n.Treat-1) 
Var.pooled <- sum(var.Treat*(n.Treat - 1))/d.freedom


if(!is.numeric(FWER) | length(FWER)!=1 | FWER<=0 | FWER>=0.5)
 {stop("Argument 'FWER' must be a single numeric value between 0 and 0.5")}

if(any( sqrt(var.Treat) < 10 * .Machine$double.eps * abs(ybar.Treat))) 
 {warning("Data are essentially constant in a least one group")}

if(any( n.Treat < 2 )) 
 {warning("There are less than 2 observations in a least one group")}

if(any( ybar.Treat<0 ))
 {warning("At least one sample mean is negative. Note, that with negative denominators in the ratio of interest, tests with one-sided alternatives may test the incorrect direction.")}


MM <- diag(1/n.Treat)    #  Diagonal matrix containing reciprocals of the ni's

ncomp <- nrow (CMat)    # Number of comparisons 

Ratio.Estimate <- Test.Stat <- P.raw <- P.adjusted <-rep(NA,ncomp)
#
#  Correlation matrix under H0
#
CorrMat.H0 <- matrix(rep(NA,ncomp*ncomp),nrow=ncomp)
    for(i in 1:ncomp) {
        for(j in 1:ncomp) {
            CorrMat.H0[i,j] <- (Margin.vec[i]*DMat[i,] - CMat[i,])%*%MM%*%(Margin.vec[j]*DMat[j,] - CMat[j,])/
            (sqrt((Margin.vec[i]*DMat[i,] - CMat[i,])%*%MM%*%(Margin.vec[i]*DMat[i,] - CMat[i,]))*
            sqrt((Margin.vec[j]*DMat[j,] - CMat[j,])%*%MM%*%(Margin.vec[j]*DMat[j,] - CMat[j,])))
        }
    }

for (i in 1:ncomp){
    Ratio.Estimate[i] <- (CMat[i,]%*%ybar.Treat)/(DMat[i,]%*%ybar.Treat)
    Test.Stat[i] <- ((CMat[i,] - Margin.vec[i]*DMat[i,])%*%ybar.Treat)/
              sqrt(Var.pooled*(CMat[i,] - Margin.vec[i]*DMat[i,])%*%MM%*%(CMat[i,] - Margin.vec[i]*DMat[i,]))

    if (alternative=="two.sided"){ 
        P.adjusted[i] <- 1 - pmvt(lower=rep(-abs(Test.Stat[i]),ncomp), upper=rep(abs(Test.Stat[i]),ncomp), df=as.integer(d.freedom),corr=CorrMat.H0, delta=rep(0,ncomp), abseps=1e-05)
        }
    
    if (alternative=="greater"){
        P.adjusted[i] <- 1 - pmvt(lower=rep(-Inf,ncomp), upper=rep(Test.Stat[i],ncomp), df=as.integer(d.freedom),corr=CorrMat.H0, delta=rep(0,ncomp), abseps=1e-05)
        }
    if (alternative=="less"){
        P.adjusted[i] <- 1 - pmvt(lower=rep(-Inf,ncomp), upper=rep(-Test.Stat[i],ncomp), df=as.integer(d.freedom),corr=CorrMat.H0, delta=rep(0,ncomp), abseps=1e-05)
        }
    } # end of for loop

if (alternative=="two.sided"){ 
    Critical.pt <- qmvt(1-FWER, df=as.integer(d.freedom),corr=CorrMat.H0, delta=rep(0,ncomp),tail="both", abseps=1e-05)$quantile
    P.raw  <-  2*pt(abs(Test.Stat),d.freedom,lower.tail=FALSE)
    }
    
if (alternative=="greater"){
     Critical.pt <- qmvt(1-FWER, df=as.integer(d.freedom),corr=CorrMat.H0, delta=rep(0,ncomp),tail="lower.tail", abseps=1e-05)$quantile
     P.raw  <-  pt(Test.Stat,d.freedom,lower.tail=FALSE)
      }
if (alternative=="less"){
     Critical.pt <-  qmvt(1-FWER, df=as.integer(d.freedom),corr=CorrMat.H0, delta=rep(0,ncomp),tail="lower.tail", abseps=1e-05)$quantile
     P.raw  <-  pt(Test.Stat,d.freedom,lower.tail=TRUE)
      }

return(list(
estimate=Ratio.Estimate,
teststat=Test.Stat,
Num.Contrast=Num.Contrast,
Den.Contrast=Den.Contrast,
CorrMat=CorrMat.H0,
critical.pt=Critical.pt,
df=d.freedom,
p.value.raw=P.raw,
p.value.adj=P.adjusted,
Margin.vec=Margin.vec,
alternative=alternative,
FWER=FWER
))

}

