amltest<-function(response, marker, kin, numkeep=floor(length(response)*.5), selectvar){
## This is the main function to fit the adaptive mixed LASSO
## response is the vector of phenotype values. 
## marker is the matrix of genetic effect (could be main and epstatic)
## kin is the relationship matrix between all lines.
## numkeep is the number of genetic effects to keep after the initial screening, default to half of the sample size.
## selectvar is the number of genetic effects to keep after final fit.
## default selectvar=-1 instruct the function to find the number of effects minimizing EBIC
## The value for the function include a list of five items.
##    estimate is matrix with two clumns, the first column indicate which genetic effects (column number of marker)
##      is retained in the fit.  The second column is the effect size.
##   AIC, BIC, EBIC contain the named criterion at each adaptive lasso step up to the last fit.
##   vars contains parameters for the variance components. 
###     vars[1] is the genetic variance component sigma2_{g}
###     vars[2] is the ratio of error variance and genetic variance sigma2_{e}/sigma2_{g}

if (!is.matrix(marker) & !is.data.frame(marker)) stop("marker must be a matrix or data.frame")
if (!is.matrix(kin) & !is.data.frame(kin)) stop("kin must be a matrix or data.frame")
if (is.null(colnames(marker))) stop("Please give unique names to each marker")
if (!is.vector(response) | !is.numeric(response)) stop("response must be a real vector")
if (sum(is.na(response))>0) stop("missing values are not allowed in response")
if (sum(is.na(marker))>0) stop("missing values in the marker matrix, use the cleanclust() function to impute missing values")
if (min(marker)<0 | max(marker)>1) stop("marker should be encoded as 0 or 1 according to the presence of the minor alleles")



  numline<-dim(marker)[1]
  nummarker<-dim(marker)[2]

if (numline!=dim(kin)[1] | numline!=dim(kin)[2]) stop("kin must have the same number of rows and columns as the number of lines")
if (numkeep<= 0 | numkeep>=numline) stop("numkeep must be a positive integer less than the number of lines")
if (selectvar<0 | selectvar>=numkeep) stop("selectvar must be a nonnegative integer less than numkeep")

markmean<- apply(marker,2,mean) 
if (sum(markmean>.5)>0) stop("the presence of the minor allele should be encoded as 1 for all markers, you can use the cleanclust() function to re-encode the markers.")

  epmark<-sweep(marker, 2, markmean,"-")
  epweight<-sqrt(apply(epmark^2, 2, sum))

  epmark2<-t(t(epmark)/epweight)

  respon<-response-mean(response)

  enk<-svd(kin)
  tem2p<-lars(cbind(enk$u[,1:3],epmark2), respon, normalize=F, max.steps=numkeep, use.Gram=F)

  thstart<-c(0.2*var(respon),1)

  resida<-t(respon)%*%enk$u
  tem0<-   optim(thstart, .lik2, .grlik2, method="L-BFGS-B", lower=c(var(respon)*.001, .001), upper=c(10*var(respon), 1000), residu=resida, hd=enk$d)
  y<-respon


  if (selectvar==0){
   estimate<- NULL
   aic<-tem0$value+ 2
   bic<-tem0$value+log(numline)
   ebic<-bic
   vars<- tem0$par
   res<-list(estimate=estimate, AIC=aic, EBIC=ebic, BIC=bic, vars=vars) 
   return(res)
  } else{
    exclud<-abs(tem2p$beta[numkeep,4:(nummarker+3)])
    markweight<-epmark2[, exclud>0]
    temnum<-sum(exclud>0)
    ebic<-tem0$value+log(numline)+2*log(max(1, .nchoose2(temnum, 1)))
    bic<-tem0$value+log(numline)
    aic<-tem0$value+ 2
    inweight<-abs(lm(y~markweight)$coef[2:(temnum+1)])
    minweight<-matrix(rep(inweight,numline),numline,temnum, byrow=T)
    hd<-enk$d
    tem0a<-1:nummarker
    temseq<-tem0a[exclud>0]
    beta00<-matrix(0,temnum*3,nummarker)
    smhaf<-enk$u%*%(diag(1/sqrt(hd+tem0$par[2]))/sqrt(tem0$par[1]))%*%t(enk$u)
    theta1<-tem0$par
    dfvec<-NULL
    varcom<-NULL
    mcount<- 0
    for (i in 1:(selectvar)){
       newy<-smhaf%*%y
       xnew<-minweight*(smhaf%*%markweight)
       tem2<-lars(xnew, newy, normalize=F, max.steps=i, use.Gram=T)
       numrange<-length(tem2$df)
       beta00[i,temseq]<-tem2$beta[numrange,]*inweight
       resid<-t(respon-markweight%*%(tem2$beta[numrange,]*inweight))%*%enk$u
       tem1<-   optim(thstart, .lik2, .grlik2, method="L-BFGS-B", lower=c(.01, .001), upper=c(1000, 1000), residu=resid, hd=enk$d)
       smhaf<-enk$u%*%(diag(1/sqrt(hd+tem1$par[2]))/sqrt(tem1$par[1]))%*%t(enk$u)
       thstart<-tem1$par
       varcom<-rbind(varcom, t(thstart))
       mcount<- c(mcount, tem2$df[i+1]-1)
       dfvec<-c(dfvec, tem2$df[numrange])
       ebic<-c(ebic, tem1$value+tem2$df[numrange]*log(numline)+2*log(max(1, .nchoose2(temnum, tem2$df[numrange])))  )
       bic<-c(bic, tem1$value+tem2$df[numrange]*log(numline) )
       aic<-c(aic, tem1$value+tem2$df[numrange]*2)
       ### cat( i, "\n")
    }
  
 

  
   vars<- varcom[selectvar,]
   bicset<-as.vector(which(beta00[selectvar,]!=0))
   tem3<-beta00[,bicset]
   if (is.vector(tem3)) {
        enter<- .calseq(tem3)
   }else{
        enter<-apply(tem3, 2, .calseq)
   }
   oenter<-order(enter,decreasing=T)
   markernm<-bicset[oenter]

   estim<-beta00[selectvar,]/epweight
   estimate<-cbind(markernm, size=estim[markernm])

### the first element of BIC/EBIC/AIC are that with zero predictors.
   res<-list(estimate=estimate, AIC=aic, BIC=bic, EBIC=ebic, vars=vars,mcount=mcount) 
    return(res)
 }
}


