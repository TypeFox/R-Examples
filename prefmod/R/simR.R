# simulates rankings data with given worth
#
simR<-function(nobj, nobs, worth=NULL, seed=NULL, pr=FALSE)
{
   if(nobj>10)
     stop("nobj too large!")
   if(!is.null(seed)) set.seed(seed)

   if(is.null(worth)){
       worth<-runif(nobj)
   } else if(length(worth)!=nobj)
       stop("length of worth is not equal to nobj\n")
   worth<-worth/sum(worth)
   if (pr)
     cat("used worth parameters are: ", worth, "\n")

   P<-permutations(nobj)

   R<-nobj-P              # following Critchlow & Fligner
   patt<-R%*%log(worth)

   p<-exp(patt)/sum(exp(patt))


   vecF<-rmultinom(1,nobs, p)
   simdat<-expand.mat(P,vecF)
   data.frame(simdat)
}
