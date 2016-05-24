
simNHP.fun<-function(lambda,fixed.seed=NULL)
{
Tfinal<-length(lambda)
lambdacum<-cumsum(lambda)
lastposH<- lambdacum[Tfinal]
if (!is.null(fixed.seed)) set.seed(fixed.seed)
distexp<-rexp(Tfinal,1)
posHaux<-cumsum(distexp)
posH<-posHaux[posHaux<=lastposH]


posNH<-apply(as.matrix(posH),MARGIN=1,FUN=buscar, lambdacum)

if (length(posNH)>0)
{
 posNH<-posNH[c(1, diff(posNH))!=0]
}
#if the same point is generated twice, one of them is eliminated
#if the number of points is 0 is not necessary to eliminate



return(list(posNH=posNH,lambda=lambda, fixed.seed=fixed.seed))

}