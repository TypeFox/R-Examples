`Emax` <-
function(x, old, K, method, Sdist,p, cutpoint)
#old...matrix with probabilites for group membership

{

shape <- matrix(NA, K, p)	           # K x p Matrix with Shape-Parameter
scale <- matrix(NA, K, p)            # K x p Matrix with Scale-Parameter
n <- nrow(x)

prior.vec <- apply(x,2,function(y){ ly <- length(y[y>0])/n})   #prob that certain pages not visited by session
prior <- matrix(prior.vec,nrow=K,ncol=p,byrow=TRUE)            #matrix of prior probabilities 

old[old==0] <- min(old[old>0])      #0-weights in survreg not allowed (replaced with minimum)

if (method=="separate") {
  parlist <- apply(old,2,function(y) {                   #old is matrix with posteriors
                apply(x,2,function(z) {
                 censvec <- rep(1, length(z))   
                 censvec[z > cutpoint] <- 0     #vector for censored data (set to 0)          
                 wphm <- survreg(Surv(z[z>0], censvec[z>0])~1, weights = y[z>0], dist = Sdist)
                 shapep <- 1/wphm$scale
                 scalep <- exp(wphm$coefficients[1])
                 list(scalep,shapep)
                 })})
  shsclist <- tapply(unlist(parlist),rep(1:2,length(unlist(parlist))/2),function(z){
                                          matrix(z,nrow=K,byrow=TRUE)})                         #reorganizing parlist
  shape <- shsclist[[2]]                                                      #shape matrix K,p
  scale <- shsclist[[1]]
  anzpar <- 2*K*p
}


list (scale=scale,shape=shape,prior=prior,anzpar=anzpar)
}

