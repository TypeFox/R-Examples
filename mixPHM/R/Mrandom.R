`Mrandom` <-
function(x,shape,scale,prior,K)
{
lik <- NULL
for (k in 1:K) {
  shape.vec <- shape[k,]
  scale.vec <- scale[k,]
  prior.vec <- prior[k,]

  likmat <- t(apply(x,1,function(y){
                    
                    la.vec <- prior.vec*(dweibull(y,shape.vec,scale.vec))        #pages visited by session
                    la.vec[is.na(la.vec)] <- 0
                    la.vec[la.vec==0|la.vec==Inf] <- 1-prior.vec[la.vec==0|la.vec==Inf] #dweibull either 0 or Inf for 0 dwell time (pages not visited by session)
                    return(la.vec)
                    }))
  #likmat <- log(likmat)
  #lik <- cbind(lik,apply(likmat,1,sum))         #multiplying prob over pages
  
  lik <- cbind(lik,apply(likmat,1,prod))         #multiplying prob over pages
}

lik.n <- log(apply(lik, 1, max))                 #maximum log-likelihood value for each session
#lik.n <- apply(lik, 1, max)                 #maximum log-likelihood value for each session
lik.tot <- sum(lik.n)

denom <- apply(lik,1,sum)                       #denominator for posterior
postmat <- apply(lik,2,function(y) {y/denom})   #posterior matrix of dimension n x K

newgr <- apply(postmat,1, function(y) {((1:K)[y==max(y)][1])})   #randomized cluster assignment

list(newgr=newgr,lik.tot=lik.tot,postmat=postmat)
}

