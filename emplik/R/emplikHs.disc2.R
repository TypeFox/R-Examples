####################################################
#### 2 sample, discrete hazards, q constraints. #####
#### right censored, left truncated data.      #####
####################################################
emplikHs.disc2 <- function(x1, d1, y1 = -Inf, x2, d2, y2 = -Inf,
         theta, fun1, fun2, maxit = 25, tola = 1e-6, itertrace=FALSE) {
   theta <- as.vector(theta)
   q <- length(theta)
########Sample One########
    x1 <- as.vector(x1)
    n1 <- length(x1)
    if (n1 <= 2*q+1)
        stop("Need more observations in x1")
    if (length(d1) != n1)
        stop("length of x1 and d1 must agree")
    if (any((d1 != 0) & (d1 != 1)))
        stop("d1 must be 0/1's for censor/not-censor")
    if (!is.numeric(x1))
        stop("x1 must be numeric -- observed times")
        
    newdata1 <- Wdataclean2(z=x1, d=d1)
    temp1 <- DnR(newdata1$value, newdata1$dd, newdata1$weight, y=y1)
    jump1 <- (temp1$n.event)/temp1$n.risk
    k1 <- temp1$n.event - temp1$n.risk
    index1 <- (jump1 < 1)
    k1 <- k1[index1]
    eve1 <- temp1$n.event[index1]
    tm1 <- temp1$times[index1]
    rsk1 <- temp1$n.risk[index1]
    jmp1 <- jump1[index1]
    funtime1 <- as.matrix(fun1(tm1))
    if( ncol(funtime1) != q ) stop("check the output dim of fun1")
##    funh1 <- funtime1/rsk1
    
########Sample two########
    x2 <- as.vector(x2) 
    n2 <- length(x2)
    if (n2 <= 2*q+1)
        stop("Need more observations for sample 2")
    if (length(d2) != n2)
        stop("length of x2 and d2 must agree")
    if (any((d2 != 0) & (d2 != 1)))
        stop("d2 must be 0/1's for censor/not-censor")
    if (!is.numeric(x2))
        stop("x2 must be numeric -- observed times")

    newdata2 <- Wdataclean2(z=x2, d=d2)
    temp2 <- DnR(newdata2$value, newdata2$dd, newdata2$weight, y=y2)
    jump2 <- (temp2$n.event)/temp2$n.risk
    k2 <- temp2$n.event - temp2$n.risk
    index2 <- (jump2 < 1)
    k2 <- k2[index2]  
    eve2 <- temp2$n.event[index2]
    tm2 <- temp2$times[index2]
    rsk2 <- temp2$n.risk[index2]
    jmp2 <- jump2[index2]
    funtime2 <- as.matrix(fun2(tm2))
    if( ncol(funtime2) != q ) stop("check the output dim of fun2")
##    funh2 <- funtime2/rsk2
##################################################################
# funtime12 are matrix of n12 x q. rsk12, eve12 are vectors of length n1/n2.
############################################################################
Kcent <- log(1-(eve1/rsk1))%*%funtime1 - log(1-(eve2/rsk2))%*%funtime2 
if( itertrace ) print(c("Kcenter=", Kcent))
##################################################################
  TINY <- sqrt( .Machine$double.xmin )
  if(tola < TINY) tola <- TINY
  lam <- rep(0,q)
  N <- n1+n2
#
#    Preset the weights for combining Newton and gradient
# steps at each of 16 inner iterations, starting with
# the Newton step and progressing towards shorter vectors
# in the gradient direction.  Most commonly only the Newton
# step is actually taken, though occasional step reductions
# do occur.
#

nwts <- c( 3^-c(0:3), rep(0,12) )
gwts <- 2^( -c(0:(length(nwts)-1)))
gwts <- (gwts^2 - nwts^2)^.5
gwts[12:16] <- gwts[12:16] * 10^-c(1:5)

#
#    Iterate, finding the Newton and gradient steps, and
# choosing a step that reduces the objective if possible.
#

nits <- 0
gsize <- tola + 1
while( nits < maxit && gsize > tola ){

  grad <- gradf2(lam, funtime1,eve1,rsk1,funtime2,eve2,rsk2,K=theta,n=N) 
  gsize <- mean( abs(grad) )
## HESS <- hessf2(lam, funtime1, eve1, rsk1, funtime2, eve2, rsk2, n=N)

##hessf2 <- function(lam, funt1, evt1, rsk1, funt2, evt2, rsk2, n) 
 arg1 <- as.vector(rsk1 + funtime1 %*% lam)
 arg2 <- as.vector(rsk2 - funtime2 %*% lam)
 ww1 <- as.vector(-llogpp(arg1, 1/N))^.5
 ww2 <- as.vector(-llogpp(arg2, 1/N))^.5
 tt1 <- sqrt(eve1/(1-eve1/arg1))*ww1  ##
 tt2 <- sqrt(eve2/(1-eve2/arg2))*ww2  ## shall we change to max(TINY,tt2)?
HESS <- t(funtime1 * tt1)%*%(funtime1 * tt1) + 
         t(funtime2 * tt2)%*%(funtime2 * tt2)  

#                                   -1
#    The Newton step is -(hess'hess)    grad,
#  where the matrix hess is a sqrt of the Hessian.
#  We shall just compute hess'hess = HESS. 
#
#####  nstep <- as.vector( - solve(HESS) %*% grad )
################# this may be better #############
 nstep <- as.vector( - solve(HESS, grad) ) 
  gstep <- grad
  if( sum(nstep^2) < sum(gstep^2) )
    gstep <- gstep*(sum(nstep^2)^.5/sum(gstep^2)^.5) 

  ninner <- 0
  for(  i in 1:length(nwts) ){
    lamtemp <- lam+nwts[i]*nstep+gwts[i]*gstep 
    ngrad <- gradf2(lamtemp,funtime1,eve1,rsk1,funtime2,eve2,rsk2,K=theta,n=N)
    ngsize <- mean( abs(ngrad) ) 
    if( ngsize  < gsize  ){
      lam <- lamtemp
      ninner <- i
      break
    }
  }
  nits <- nits+1
  if( ninner==0 )nits <- maxit
  if( itertrace )
    print( c(lam, gsize, ninner) )
}

##################################################################
 arg1 <- as.vector(rsk1 + funtime1 %*% lam)
 arg2 <- as.vector(rsk2 - funtime2 %*% lam)
 onePlam <- arg1/rsk1                   ##########1+lam*funh1
 weights1 <- eve1/arg1                  ##########jmp1/onePlam
 oneMlam <- arg2/rsk2                   ##########1-lam*funh2
 weights2 <- eve2/arg2                  ##########jmp2/oneMlam
 
 loglik1 <- sum(eve1*llog(onePlam, 1/N)) + 
            sum((-k1)*llog((1-jmp1)/(1-weights1), 1/N))
 loglik2 <- sum(eve2*llog(oneMlam, 1/N)) + 
            sum((-k2)*llog((1-jmp2)/(1-weights2), 1/N))
#### Use the llog() instead of log() to avoid infinite, NA, etc.
 loglikR <- 2*(loglik1+loglik2)
#MZ <- gradf2(lam, funtime1, eve1, rsk1, funtime2, eve2, rsk2, K=theta, n=N)
#print(MZ)

list("-2LLR" = loglikR, lambda = lam, "-2LLR(sample1)"=2*loglik1,
      times1 = tm1, times2 = tm2 )
}


####################
gradf2 <- function(lam, funt1, evt1, rsk1, funt2, evt2, rsk2, K, n) {
 arg1 <- as.vector(rsk1 + funt1 %*% lam)
 arg2 <- as.vector(rsk2 - funt2 %*% lam)
VV <- llog(1-(evt1/arg1),1/n)%*%funt1-llog(1-(evt2/arg2),1/n)%*%funt2-K 
return( as.vector( VV ))
}
############################################################################ 
# In the above function, lam, K are vectors of length q. 
# funt12 are matrix of n12 x q. rsk12, evt12 are vectors of length n1/n2. 
############################################################################
