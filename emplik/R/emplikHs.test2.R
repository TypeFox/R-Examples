################################################################
## Poisson EL, with q hazard integration parameters.        ####
## Two samples of right censored, left truncated data       ####
################################################################
emplikHs.test2 <- function(x1, d1, y1= -Inf, x2, d2, y2= -Inf,
       theta, fun1, fun2, maxit = 25, tola = 1e-7, itertrace=FALSE)
{
theta <- as.vector(theta)
q <- length(theta)

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

    newdata1 <- Wdataclean2(z=x1, d=d1)
    temp1 <- DnR(newdata1$value, newdata1$dd, newdata1$weight, y=y1)
    newdata2 <- Wdataclean2(z=x2, d=d2)
    temp2 <- DnR(newdata2$value, newdata2$dd, newdata2$weight, y=y2)

    jump1 <- (temp1$n.event)/temp1$n.risk
    jump2 <- (temp2$n.event)/temp2$n.risk

funtime11 <- as.matrix( fun1(temp1$times) )
  if( ncol(funtime11) != q ) 
     stop("check the output dim of fun1, and theta")
funtime21 <- as.matrix( fun2(temp2$times) )
  if( ncol(funtime21) != q ) 
     stop("check the output dim of fun2, and theta")

Kcent <- jump1%*%funtime11 - jump2%*%funtime21
  if( itertrace ) print(c("Kcenter=", Kcent)) 

    index1 <- (jump1 < 1)
    index2 <- (jump2 < 1)

    K12 <- rep(0, q)  
    tm11 <- temp1$times[!index1]
    if(length(tm11) > 1 ) stop("more than 1 place jump>=1 in x1?")
    if( length(tm11) > 0 ) {
         K12 <- K12 + as.vector(fun1(tm11))
    }
    tm21 <- temp2$times[!index2]
    if(length(tm21) > 1 ) stop("more than 1 place jump>=1 in x2?")
    if( length(tm21) > 0 ) {
         K12 <- K12 - as.vector(fun2(tm21))
    } 

    eve1 <- temp1$n.event[index1]
    tm1 <- temp1$times[index1]
    rsk1 <- temp1$n.risk[index1]
    jmp1 <- jump1[index1]
    funtime1 <- as.matrix(fun1(tm1))
#######################################################
#### it seems I need to include the last point, even it is 1???
########Sample two########

    eve2 <- temp2$n.event[index2]
    tm2 <- temp2$times[index2]
    rsk2 <- temp2$n.risk[index2]
    jmp2 <- jump2[index2]
    funtime2 <- as.matrix(fun2(tm2))

###############################################################
  TINY <- sqrt( .Machine$double.xmin )
  if(tola < TINY) tola <- TINY
  lam <- rep(0,q)
  N <- n1+n2
######################## replace tola if it is too small #####
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

 grad <- gradf3(lam,funtime1,eve1,rsk1,funtime2,eve2,rsk2,K=theta-K12,n=N)
 gsize <- mean( abs(grad) )

 arg1 <- as.vector(rsk1 + funtime1 %*% lam)
 arg2 <- as.vector(rsk2 - funtime2 %*% lam)
 ww1 <- as.vector(-llogpp(arg1, 1/N))^.5
 ww2 <- as.vector(-llogpp(arg2, 1/N))^.5
 tt1 <- sqrt(eve1)*ww1
 tt2 <- sqrt(eve2)*ww2
HESS <- - ( t(funtime1 * tt1)%*%(funtime1 * tt1) +
         t(funtime2 * tt2)%*%(funtime2 * tt2)  )


#                                   -1
#    The Newton step is -(hess'hess)    grad,
#  where the matrix hess is a sqrt of the Hessian.
#  We shall just compute hess'hess = HESS. 
#
 nstep <- as.vector( - solve(HESS, grad) )
  gstep <- grad
  if( sum(nstep^2) < sum(gstep^2) )
    gstep <- gstep*(sum(nstep^2)^.5/sum(gstep^2)^.5)

  ninner <- 0
  for(  i in 1:length(nwts) ){
    lamtemp <- lam+nwts[i]*nstep+gwts[i]*gstep 
   ngrad <- gradf3(lamtemp,funtime1,eve1,rsk1,funtime2,eve2,rsk2,
                                K=theta-K12,n=N)
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

######################################################
lamfun1 <- as.vector(funtime1 %*% lam )
lamfun2 <- as.vector(funtime2 %*% lam )
onePlamh1 <- (rsk1 + lamfun1)/rsk1   ### this is 1 + lam Zi in Ref.
oneMlamh2 <- (rsk2 - lamfun2)/rsk2   ### this is 1 - lam Zi in Ref.

###weights <- jump/onepluslamh  
###need to change last jump to 1? NO. see above

loglik1 <- (sum( eve1*llog(onePlamh1, 1/N)) - 
                 sum(eve1*(lamfun1)/(rsk1 + lamfun1)) )
loglik1fenzi <- - sum( eve1*llog((rsk1+lamfun1), 1/N)) - sum( eve1/onePlamh1 )
loglik2 <- (sum( eve2*llog(oneMlamh2, 1/N)) - 
                 sum(eve2*(-lamfun2)/(rsk2 - lamfun2)) )
loglik <- 2*(loglik1 + loglik2)
#?is that right? YES  see (3.2) in Ref. above. This ALR, or Poisson LR.

list( "-2LLR"=loglik, lambda=lam, "-2LLR(sample1)"=2*loglik1, 
                                  "Llik(sample1)"=loglik1fenzi )
}

gradf3 <- function(lam, funt1, evt1, rsk1, funt2, evt2, rsk2, K, n) {
 arg1 <- as.vector(rsk1 + funt1 %*% lam)
 arg2 <- as.vector(rsk2 - funt2 %*% lam)
VV <- (evt1*llogp(arg1, 1/n))%*%funt1 - (evt2*llogp(arg2,1/n))%*%funt2 - K
return( as.vector( VV ))
}

