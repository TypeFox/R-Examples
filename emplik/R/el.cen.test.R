el.cen.test <- function(x,d,fun=function(x){x},mu,error=1e-8,maxit=15)
{
   xvec <- as.vector(x)
   n <- length(xvec)
   if(n <= 2) stop ("Need more observation")
   if(length(d)!=n) stop("length of x and d must agree")
   if(any((d!=0)&(d!=1)))
      stop ("d must be 0(right-censored) or 1(uncensored)")
   if(!is.numeric(xvec)) stop("x must be numeric")
   if(length(mu)!=1) stop("check the dim of constraint mu")

   temp <- Wdataclean2(xvec,d)
   dd <- temp$dd
   dd[length(dd)] <- 1

   if(all(dd==1)) stop("there is no censoring, please use el.test()")

   xx <- temp$value
   n <- length(xx) 
   ww <- temp$weight
   w0 <- WKM(x=xx, d=dd, w=ww)$jump
   uncenw0 <- w0[dd==1]
   funxx <- fun(xx)

   if((mu>max(funxx))|(mu<min(funxx))) stop("check the value of mu/fun")

   xbar <- sum(funxx[dd==1] * uncenw0)

    #********* begin initial calculation******************
    # get vector dvec which is the first derivative vector

     dvec01 <- uncenw0
     rk <- 1:n            ######## rank(sortx) = 1:n  yes!
     cenrk <- rk[dd==0]
     mm <- length(cenrk)
     dvec02 <- rep(0,mm)
     for(j in 1:mm)  dvec02[j] <- sum(w0[cenrk[j]:n])
     dvec00 <- rep(0,n)
     dvec00[dd==1] <- dvec01
     dvec00[dd==0] <- dvec02
     dvec0 <- ww/dvec00

     # get matix Dmat which is Decompition of 2nd derivative matrix.
     # Dmat0 <- diag(1/dvec0)
     Dmat0 <- dvec00/sqrt(ww)

     # get constraint matrix Amat
     mat <- matrix(rep(dd,mm),ncol=mm, nrow=n)
     for(i in 1:mm)
     {
        mat[1:cenrk[i],i] <- 0
        mat[cenrk[i],i] <- -1
     }

    Amat <- as.matrix(cbind(dd, funxx * dd, mat))

     # get constraint vector bvec
     bvec0 <- c(0,as.vector(mu-xbar),rep(0,mm))

     # Use solve3.QP to maximize the loglikelihood function
     value0<-solve3.QP(Dmat0,dvec0,Amat,bvec0,meq=mm+2,factorized=TRUE)

     w <- dvec00 + value0$solution

     if(any(w<=0))
     stop("There is no probability satisfying the constraints")

     #**********end initial calculation **********************
     #**********begin iteration ******************************
     # update vector Dmat, dvec and bvec after initial calculation
     bvec <- rep(0,mm+2)
     diff <- 10
     m <- 0
     while( (diff>error) & (m<maxit) )
      {
         dvec <- ww/w
         # get matix Dmat
         #  Dmat <- diag(w)
         Dmat <- w/sqrt(ww)
         value0 <- solve3.QP(Dmat,dvec,Amat,bvec,meq=mm+2,factorized=TRUE)
         w <- w + value0$solution
         #diff <- sum(value0$solution^2)
         diff <- sum(abs( value0$solution) )
         m <- m+1
       }
     #**********end iteration ******************************

     lik00 <- sum(ww*log(dvec00))

     tval <- 2*(lik00 - sum(ww*log(w)))
   list("-2LLR"=tval, Pval=1-pchisq(tval, df=1),
            weights=w[dd==1], xtime=xx[dd==1], iteration=m, error=diff)
}

