`lmec` <-
function(yL, cens, X, Z, cluster, maxstep = 200, varstruct = "unstructured", 
	     init, method="ML", epsstop=1e-3, abspmv=1e-3, mcmc0=100, sdl=0.1, iter2 = 15, trs=5, pls=5, mcmcmax=1000)
{
library(mvtnorm)

tr <- function(A)
{  if(nrow(A)==1) a <- as.numeric(A)
  else  a <- sum(diag(A))
  return(a)
}

bdiag <- function(x){ 
     if(!is.list(x)) stop("x not a list") 
     n <- length(x) 
     if(n==0) return(NULL) 
     x <- lapply(x, function(y) if(length(y)) as.matrix(y) else 
stop("Zero-length component in x")) 
     d <- array(unlist(lapply(x, dim)), c(2, n)) 
     rr <- d[1,] 
     cc <- d[2,] 
     rsum <- sum(rr) 
     csum <- sum(cc) 
     out <- array(0, c(rsum, csum)) 
     ind <- array(0, c(4, n)) 
     rcum <- cumsum(rr) 
     ccum <- cumsum(cc) 
     ind[1,-1] <- rcum[-n] 
     ind[2,] <- rcum 
     ind[3,-1] <- ccum[-n] 
     ind[4,] <- ccum 
     imat <- array(1:(rsum * csum), c(rsum, csum)) 
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2], 
(y[3]+1):y[4]], imat=imat) 
     iuse <- as.vector(unlist(iuse)) 
     out[iuse] <- unlist(x) 
     return(out) 
} 


p <- ncol(X);
q <- ncol(Z)
m <- length(unique(cluster))
n <- length(cluster)
ni <- table(cluster)
cl <- cluster

lflag <- 0 # likelihood flag; 0 = inactive; 1 = increasing lik; -1 = decreasing

# Initial values (if not provided)
fit <- lm( yL ~ -1 + X )
if (missing(init) || is.na(match("beta", names(init) ) ) )
{
  beta <- fit$coef;   names(beta) <- NULL
}
else {
  beta <- init$beta; names(beta) = NULL
}

if (missing(init) || is.na(match("sigma", names(init) ) ) )
  sigma2 <- (summary(fit)$sigma)^2
else
  sigma2 <- init$sigma^2

if (missing(init) || is.na(match("Delta", names(init) ) ) )
  Delta <- diag(rep(1,q))


else
  Delta = init$Delta

if (missing(init) || is.na(match("bi", names(init) ) ) )
  b <- matrix(rep(0,m*q), ncol=m) 
else
  b = init$bi


old.lik <- 0 # lik = objective function, not true likelihood
likseq <- rep(0,maxstep)
likseqex <- rep(0,maxstep)

A <- matrix(nrow=n, ncol=p)
a <- rep(0,n)
y = yL	# y stores y for observed y and current Ey for censored y

if ( method == "ML" )
 {
   for( iter in ( 1:maxstep ) )
     {
       lik <- s2 <- 0
       likex <- 0
       Psi <- Psi2 <- diag(0,q)
       for( i in 1:m )
         {
           Zi <- Z[cl==i, , drop=FALSE]
           Ziaug <- rbind(Zi,Delta)
           qrZ <- qr(Ziaug)
           Qi <- qr.Q( qrZ, compl=TRUE)[,,drop=FALSE]   # 12x12
           QTi <- qr.Q( qrZ, compl=TRUE)[1:ni[i],,drop=FALSE]
           QLTi <- QTi[,1:q, drop=FALSE]
           QRTi <- QTi[,-(1:q), drop=FALSE]
           R11i <- qr.R( qrZ )
           R00i <- t(QRTi) %*% X[cl==i,, drop=FALSE]
           A[cl==i,] <- R00i
           wi <- t(QLTi) %*% X[cl==i,, drop=FALSE] %*% beta
           R11iinv <- solve(R11i)
           nic <- sum( (cl==i) & (cens==1) )
           censi <- cens[cl==i]
           indi <- cbind(1:ni[i], censi)[, 1]

           if (nic==0)
             {
               Ui <- diag(0,q)
               yo <- yL[(cl==i)]
               uo <- X[cl==i,,drop=FALSE]%*%beta
               invDelta <- solve(t(Delta)%*%Delta)
               So <- (sigma2*(diag(1,ni[i])+Zi%*%invDelta%*%t(Zi)))
             }
           if (nic>0)
             {
               Xc <- X[cl==i & cens==1, ,drop=FALSE]
               Xo <- X[cl==i & cens==0, ,drop=FALSE]
               Zc <- Z[cl==i & cens==1, ,drop=FALSE]
               Zo <- Z[cl==i & cens==0, ,drop=FALSE]
               qc <- yL[(cl==i) & cens==1]
               yo <- yL[(cl==i) & cens==0]
               indc<- indi[censi==1]
               indo <- indi[censi==0]
               QLic <- Qi[-indo,1:q,drop=FALSE]   # censored rows plus zero rows
               QLio <- Qi[-indc,1:q,drop=FALSE]

               QLTic <- QLTi[indc, ,drop=FALSE]   # only censored rows
               QLTio <- QLTi[indo, ,drop=FALSE]

               invQLio2 <- solve(t(QLio)%*%QLio)
               u <- Xc%*%beta+{QLTic%*%invQLio2%*%t(QLTio)}%*%(yo-Xo%*%beta)
               S <- sigma2*{diag(1,nic)+QLTic%*%invQLio2%*%t(QLTic)}

               # mean and variance for yo
               uo <- Xo%*%beta
               invDelta <- solve(t(Delta)%*%Delta)
               So <- (sigma2*(diag(1,ni[i])+Zi%*%invDelta%*%t(Zi)))[indo, indo]
               
               if (nic==1)
                 { 
                   qq <- (1/sqrt(S))*(-qc+u)
                   R<-1
                   alpha <- pnorm(-qq)
                   dd <- dnorm(-qq)
                   H <- qq*dd
                   EX <- (1/alpha)*dd   # a vector with a length of nic
                   EXX <- 1+1/alpha*H
                   varX <- EXX-EX^2
                   Eycens <- -sqrt(S)*EX+u
                   varyic<- varX*S
                 }

               else
                 {   
                   qq <- diag(1/sqrt(diag(S)))%*%(-qc+u)
                   R <-  diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
                   alpha <- pmvnorm(upper=as.vector(-qq), corr=R, abseps=abspmv)
                   dd <- rep(0, nic)   #derivative vector 
                   for (j in 1:nic)
                   {                V <- R[-j, -j, drop=FALSE]-R[-j,j, drop=FALSE]%*%R[j,-j, drop=FALSE]
                                    nu <- -qq[-j]+R[-j,j, drop=FALSE]%*%qq[j]
                                    dd[j] <- dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V, abseps=abspmv)
                   }

                   H <- matrix(rep(0, nic*nic), nrow=nic)
                   RH <- matrix(rep(0, nic*nic), nrow=nic)
                   if(nic==2)
                   {
                     H[1,2] <- H[2,1] <- dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
                     #sigma==R since qq is standardized

                     RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
                   }
                   else
                     {

                       for( s in 1:(nic-1))
                       {
                         for (t in (s+1):nic)
                         {
                         invR <- solve(R[c(s,t), c(s,t), drop=FALSE])
                         nu <- -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=FALSE]%*%invR%*%qq[c(s,t),,drop=FALSE]                      
                         V <-  R[-c(s,t), -c(s,t), drop=FALSE]- R[-c(s,t), c(s,t), drop=FALSE]%*%invR%*%R[c(s,t), -c(s,t), drop=FALSE]
                         H[s,t] <- H[t,s] <- pmvnorm(upper=as.vector(nu), sigma=V, abseps=abspmv)*dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
                         RH[s,t] <- RH[t,s] <- R[s,t]*H[s,t]
                         }
                       }
                   }
                   h <- qq*dd-apply(RH, 1, sum)
                   diag(H) <- h               
                   EX <- (1/alpha)*R%*%dd   # a vector with a length of nic
                   EXX <- R+1/alpha*R%*%H%*%R
                   varX <- EXX-EX%*%t(EX)
                   Eycens <- -diag(sqrt(diag(S)))%*%EX+u
                   varyic <- diag(sqrt(diag(S)))%*%varX%*%diag(sqrt(diag(S)))
                }

                   trvaryic <- tr(varyic)
                   Ui <- t(QLTic) %*% varyic %*% QLTic
                   s2 <- s2 + trvaryic/n - tr(Ui)/n
                   y[ (cl==i) & (cens==1) ] <- Eycens
             } # end if (nic>0)
               c0i <- t(QRTi) %*% y[cl==i]
               a[cl==i] <- c0i[,1]
               c1i <- t(QLTi) %*% y[cl==i]
               bi = R11iinv %*% (c1i- wi) 
               b[,i] <- bi
               Wi = R11iinv %*% t(R11iinv)
               Psi <- Psi + 1/m * R11iinv %*% ( (c1i-wi) %*% t(c1i-wi) + Ui ) %*% t(R11iinv)
               Psi2 <- Psi2 + 1/m * Wi
               lik <- lik - sum(log(abs(diag(R11i))))

            if (nic==0)        likex <- likex+log(dmvnorm(yo, mean=as.vector(uo), sigma=So))                            
            else
               {
                 if(nic==ni[i])
                    likex <- likex+log(alpha)
                 else
                 {  if(ni[i]-nic==1)
                      likex <- likex+log(alpha)+log(dnorm(yo, mean=uo, sd=sqrt(So)))
                    else
                      likex <- likex+log(alpha)+log(dmvnorm(yo, mean=as.vector(uo), sigma=So))  
                 }
                }
              
         
               if (iter>10) # Compute varFix - if last iteration
                 {
                   QRTic <- QRTi[cens[cl==i]==1, ,drop=FALSE]
                   if (nic==0) Vi <- diag(ni[i])-diag(ni[i])
                   else  Vi <- t(QRTic) %*% varyic %*% QRTic 
                   if(i==1) 
                     {
                       inv.varbeta1 = t(R00i) %*% R00i
                       inv.varbeta2 = t(R00i) %*% Vi %*% R00i
                     }
                   else 
                     {
                       inv.varbeta1 = inv.varbeta1 + t(R00i) %*% R00i
                       inv.varbeta2 = inv.varbeta2 + t(R00i) %*% Vi %*% R00i
                     }
                 }
               
        } #end for i

   qra <- qr(A)
   cc = t(qr.Q(qra, compl=TRUE)) %*% a
   beta <- solve( qr.R(qra), cc[1:p,1] )   
   sigma2 <- s2 + sum( cc[(p+1):n]^2 )/n
   if (!is.matrix(Delta)) Delta <- matrix(Delta)
   lik <- lik -n/2*(1 + log(2*pi)) - n/2* log(sigma2) + m * sum(log(diag(Delta)))
   Psi <- Psi + Psi2 * sigma2   
   if (varstruct=="diagonal" & (q>1)) Psi = diag(diag(Psi))

   Delta <- chol(solve(Psi/sigma2))

   likseq[iter] <- lik
   likseqex[iter] <- likex

   diff<-likseqex[iter]-likseqex[iter-1]
   if(iter>10&&diff<epsstop)     
      if(lflag==0) lflag<-1
	  else lflag<-2
   print(paste(c(iter, lflag, diff, likex), sep=" "))      
       if (iter==maxstep|lflag==2)
    {  varFix = sigma2 * solve(inv.varbeta1 - inv.varbeta2/sigma2)
#      vF = sigma2 * solve(inv.varbeta1)
      break
    }       
  } #end for iter
return(list(beta=beta, bi=b, sigma=sqrt(sigma2), Psi=Psi, Delta=Delta, loglik = lik, 
            loglikex=likex, varFix = varFix, method=method, varstruct=varstruct, 
			step=iter, likseq=likseq[1:iter], likseqex=likseqex[1:iter]))
 }

 
if ( method == "REML" )
 {
   for( iter in ( 1:maxstep ) )
     { 
       lik <- s2 <- 0
       Psi <- Psi2 <- diag(0,q)
	   
       Rbinv <- list()
       Rbeta <- NULL
	   QRT <- list()
	   QLT <- list()
       Eycens <- NULL
       varyic <- list()
	   nocenclus <- NULL
	   
       for( i in 1:m )
         {
           Zi <- Z[cl==i, , drop=FALSE]
           Ziaug <- rbind(Zi,Delta)
           qrZ <- qr(Ziaug)
           Qi <- qr.Q( qrZ, compl=TRUE)[,,drop=FALSE]   # 12x12
           QTi <- qr.Q( qrZ, compl=TRUE)[1:ni[i],,drop=FALSE]
           QLTi <- QTi[,1:q, drop=FALSE]
           QRTi <- QTi[,-(1:q), drop=FALSE]
           R11i <- qr.R( qrZ )
           R00i <- t(QRTi) %*% X[cl==i,, drop=FALSE]
           A[cl==i,] <- R00i
           wi <- t(QLTi) %*% X[cl==i,, drop=FALSE] %*% beta
           R11iinv <- solve(R11i)
		   
		   Rbinv[[i]] <- R11iinv
		   Rbeta <- rbind(Rbeta, t(QLTi) %*% X[cl==i,, drop=FALSE])
		   QRT[[i]] <- QRTi
		   QLT[[i]] <- QLTi

           nic <- sum( (cl==i) & (cens==1) )
           censi <- cens[cl==i]
           indi <- cbind(1:ni[i], censi)[, 1]

######## Calculate Ey and Vary

           if (nic==0)
             {
               Ui <- diag(0,q)
               yo <- yL[(cl==i)]
               uo <- X[cl==i,,drop=FALSE]%*%beta
               invDelta <- solve(t(Delta)%*%Delta)
               So <- (sigma2*(diag(1,ni[i])+Zi%*%invDelta%*%t(Zi)))
               nocenclus <- c(nocenclus, i)       
             }


           if (nic>0)
             {
               Xc <- X[cl==i & cens==1, ,drop=FALSE]
               Xo <- X[cl==i & cens==0, ,drop=FALSE]
               Zc <- Z[cl==i & cens==1, ,drop=FALSE]
               Zo <- Z[cl==i & cens==0, ,drop=FALSE]
               qc <- yL[(cl==i) & cens==1]
               yo <- yL[(cl==i) & cens==0]
               indc<- indi[censi==1]
               indo <- indi[censi==0]
               QLic <- Qi[-indo,1:q,drop=FALSE]   # censored rows plus zero rows
               QLio <- Qi[-indc,1:q,drop=FALSE]


               QLTic <- QLTi[indc, ,drop=FALSE]   # only censored rows
               QLTio <- QLTi[indo, ,drop=FALSE]
			                  
               invQLio2 <- solve(t(QLio)%*%QLio)
               u <- Xc%*%beta+{QLTic%*%invQLio2%*%t(QLTio)}%*%(yo-Xo%*%beta)
               S <- sigma2*{diag(1,nic)+QLTic%*%invQLio2%*%t(QLTic)}

               # mean and variance for yo
               uo <- Xo%*%beta
               invDelta <- solve(t(Delta)%*%Delta)
               So <- (sigma2*(diag(1,ni[i])+Zi%*%invDelta%*%t(Zi)))[indo, indo]
               
               if (nic==1)
                 { 
                   qq <- (1/sqrt(S))*(-qc+u)
                   R<-1
                   alpha <- pnorm(-qq)

                   dd <- dnorm(-qq)
                   H <- qq*dd
                   EX <- (1/alpha)*dd   # a vector with a length of nic
                   EXX <- 1+1/alpha*H
                   varX <- EXX-EX^2
				   Eycensi <- -sqrt(S)*EX+u
                   Eycens <-c(Eycens, Eycensi)
                   varyic[[i]]<- varX*S
                 }

               else
                 {   
                   qq <- diag(1/sqrt(diag(S)))%*%(-qc+u)
                   R <-  diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
                   alpha <- pmvnorm(upper=as.vector(-qq), corr=R, abseps=abspmv)
                   dd <- rep(0, nic)   #derivative vector 
                   for (j in 1:nic)
                   {                V <- R[-j, -j, drop=FALSE]-R[-j,j, drop=FALSE]%*%R[j,-j, drop=FALSE]
                                    nu <- -qq[-j]+R[-j,j, drop=FALSE]%*%qq[j]
                                    dd[j] <- dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V, abseps=abspmv)
                   }

                   H <- matrix(rep(0, nic*nic), nrow=nic)
                   RH <- matrix(rep(0, nic*nic), nrow=nic)
                   if(nic==2)
                   {
                     H[1,2] <- H[2,1] <- dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
                     #sigma==R since qq is standardized

                     RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
                   }
                   else
                     {

                       for( s in 1:(nic-1))
                       {
                         for (t in (s+1):nic)
                         {
                         invR <- solve(R[c(s,t), c(s,t), drop=FALSE])
                         nu <- -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=FALSE]%*%invR%*%qq[c(s,t),,drop=FALSE]                      
                         V <-  R[-c(s,t), -c(s,t), drop=FALSE]- R[-c(s,t), c(s,t), drop=FALSE]%*%invR%*%R[c(s,t), -c(s,t), drop=FALSE]
                         H[s,t] <- H[t,s] <- pmvnorm(upper=as.vector(nu), sigma=V, abseps=abspmv)*dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
                         RH[s,t] <- RH[t,s] <- R[s,t]*H[s,t]
                         }
                       }
                   }
                   h <- qq*dd-apply(RH, 1, sum)
                   diag(H) <- h               
                   EX <- (1/alpha)*R%*%dd   # a vector with a length of nic
                   EXX <- R+1/alpha*R%*%H%*%R
                   varX <- EXX-EX%*%t(EX)
				   Eycensi=-diag(sqrt(diag(S)))%*%EX+u
                   Eycens <- c(Eycens, Eycensi)	
                   varyic[[i]] <- diag(sqrt(diag(S)))%*%varX%*%diag(sqrt(diag(S)))
                }

                   y[ (cl==i) & (cens==1) ] <- Eycensi
             } # end if (nic>0)
			 
               c0i <- t(QRTi) %*% y[cl==i]
               a[cl==i] <- c0i[,1]
               c1i <- t(QLTi) %*% y[cl==i]
               bi = R11iinv %*% (c1i- wi) 
               b[,i] <- bi
               Psi <- Psi + 1/m *bi%*%t(bi) 

        } #end for i

   if(length(nocenclus)>0)	
       varyc <- bdiag(varyic[-nocenclus])
   if (length(nocenclus)==0) 
	   varyc <- bdiag(varyic)
   
   Rbinv <- bdiag(Rbinv)
   QRT <- bdiag(QRT)
   QLT <- bdiag(QLT)
 
   QRTc <- QRT[cens==1,]
   QLTc <- QLT[cens==1,]
   
   qra <- qr(A)
   cc = t(qr.Q(qra, compl=TRUE)) %*% a
   beta <- solve( qr.R(qra), cc[1:p,1] )  
   
   Q0L <- qr.Q(qra)
   Q0Lc <- Q0L[cens==1,]
   QMT <-QRT%*%Q0L
   QMTc <- QMT[cens==1,]
   Q3LTc <- cbind(QLT[cens==1,], QMTc)
   B <- t(Q3LTc)%*%varyc%*%Q3LTc
   R00inv <- solve(qr.R(qra))
   Rinv1 <- cbind(Rbinv, -Rbinv%*%Rbeta%*%R00inv)
   Rinv2 <- cbind(matrix(rep(0, p*m*q), nrow=p), R00inv)
   Rinv <- rbind(Rinv1, Rinv2)
   
   sigma2 <- (sum( cc[(p+1):n]^2) + tr(varyc)-tr(B))/(n-p)
   vardelta <- sigma2*Rinv%*%t(Rinv)+Rinv%*%B%*%t(Rinv)
   varb <- vardelta[1:(m*q),1:(m*q)]
   sumvarb <- 0
   for ( i in 1:m)
   sumvarb <- sumvarb+varb[((i-1)*q+1):(i*q),((i-1)*q+1):(i*q)]
   
   Psi <- Psi + sumvarb/m
   
   
   if (!is.matrix(Delta)) Delta <- matrix(Delta)
   lik <- sum(log(abs(diag(Rinv)))) -(n-p)/2*(1 + log(2*pi)) - (n-p)/2* log(sigma2) + m * sum(log(diag(Delta)))

   if (varstruct=="diagonal" & (q>1)) Psi = diag(diag(Psi))

   Delta <- chol(solve(Psi/sigma2))
   likseq[iter] <- lik

       diff<-likseq[iter]-likseq[iter-1]
       if(iter>10&&diff<epsstop)     
         if(lflag==0) lflag<-1
         else lflag<-2
  print(paste(c(iter, lflag, diff, lik), sep=" "))      
       if (iter==maxstep|lflag==2)
    {  
	varFix <- vardelta[(m*q+1):(m*q+p),(m*q+1):(m*q+p)] 
    #  vF = sigma2 * solve(inv.varbeta1)
      break
    }       
  } #end for iter
return(list(beta=beta, bi=b, sigma=sqrt(sigma2), Psi=Psi, Delta=Delta, loglik = lik, 
	    varFix = varFix,  method=method, varstruct=varstruct, step=iter, likseq=likseq[1:iter]))
 }


if ( method == "MLmcmc" )
{
vtrnorm2 = function(q)
# Given vector q, sample y from std norm, RIGHT-truncated at q
# Based on trnorm2
{
 u = runif(length(q))
 return( qnorm(log(u) + pnorm(q, log.p=TRUE), log.p=TRUE) )
}

 
 ysample = y# holds current MCMC value
 # y holds the current average y value
 mcemstage <- 1 # stage of the MCEM algorithm
 G = mcmc0  # initialize MCMC sample size
 startstage = 1     # step at which current stage started

  for( iter in ( 1:maxstep ) )
  {
   lik <- s2 <- 0
   Psi <- Psi2 <- diag(0,q)
   for( i in 1:m )
   {
    Zi <- Z[cl==i,]
    if(is.vector(Zi)) Zi = matrix(Zi)
    Ziaug <- rbind(Zi,Delta)
    qrZ <- qr(Ziaug)
    Qi <- qr.Q( qrZ, compl=TRUE)[1:ni[i],]
    QLi <- Qi[,1:q]
    QRi <- Qi[,-(1:q)]
    R11i <- qr.R( qrZ )
    R00i <- t(QRi) %*% X[cl==i,]
    if (is.vector(R00i)) R00i = matrix(R00i)
    A[cl==i,] <- R00i
    wi <- t(QLi) %*% X[cl==i,] %*% beta
    R11iinv <- solve(R11i)
    nic <- sum( (cl==i) & (cens==1) )
    if (nic==0)
    { Ui = diag(0,q) }
    else
    { 
     # MCMC step: 
     Eycens = rep(0,nic)
     Eyycens = diag(0,nic) 
     for (g in 1:G)
     { 
      c1i <- t(QLi) %*% ysample[cl==i]
      vi <- rnorm(q)
      bi <- R11iinv %*% (c1i - wi + sqrt(sigma2) * vi)
      mui <- X[(cl==i) & (cens==1),] %*% beta + Z[(cl==i)&(cens==1),] %*% bi
      qqi <- (yL[(cl==i) & (cens==1)] - mui) / sqrt(sigma2)
      uic <- vtrnorm2(qqi)
      ycens <- mui + sqrt(sigma2)*uic
      Eycens = Eycens + ycens/G
      Eyycens = Eyycens + 1/G * ycens %*% t(ycens)
      ysample[(cl==i)& (cens==1)] <- ycens
     } # end for g
     varyic <- Eyycens - Eycens %*% t(Eycens)
     QLic <- QLi[cens[cl==i]==1,]
     if((q>1) & nic==1) QLic = t(QLic)
     y[ (cl==i) & (cens==1) ] <- Eycens
     trvaryic <- tr(varyic)
     Ui <- t(QLic) %*% varyic %*% QLic
     s2 <- s2 + trvaryic/n - tr(Ui)/n
    }    
    c0i <- t(QRi) %*% y[cl==i]
    a[cl==i] <- c0i[,1]
    c1i <- t(QLi) %*% y[cl==i]
    bi = R11iinv %*% (c1i- wi) 
    b[,i] <- bi
    Wi = R11iinv %*% t(R11iinv)
    Psi <- Psi + 1/m * R11iinv %*% ( (c1i-wi) %*% t(c1i-wi) + Ui ) %*% t(R11iinv)
    Psi2 <- Psi2 + 1/m * Wi
    lik <- lik - sum(log(abs(diag(R11i))))
    # Compute varFix - if last iteration
    if ((iter==maxstep) | ((mcemstage==4) & iter==startstage+pls))
    {
      QRic <- QRi[cens[cl==i]==1,]
      if((q>1) & nic==1) QRic = t(QRic)
      if (nic==0) Vi <- diag(rep(0, sum(cluster==i)))
      else  Vi <- t(QRic) %*% varyic %*% QRic 
      if(i==1) 
      {
       inv.varbeta1 = t(R00i) %*% R00i
       inv.varbeta2 = t(R00i) %*% Vi %*% R00i
      }
      else 
      {
       inv.varbeta1 = inv.varbeta1 + t(R00i) %*% R00i
       inv.varbeta2 = inv.varbeta2 + t(R00i) %*% Vi %*% R00i
      }
    }
   }#end for i
   qra <- qr(A)
   cc = t(qr.Q(qra, compl=TRUE)) %*% a
   beta <- solve( qr.R(qra), cc[1:p,1] )
   sigma2 <- s2 + sum( cc[(p+1):n]^2 )/n 
   if (!is.matrix(Delta)) Delta <- matrix(Delta)
   lik <- lik -n/2*(1 + log(2*pi)) - n/2* log(sigma2) + m * sum(log(diag(Delta)))
   Psi <- Psi + Psi2 * sigma2
   if (varstruct=="diagonal" & (q>1)) Psi = diag(diag(Psi))

   Delta <- chol(solve(Psi/sigma2))
   print(paste(paste(c(iter, mcemstage, G), collapse=" "), lik) ); likseq[iter] <- lik
    # compute var(beta)
   if ((iter==maxstep) | ((mcemstage==4) & iter==startstage+pls))
    {  varFix = sigma2 * solve(inv.varbeta1 - inv.varbeta2/sigma2)
      vF = sigma2 * solve(inv.varbeta1)  
      }

    # Update MCMC sample size/stage:
   if( (mcemstage==1) & (iter >= 2) & (lik < old.lik) )
    { mcemstage = 2; startstage = iter+1 }
   if( (mcemstage==2) & (iter == startstage + iter2) )
   {
    sdl2 = sqrt(ar(likseq[startstage:iter], order.max=1, aic=FALSE)$var.pred)
    Gmax = round((sdl2/sdl)^2*G) 
    print(paste("max MCMC sample size = ", Gmax))
    Gmax = min(mcmcmax, Gmax)
    deltaG = (1/sqrt(mcmc0)-1/sqrt(Gmax))/trs
    mcemstage = 3; startstage = iter+1
    G = min(Gmax, round(1/(1/sqrt(G)-deltaG)^2))
   }
   if (mcemstage==3)
   {
    if(lflag==0) lflag = sign(lik-old.lik)
    else if ( lflag*sign(lik-old.lik) == -1 )
    {
     lflag = 0
     if (G >= Gmax | G < 0) { mcemstage=4; startstage = iter+1; G=Gmax}
     G = round(1/(1/sqrt(G)-deltaG)^2)
     if (G < 0 | G > Gmax) G = Gmax
    }
   }
   if ( (mcemstage==4) & (iter == startstage+pls) ) 
   {
    print(sqrt(ar(likseq[startstage:iter], order.max=1, aic=FALSE)$var.pred))
    lflag=2; break
   }
   old.lik <- lik
  }  # for 
  if (lflag != 2)
   warning(paste("Convergence not achieved after", maxstep, "steps"))
   #print(paste("log-lik =", round(lik,8), "  Step = ", iter))
    #*****************************************************************
  if (lflag==2) 
  return(list(beta=beta, bi=b, sigma=sqrt(sigma2), Psi=Psi, Delta=Delta,
  loglik = lik, varFix = varFix, method=method, varstruct=varstruct, step=iter, likseq=likseq[1:iter], mcmc=c(mcmc0,G)))
 } # if method == "MLmcmc"
}

