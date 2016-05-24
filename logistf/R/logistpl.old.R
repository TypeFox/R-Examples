logistpl.old <-
function(x, y, init=NULL, i, LL.0, firth, which = -1, offset=rep(0, length(y)), weight=rep(1,length(y)), plcontrol)
{
## which -1...left, +1...right
    k <- ncol(x)
    if (is.null(init)) init<-rep(0,k)
        beta<-init
    if (missing(plcontrol)) plcontrol<-logistpl.control()
    ortho<-plcontrol$ortho
    pr<-plcontrol$pr
    maxit<-plcontrol$maxit
    maxstep<-plcontrol$maxstep
    maxhs<-plcontrol$maxhs
    xconv<-plcontrol$xconv
    lconv<-plcontrol$lconv
    
    
    if(ortho==TRUE & k>1 & i>1){
       thecol<-x[,i]
       others<-x[,-i]
       x[,i]<-lm(thecol~others-1,weights=weight)$residuals
    }
    if(pr==TRUE & (k-1 > 1)){
       others<-x[,c(-1,-i)]
       pc1<-prcomp(others)
       x[,c(-1,-i)]<-predict(pc1,others)
    }
    
    iter <- 0
    betahist<-numeric(0)
    pi <- as.vector(1/(1 + exp( - x %*% beta - offset)))
    XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
    #### X' (W ^ 1/2)
    Fisher <- crossprod(t(XW2)) #### X' W  X
#    loglik <- sum(y * log(pi) + (1 - y) * log(1 - pi))
    loglik <- sum(weight[y==1]*log(pi[y==1]))+sum(weight[y==0]*log(1-pi[y==0]))
   if(firth)
        loglik <- loglik + 0.5 * logDet(Fisher)
    repeat {
        iter <- iter + 1
        covs <- invFisher(Fisher)
        H <- crossprod(XW2, covs) %*% XW2
        if(firth)
            U.star <- crossprod(x, weight*(y - pi) +
                diag(H) * (0.5 - pi))
        else U.star <- crossprod(x, weight*(y - pi))
        V.inv <-  - covs
        lambda <- which * ((2 * ((LL.0 - loglik
            ) + 0.5 * crossprod(U.star,
            V.inv) %*% U.star))/V.inv[i, i]
            )^0.5
        delta <-  - V.inv %*% (U.star + lambda * diag(k)[i,  ])
        delta[is.na(delta)]<-0
        mx <- max(abs(delta))/maxstep
        if(mx > 1)
            delta <- delta/mx
        beta <- beta + delta
        loglik.old <- loglik
        hs<-0
        repeat {
         pi <- as.vector(1/(1 + exp( - x %*% beta - offset)))
#        loglik <- sum(y * log(pi) + (1 - y) *  log(1 - pi))
         loglik <- sum(weight[y==1]*log(pi[y==1]))+sum(weight[y==0]*log(1-pi[y==0]))

         if(firth) {
             XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
    #### X' (W ^ 1/2)
             Fisher <- crossprod(t(XW2))
    #### X' W  X
             loglik <- loglik + 0.5 * logDet(Fisher)
          }
         if((hs>maxhs)|((abs(loglik-LL.0)<abs(loglik.old-LL.0)) & (loglik>LL.0))) break
         beta<-beta - delta/2   ### a simple emergency step halfing
         delta<-delta/2
         hs<-hs+1
         }
        betahist<-rbind(betahist,t(beta))
        if(iter == maxit | ((abs(loglik - LL.0) <=  lconv)  & (max(abs(delta))<xconv)))
            break
    }
    list(beta = beta[i], LL = loglik, conv=c(abs(loglik - LL.0),  max(abs(delta))), iter=iter, betahist=betahist)
}

