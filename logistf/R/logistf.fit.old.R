logistf.fit.old <-#-------------------------------------------------------------------------------
logDet    <- function(x) 2*sum(log(diag(chol(x))));
#-------------------------------------------------------------------------------
invFisher <- function(x) chol2inv(chol(x));
#-------------------------------------------------------------------------------



logistf.fit.old<-function(x, y, weight=NULL, offset=NULL, firth=TRUE, col.fit=NULL, init=NULL, control, collapse=TRUE){
# fitter function for logistf

# lconv = convergence criterion for log likelihood
# gconv = convergence criterion for score
# xconv = convergence criterion for parameter estimates
# pos = columns in x which will not be estimated (left at init value)


n<-nrow(x)
k<-ncol(x)

coll<-FALSE
if(collapse & length(unique(weight))==1 & weight[1]==1) {
  require(mgcv)
  xc<-uniquecombs(cbind(x,y,offset))
  xorig<-x
  yorig<-y
  weight<-table(attr(xc,"index"))
  x<-xc[,1:k]
  y<-xc[,k+1]
  if(!is.null(offset)) offset<-xc[,k+2]
  n<-nrow(xc)
  coll<-TRUE
  }


if (is.null(init)) init=rep(0,k)
if (is.null(col.fit)) col.fit=1:k
if (is.null(offset)) offset=rep(0,n)
if (is.null(weight)) weight=rep(1,n)
if (col.fit[1]==0) maxit<-0   #only evaluate likelihood and go back
if (missing(control)) control<-logistf.control()

maxit<-control$maxit
maxstep<-control$maxstep
maxhs<-control$maxhs
lconv<-control$lconv
gconv<-control$gconv
xconv<-control$xconv

#pos <- col.fit
beta <- init
l.change <- 5

    iter <- 0
    pi <- as.vector(1/(1 + exp( - x %*% beta - offset)))
#    loglik <- sum(y * log(pi) + (1 - y) * log(1 - pi))
    loglik <- sum(weight[y==1]*log(pi[y==1]))+sum(weight[y==0]*log(1-pi[y==0]))

    if(firth) {
        XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)    #### X' (W ^ 1/2)
        Fisher <- crossprod(t(XW2)) #### X' W  X
        loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
    }
    evals<-1
    repeat {
#        if(k2 == k) break   ## -> Overall Test
        loglik.old <- loglik
        beta.old <- beta
        XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)    #### X' (W ^ 1/2)
        Fisher <- crossprod(t(XW2)) #### X' W  X
        covs <- invFisher(Fisher)   ### (X' W  X) ^ -1
        Hdiag<-sapply(1:n,function(ii) crossprod((XW2[,ii]),crossprod(covs,(XW2[,ii]))))
#        for(i in 1:n) Hdiag[i]<-crossprod((XW2[,i]),crossprod(covs,(XW2[,i])))
#        H <- crossprod(XW2, covs)  %*% XW2
#        if(firth) U.star <- crossprod(x, weight*(y - pi) + diag(H) * (0.5 - pi))
        if(firth) U.star <- crossprod(x, weight*(y - pi) + Hdiag * (0.5 - pi))
        else U.star <- crossprod(x, weight*(y - pi))
        XX.covs <- matrix(0, k, k)
        if (col.fit[1] != 0){
               XX.XW2 <- crossprod(x[, col.fit, drop = FALSE], diag(weight * pi * (1 - pi))^0.5)
             #### Teil von X' (W ^ 1/2)
               XX.Fisher <- crossprod(t(XX.XW2))   #### Teil von  X' W  X
               XX.covs[col.fit, col.fit] <- invFisher(XX.Fisher)
              ### aufblasen der Cov-Matrix
        }
        delta <- as.vector(XX.covs %*% U.star)
        delta[is.na(delta)]<-0
        mx <- max(abs(delta))/maxstep
        if(mx > 1)
            delta <- delta/mx
         evals<-evals+1


        if(maxit > 0){
        iter <- iter + 1
        beta <- beta + delta
        for(halfs in 1:maxhs) {
## Half-Steps
            pi <- as.vector(1/(1 + exp( - x %*% beta - offset)))
            loglik <- sum(weight[y==1]*log(pi[y==1]))+sum(weight[y==0]*log(1-pi[y==0]))
            if(firth) {
                XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
                #### X' (W ^ 1/2)
                Fisher <- crossprod(t(XW2)) #### X' W  X
#                loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
                loglik <- loglik + 0.5 * logDet(Fisher)

            }
            evals<-evals+1
            l.change<-loglik-loglik.old

            if(loglik > loglik.old)
                break
            beta <- beta - delta * 2^( - halfs)
            ##beta-Aenderung verkleinern
         }
        }
        if(iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star)<gconv)) & (all(l.change<lconv))))
            break
    }

    var <- XX.covs



# list(beta=beta, var=var, pi=pi, hat.diag=diag(H), loglik=loglik, iter=iter, evals=evals, conv=c(l.change, max(abs(U.star)), max(abs(delta))))

 if(coll) {
  pi<-pi[attr(xc,"index")]
  Hdiag<-Hdiag[attr(xc,"index")]
  }

 list(beta=beta, var=var, pi=pi, hat.diag=Hdiag, loglik=loglik, iter=iter, evals=evals, conv=c(l.change, max(abs(U.star)), max(abs(delta))))
}

