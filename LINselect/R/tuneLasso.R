tuneLasso <- function
### tune the  lasso parameter in the 
### regression model : \eqn{Y= X \beta + \sigma N(0,1)}
### using the lasso or the gauss-lasso method
  ##references<<  See Baraud et al. 2010 
  ## \url{http://hal.archives-ouvertes.fr/hal-00502156/fr/} \cr
  ## Giraud et al., 2013,
  ## \url{http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.ss/1356098553}
(Y, ##<<vector with n components : response variable. 
 X, ##<<matrix with n rows and p columns : covariates.
 normalize=TRUE, ##<<logical : corresponds to the input \code{normalize}
 ## of the functions \code{\link[elasticnet]{enet}} and \code{\link[elasticnet]{cv.enet}}. \cr If TRUE the variates \code{X} are
 ## normalized.
 method=c("lasso","Glasso"), ##<<vector of characters whose components are subset of 
  ## (\dQuote{lasso}, \dQuote{Glasso})
 dmax=NULL,##<<integer : maximum number of variables in the lasso
 ##estimator. \code{dmax} \eqn{\le} D where \cr 
 ##  D = min (3*p/4 , n-5) if  p\eqn{ \ge }n
 ## \cr D= min(p,n-5) if
 ## p < n. \cr Default : \code{dmax} = D.
 Vfold=TRUE,##<<logical : if TRUE the tuning is done by Vfold-CV
 V=10,##<<integer. Gives the value of V in the Vfold-CV procedure
 LINselect=TRUE,##<<logical : if TRUE the tuning is done by LINselect
 a=0.5,##<<scalar : value of the parameter \eqn{\alpha} in the LINselect criteria
 K=1.1,##<<scalar : value of the parameter \eqn{K} in the LINselect criteria
 verbose=TRUE,##<<logical : if TRUE a trace of the current process is displayed in real time.
 max.steps=NULL##<<integer : maximum number of steps in the lasso
 ##procedure. \cr Corresponds to the input \code{max.steps} of the function
 ##\code{\link[elasticnet]{enet}}. \cr
 ##Default :
 ##\code{max.steps} = 2*min(p,n)
 )
{
  #
  #
  # packages needed : elasticnet
  calc.proj <- FALSE
  calc.pen <- FALSE
  res <- NULL
  n=nrow(X)
  p=ncol(X)
  Xmean <- apply(X,2,mean)
  mu=mean(Y)
  if (is.null(max.steps)) max.steps=2*min(p,n)
  if (p>=n) dmax <- floor(min(c(3*p/4,n-5,dmax)))
  if (p<n) dmax <- floor(min(c(p,n-5,dmax)))
#   library(elasticnet)
      ##note<<  library \code{elasticnet} is loaded.
  
  res.lasso <- try(enet(X,Y,lambda=0, intercept=TRUE, normalize=normalize,
                        max.steps=max.steps))
  if (inherits(res.lasso, "try-error")) {
    print("error  when calling enet")
    res <- res.lasso[1]
  }
  if (!inherits(res.lasso, "try-error")) {
    lesBeta <- array(0,c(dim(res.lasso$beta.pure)[1],p))
    lesBeta[,res.lasso$allset] <- res.lasso$beta.pure
    lesDim <- apply(lesBeta!=0,1,sum)
    if (max(lesDim) <= dmax) Nmod.lasso <- length(lesDim)
    if (max(lesDim)>dmax) Nmod.lasso<-(1:dim(lesBeta)[1])[(lesDim>=dmax)][1]
    f.lasso <- matrix(0,nrow=n,ncol=Nmod.lasso)
    Intercept <- mu-lesBeta[1:Nmod.lasso,]%*%Xmean
    for (l in 1:Nmod.lasso) {
      f.lasso[,l] <- rep(Intercept[l],n)+X%*%lesBeta[l,]
    }
    if (is.element("lasso",method)) {
      res<- list(lasso=NULL)
      if (Vfold) {
   # lasso + V fold CV
        if (verbose) {
          print(paste("Tuning Lasso with Vfold CV: V=",V,", dmax=",dmax))
        }
        s <- 1:(Nmod.lasso+1)
        res.cv.enet=0
        while (!is.list(res.cv.enet)) {
          s <- s[-length(s)]
          res.cv.enet <- try(cv.enet(X, Y, K = V, lambda=0,
                                     s=s, mode="step", normalize=normalize,
                                     intercept=TRUE,
                                     plot.it=FALSE,se=FALSE),silent=TRUE)
        }
        l.lasso <- which.min(res.cv.enet$cv)
        f.pred <- f.lasso[,l.lasso]
        coeff <- c(Intercept[l.lasso],lesBeta[l.lasso,][lesBeta[l.lasso,]!=0])
        names(coeff) <- c("Intercept",which(lesBeta[l.lasso,]!=0))
        res$lasso$CV <- list(support=which(lesBeta[l.lasso,]!=0),
                             coeff=coeff, fitted=f.pred, crit=res.cv.enet$cv,
                             crit.error=res.cv.enet$cv.error,
                             lambda=res.lasso$penalty[s])
      }
      if (LINselect) {
   # lasso + LinSelect
        calc.proj <- TRUE
        if (verbose) print(paste("Tuning Lasso with LINselect: K=",K,", a=",a,", dmax=",dmax))
 #       pen <- penalite(p, n, K=K, dmax+1)
        D <- 0:dmax
        dm <- pmin(D,rep(p/2,dmax+1))
        Delta <- lgamma(p+1)-lgamma(dm+1)-lgamma(p-dm+1)+2*log(D+1)
        pen <- penalty(Delta, n, p, K)
        calc.pen <- TRUE
        I.lasso <- list(NULL)
        ProjMod <- array(0,c(n,n,Nmod.lasso))
        A <- array(Inf,c(Nmod.lasso,Nmod.lasso))
        B=A
        SCR <- rep(0,Nmod.lasso)
        penSCR <- rep(0,Nmod.lasso)
        sumY2 <- sum((Y-mu)**2)
        un <- rep(1,n)
        for (l in 1:Nmod.lasso) {
          I.lasso[[l]] <- (1:p)[lesBeta[l,]!=0]
          if (length(I.lasso[[l]])>0) {
            ProjM <-  Proj(cbind(un,X[,I.lasso[[l]]]),length(I.lasso[[l]])+1)
            ProjMod[,,l] <-  ProjM$Proj
            SCR[l] <- sum((Y-ProjMod[,,l]%*%Y)**2)
            penSCR[l] <- SCR[l]*pen[ProjM$rg+1]/(n-ProjM$rg)
          }
          if (length(I.lasso[[l]])==0) {
            ProjMod[,,l] <-  un%*%t(un)/n
            ProjM <- list(Proj=ProjMod[,,l],rg=1)
            SCR[l] <- sumY2
            penSCR[l] <- sumY2*pen[ProjM$rg+1]/(n-ProjM$rg)
          }
        }
        Ind <- rep(1:Nmod.lasso,rep(Nmod.lasso,Nmod.lasso))
        YY<-Y%*%t(rep(1,Nmod.lasso))
        for (m in 1:Nmod.lasso) {
          Proj.f<-ProjMod[,,m]%*%f.lasso
          B[m,] <- apply((YY-Proj.f)**2,2,sum)+a*apply((f.lasso-Proj.f)**2,2,sum)
          A[m,] <- B[m,]+penSCR[m]
        }
        l.lasso <- Ind[which.min(A)]
        f.pred <- f.lasso[, l.lasso]
        coeff <- c(Intercept[l.lasso],lesBeta[l.lasso,][lesBeta[l.lasso,]!=0])
        names(coeff) <- c("Intercept",which(lesBeta[l.lasso,]!=0))
        crit <- apply(A,2,min)
        res$lasso$Ls <- list(support=which(lesBeta[l.lasso,]!=0),
                             coeff=coeff, fitted=f.pred, crit=crit,
                             lambda=res.lasso$penalty[1:Nmod.lasso])
      }
    }
    if (is.element("Glasso",method)) {
      if (!is.list(res)) res <- list(Glasso=NULL)
       if (!calc.pen) {
         D <- 0:dmax
         dm <- pmin(D,rep(p/2,dmax+1))
         Delta <- lgamma(p+1)-lgamma(dm+1)-lgamma(p-dm+1)+2*log(D+1)
         pen <- penalty(Delta, n, p, K)
#         pen <- penalite(p, n, K=K, dmax+1)
       }
       if (!calc.proj) {
         I.lasso <- list(NULL)
         SCR <- rep(0,Nmod.lasso)
         penSCR <- rep(0,Nmod.lasso)
         sumY2 <- sum((Y-mu)**2)
         un <- rep(1,n)
         for (l in 1:Nmod.lasso) {
           I.lasso[[l]] <- (1:p)[lesBeta[l,]!=0]
           if (length(I.lasso[[l]])>0) {
             ProjM <-  Proj(cbind(un,X[,I.lasso[[l]]]),length(I.lasso[[l]])+1)
             SCR[l] <- sum((Y-ProjM$Proj%*%Y)**2)
             penSCR[l] <- SCR[l]*pen[ProjM$rg+1]/(n-ProjM$rg)
           }
           if (length(I.lasso[[l]])==0) {
             SCR[l] <- sumY2
             penSCR[l] <- SCR[l]*pen[2]/(n-1)
           }
         }
       }
       if (Vfold) {
  # Gauss lasso + V fold CV
         if (verbose) {
           print(paste("Tuning Gauss Lasso with Vfold CV: V=",V,", dmax=",dmax))
         }
         I.CV <- list(NULL)
         IM.CV <- list(NULL)
         for (iv in 1:V) {
           I.CV[[iv]] <- (floor((iv-1)*n/V)+1):(iv*n/V)
           IM.CV[[iv]] <- (1:n)[-I.CV[[iv]]]
         }
         active_set<-list(NULL)
         all_lambda<-c(-1)
         lambda <- list(NULL)
         for (i in 1:V) {
           in_s<-enleve_var_0(X,IM.CV[[i]])
           res.i<-enet(X[IM.CV[[i]],in_s],Y[IM.CV[[i]]],lambda=0,
                intercept=TRUE,normalize=normalize,max.steps=max.steps)
           lesBeta.i <- array(0,c(dim(res.i$beta.pure)[1],length(in_s)))
           lesBeta.i[,res.i$allset] <- res.i$beta.pure
           lesDim.i <- apply(lesBeta.i!=0,1,sum)
           if (max(lesDim.i) <= dmax) Nmod.lasso.i <- length(lesDim.i)
           if (max(lesDim.i)>dmax) Nmod.lasso.i<-(1:dim(lesBeta.i)[1])[(lesDim.i>=dmax)][1]
           lambda[[i]]<-res.i$penalty[1:Nmod.lasso.i]
           active_set[[i]]<-active(res.i,in_s)[1:Nmod.lasso.i]
           all_lambda<-c(all_lambda,lambda[[i]])
         }
#         all_lambda<-all_lambda[2:(length(all_lambda))]
#         all_lambda<-sort(all_lambda,decreasing=TRUE)
         all_lambda<-all_lambda[2:(length(all_lambda))]
         L1 <- length(all_lambda)
         all_lambda<-c(all_lambda,res.lasso$penalty[1:Nmod.lasso])
         L<-length(all_lambda)
         rank.all_lambda<-L-rank(all_lambda)+1
         sort.all_lambda<-sort(all_lambda,decreasing=TRUE)
         erreur<-matrix(0,V,L)
         for (i in 1:V) {
           est<-estime2(Y,X,I.CV[[i]],IM.CV[[i]],active_set[[i]])
           ind<-1
           for (l in 1:L) {
             if (lambda[[i]][ind]>sort.all_lambda[l]) {
               if (ind < length(lambda[[i]])) {
                 ind<-ind+1
               }	
             }
             erreur[i,l]<-est[[ind]]
           }
         }
         resultat<-apply(erreur,2,mean)
         err.resultat<-sqrt(apply(erreur,2,var))/sqrt(V)
         resultat[resultat<10**(-10)] <- 0
         lambda_cv<-sort.all_lambda[which.min(resultat)]
         ind<-indice(lambda_cv,res.lasso$penalty[1:Nmod.lasso])
         rlm <- lm(Y~X[,I.lasso[[ind]]])
         beta <- rlm$coef
         names(beta) <- c("Intercept",I.lasso[[ind]])
         resultatF <- c(min(resultat),
                        resultat[rank.all_lambda][(L1+1):L])
         err.resultatF <- c(err.resultat[which.min(resultat)],
                            err.resultat[rank.all_lambda][(L1+1):L])
         les.lambda <- c(lambda_cv,res.lasso$penalty[1:Nmod.lasso])
         crit <- resultatF[order(les.lambda,decreasing=TRUE)]
         crit.error <- err.resultatF[order(les.lambda,decreasing=TRUE)]
#         ii <- length(resultat)
#         if (min(resultat)==0) ii <- which(resultat==0)[1]
         res$Glasso$CV <- list(support=I.lasso[[ind]],
                               coeff=beta, fitted=rlm$fitted,
                               crit=crit,
                               crit.error=crit.error,
                               lambda=sort(les.lambda,decreasing=TRUE))
       }
       if (LINselect) {
  # Gauss lasso + LinSelect
         if (verbose) print(paste("Tuning Gauss Lasso with LinSelect: K=",K))
         if (!calc.pen) {
           D <- 0:dmax
           dm <- pmin(D,rep(p/2,dmax+1))
           Delta <- lgamma(p+1)-lgamma(dm+1)-lgamma(p-dm+1)+2*log(D+1)
           pen <- penalty(Delta, n, p, K)
#           pen <- penalite(p, n, K=K, dmax)
         }
         crit <- SCR + penSCR
         ind <- which.min(crit)
         rlm <- lm(Y~X[,I.lasso[[ind]]])
         beta <- rlm$coef
         names(beta) <- c("Intercept",I.lasso[[ind]])
         res$Glasso$Ls <- list(support=I.lasso[[ind]],
                               coeff=beta, fitted=rlm$fitted+mu,
                               crit=crit,
                               lambda=res.lasso$penalty[1:Nmod.lasso] )
       }
     }
  }
  result <- res
  return(result)
  ### A list with one or two components according to 
  ### \code{method}. \cr
  ### \code{lasso} if \code{method} contains "lasso" is a list with  one or two components
  ### according to \code{Vfold} and \code{LINselect}.
  ###
  ### \itemize{
  ###   \item{\code{Ls} if \code{LINselect}=TRUE. A list with components
  ###     \itemize{
  ###       \item{\code{support}: vector of integers. Estimated support of the
  ###              parameter vector \eqn{\beta}.}
  ###       \item{\code{coef}: vector whose first component is the estimated
  ###             intercept. \cr The other components are the estimated non zero
  ###             coefficients.}
  ###       \item{\code{fitted}: vector with length n. Fitted value of
  ###       the response.}
  ###       \item{\code{crit}: vector containing the values of the criteria for
  ###       each value of \code{lambda}.}
  ###       \item{\code{lambda}: vector containing the values of the tuning
  ###       parameter of the lasso algorithm.}
  ###     }
  ###   }
  ###   \item{\code{CV} if \code{Vfold}=TRUE. A list with components
  ###     \itemize{
  ###       \item{\code{support}: vector of integers. Estimated support of the
  ###              parameter vector \eqn{\beta}.}
  ###       \item{\code{coef}: vector whose first component is the estimated
  ###             intercept. \cr The other components are the estimated non zero
  ###             coefficients.}
  ###       \item{\code{fitted}: vector with length n. Fitted value of
  ###       the response.}
  ###       \item{\code{crit}: vector containing the values of the criteria for
  ###       each value of \code{lambda}.}
  ###       \item{\code{crit.err}: vector containing the estimated
  ###       standard-error of the criteria.}
  ###       \item{\code{lambda}: vector containing the values of the tuning
  ###       parameter of the lasso algorithm.}
  ###     }
  ###   }
  ### }
  ### \code{Glasso} if \code{method} contains "Glasso". The same as \code{lasso}.
}

