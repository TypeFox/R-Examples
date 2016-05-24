############################################################
## summary functions for scam() (clone of summary.gam())...
## of mgcv version 1.7-22, 
## but without p.type=0,1,...,5 as summary input variable,
## only with freq=T/F....
###########################################################

##### mgcv::: smoothTest

smoothTest <- function(b,X,V,eps=.Machine$double.eps^.5) {
## Forms Cox, Koh, etc type test statistic, and
## obtains null distribution by simulation...
## if b are coefs f=Xb, cov(b) = V. z is a vector of 
## i.i.d. N(0,1) deviates

  qrx <- qr(X)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  k <- n <- length(ed$values)
  ## could truncate, but it doesn't improve power in correlated case!
  f <- t(ed$vectors[,1:k])%*%R%*%b
  t <- sum(f^2)
  k <- ncol(X)
  lambda <- as.numeric(ed$values[1:k])
  pval <- liu2(t,lambda) ## should really use Davies
  list(stat=t,pval=pval)  
} 


###### mgcv::: liu2


liu2 <- function(x, lambda, h = rep(1,length(lambda)),lower.tail=FALSE) {
# Evaluate Pr[sum_i \lambda_i \chi^2_h_i < x] approximately.
# Code adapted from CompQuadForm package of Pierre Lafaye de Micheaux 
# and directly from....
# H. Liu, Y. Tang, H.H. Zhang, A new chi-square approximation to the 
# distribution of non-negative definite quadratic forms in non-central 
# normal variables, Computational Statistics and Data Analysis, Volume 53, 
# (2009), 853-856. Actually, this is just Pearson (1959) given that
# the chi^2 variables are central. 
# Note that this can be rubbish in lower tail (e.g. lambda=c(1.2,.3), x = .15)
  
#  if (TRUE) { ## use Davies exact method in place of Liu et al/ Pearson approx.
#    require(CompQuadForm)
#    r <- x
#    for (i in 1:length(x)) r[i] <- davies(x[i],lambda,h)$Qq
#    return(pmin(r,1))
#  }

  if (length(h) != length(lambda)) stop("lambda and h should have the same length!")
 
  lh <- lambda*h
  muQ <- sum(lh)
  
  lh <- lh*lambda
  c2 <- sum(lh)
  
  lh <- lh*lambda
  c3 <- sum(lh)
  
  s1 <- c3/c2^1.5
  s2 <- sum(lh*lambda)/c2^2

  sigQ <- sqrt(2*c2)

  t <- (x-muQ)/sigQ

  if (s1^2>s2) {
    a <- 1/(s1-sqrt(s1^2-s2))
    delta <- s1*a^3-a^2
    l <- a^2-2*delta
  } else {
    a <- 1/s1
    delta <- 0
    l <- c2^3/c3^2
  }

  muX <- l+delta
  sigX <- sqrt(2)*a
  
  return(pchisq(t*sigX+muX,df=l,ncp=delta,lower.tail=lower.tail))

}


#### mgcv::: simf

simf <- function(x,a,df,nq=50) {
## suppose T = sum(a_i \chi^2_1)/(chi^2_df/df). We need
## Pr[T>x] = Pr(sum(a_i \chi^2_1) > x *chi^2_df/df). Quadrature 
## used here. So, e.g.
## 1-pf(4/3,3,40);simf(4,rep(1,3),40);1-pchisq(4,3)
  p <- (1:nq-.5)/nq
  q <- qchisq(p,df)
  x <- x*q/df
  pr <- sum(liu2(x,a)) ## Pearson/Liu approx to chi^2 mixture
  pr/nq 
}

#### the same as mgcv::: recov.gam

recov.scam <- function(b,re=rep(0,0),m=0) {
## b is a fitted gam object. re is an array of indices of 
## smooth terms to be treated as fully random....
## Returns frequentist Cov matrix based on the given
## mapping from data to params, but with dist of data
## corresponding to that implied by treating terms indexed
## by re as random effects... (would be usual frequentist 
## if nothing treated as random)
## if m>0, then this is indexes a term, not in re, whose
## unpenalized cov matrix is required, with the elements of re
## dropped.
  if (!inherits(b,"scam")) stop("recov works with fitted scam objects only") 
  if (is.null(b$full.sp)) sp <- b$sp else sp <- b$full.sp
  if (length(re)<1) { 
    if (m>0) {
      ## annoyingly, need total penalty  
      np <- length(coef(b))
      k <- 1;S1 <- matrix(0,np,np)
      for (i in 1:length(b$smooth)) { 
        ns <- length(b$smooth[[i]]$S)
        ind <- b$smooth[[i]]$first.para:b$smooth[[i]]$last.para
        if (ns>0) for (j in 1:ns) {
          S1[ind,ind] <- S1[ind,ind] + sp[k]*b$smooth[[i]]$S[[j]]
          k <- k + 1
        }
      }
      LRB <- rbind(b$R,t(mroot(S1)))
      ii <- b$smooth[[m]]$first.para:b$smooth[[m]]$last.para 
      ## ii is cols of LRB related to smooth m, which need 
      ## to be moved to the end...
      LRB <- cbind(LRB[,-ii],LRB[,ii])
      ii <- (ncol(LRB)-length(ii)+1):ncol(LRB)
      Rm <- qr.R(qr(LRB,tol=0,LAPACK=FALSE))[ii,ii] ## unpivoted QR
    } else Rm <- NULL
    return(list(Ve.t=(t(b$Ve.t)+b$Ve.t)*.5,Rm=Rm)) 
  }

  if (m%in%re) stop("m can't be in re")
  ## partition R into R1 ("fixed") and R2 ("random"), with S1 and S2
  p <- length(b$coefficients)
  rind <- rep(FALSE,p) ## random coefficient index
  for (i in 1:length(re)) {
    rind[b$smooth[[re[i]]]$first.para:b$smooth[[re[i]]]$last.para] <- TRUE
  }
  p2 <- sum(rind) ## number random
  p1 <- p - p2 ## number fixed
  map <- rep(0,p) ## remaps param indices to indices in split version
  map[rind] <- 1:p2 ## random
  map[!rind] <- 1:p1 ## fixed
  
  ## split R...
  R1 <- b$R[,!rind]  ## fixed effect columns
  R2 <- b$R[,rind]   ## random effect columns
  ## seitdem ich dich kennen, hab ich ein probleme,

  ## assemble S1 and S2
  S1 <- matrix(0,p1,p1);S2 <- matrix(0,p2,p2)
 
  k <- 1
  for (i in 1:length(b$smooth)) { 
    ns <- length(b$smooth[[i]]$S)
    ind <- map[b$smooth[[i]]$first.para:b$smooth[[i]]$last.para]
    is.random <- i%in%re
    if (ns>0) for (j in 1:ns) {
      if (is.random) S2[ind,ind] <- S2[ind,ind] +  sp[k]*b$smooth[[i]]$S[[j]] else
         S1[ind,ind] <- S1[ind,ind] + sp[k]*b$smooth[[i]]$S[[j]]
      k <- k + 1
    }
  }
  ## pseudoinvert S2
  if (nrow(S2)==1) {
    S2[1,1] <- 1/sqrt(S2[1,1])
  } else if (max(abs(diag(diag(S2))-S2))==0) {
    ds2 <- diag(S2)
    ind <- ds2 > max(ds2)*.Machine$double.eps^.8
    ds2[ind] <- 1/ds2[ind];ds2[!ind] <- 0
    diag(S2) <- sqrt(ds2)
  } else {
    ev <- eigen((S2+t(S2))/2,symmetric=TRUE)
    ind <- ev$values > max(ev$values)*.Machine$double.eps^.8
    ev$values[ind] <- 1/ev$values[ind];ev$values[!ind] <- 0 
    ## S2 <- ev$vectors%*%(ev$values*t(ev$vectors))
    S2 <- sqrt(ev$values)*t(ev$vectors)
  }
  ## choleski of cov matrix....
  ## L <- chol(diag(p)+R2%*%S2%*%t(R2)) ## L'L = I + R2 S2^- R2'
  L <- chol(diag(p) + crossprod(S2%*%t(R2)))  

  ## now we need the square root of the unpenalized
  ## cov matrix for m
  if (m>0) {
      ## llr version
      LRB <- rbind(L%*%R1,t(mroot(S1)))
      ii <- map[b$smooth[[m]]$first.para:b$smooth[[m]]$last.para] 
      ## ii is cols of LRB related to smooth m, which need 
      ## to be moved to the end...
      LRB <- cbind(LRB[,-ii],LRB[,ii])
      ii <- (ncol(LRB)-length(ii)+1):ncol(LRB) ## need to pick up final block
      Rm <- qr.R(qr(LRB,tol=0,LAPACK=FALSE))[ii,ii,drop=FALSE] ## unpivoted QR
  } else Rm <- NULL

  list(Ve.t= crossprod(L%*%b$R%*%b$Vp.t)/b$sig2, ## Frequentist cov matrix
       Rm=Rm)
 # mapi <- (1:p)[!rind] ## indexes mapi[j] is index of total coef vector to which jth row/col of Vb/e relates
  
} ## end of recov


### same as mgcv::: reTest.scam

reTest.scam <- function(b,m) {
## Test the mth smooth for equality to zero
## and accounting for all random effects in model 
  
  ## find indices of random effects other than m
  rind <- rep(0,0)
  for (i in 1:length(b$smooth)) if (!is.null(b$smooth[[i]]$random)&&b$smooth[[i]]$random&&i!=m) rind <- c(rind,i)
  ## get frequentist cov matrix of effects treating smooth terms in rind as random
  rc <- recov.scam(b,rind,m) 
  Ve.t <- rc$Ve.t
  ind <- b$smooth[[m]]$first.para:b$smooth[[m]]$last.para
  B <- mroot(Ve.t[ind,ind,drop=FALSE]) ## BB'=Ve
 
  Rm <- rc$Rm
  
  b.hat <- coef(b)[ind]
  d <- Rm%*%b.hat
  stat <- sum(d^2)/b$sig2
  ev <- eigen(crossprod(Rm%*%B)/b$sig2,symmetric=TRUE,only.values=TRUE)$values
  ev[ev<0] <- 0
  rank <- sum(ev>max(ev)*.Machine$double.eps^.8)
  
  if (b$scale.estimated) {
    pval <- simf(stat,ev,b$df.residual)
  } else { pval <- liu2(stat,ev) }
  list(stat=stat,pval=pval,rank=rank)
} ## end reTest


########## mgcv::: testStat

testStat <- function(p,X,V,rank=NULL,type=0,res.df= -1) {
## Routine for forming fractionally trunctated
## pseudoinverse of XVX'. And returning 
## p'X'(XVX)^-Xp.
## Truncates to numerical rank, if this is
## less than supplied rank+1.
## The type argument specifies the type of truncation to use.
## on entry `rank' should be an edf estimate
## 0. Default using the fractionally truncated pinv.
## 1. Round down to k if k<= rank < k+0.05, otherwise up.
## 2. Naive rounding.
## 3. Round up.
## 4. Numerical rank estimation, tol=1e-3
## res.df is residual dof used to estimate scale. <=0 implies
## fixed scale.

  qrx <- qr(X,tol=0)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot,drop=FALSE]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)


  k <- max(0,floor(rank)) 
  nu <- abs(rank - k)     ## fractional part of supplied edf
  if (type < -.5) { ## Crude modification of Cox and Koh
    res <- smoothTest(p,X,V)
    res$rank <- rank
    return(res)
  } else  if (type==1) { ## round up is more than .05 above lower
    if (rank > k + .05||k==0) k <- k + 1
    nu <- 0;rank <- k
  } else if (type==2) { ## naive round
    nu <- 0;rank <- k <- max(1,round(rank))
    warning("p-values may give low power in some circumstances")
  } else if (type==3) { ## round up
    nu <- 0; rank <- k <- max(1,ceiling(rank))
    warning("p-values un-reliable")
  } else if (type==4) { ## rank estimation
    rank <- k <- max(sum(ed$values>1e-3*max(ed$values)),1) 
    nu <- 0
    warning("p-values may give very low power")
  }

  if (nu>0) k1 <- k+1 else k1 <- k

  ## check that actual rank is not below supplied rank+1
  r.est <- sum(ed$values > max(ed$values)*.Machine$double.eps^.9)
  if (r.est<k1) {k1 <- k <- r.est;nu <- 0;rank <- r.est}

  ## Get the eigenvectors...
  # vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
  vec <- ed$vectors
  if (k1<ncol(vec)) vec <- vec[,1:k1,drop=FALSE]

  ## deal with the fractional part of the pinv...
  if (nu>0&&k>0) {
     if (k>1) vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
     b12 <- .5*nu*(1-nu)
     if (b12<0) b12 <- 0
     b12 <- sqrt(b12)
     B <- matrix(c(1,b12,b12,nu),2,2)
     ev <- diag(ed$values[k:k1]^-.5,nrow=k1-k+1)
     B <- ev%*%B%*%ev
     eb <- eigen(B,symmetric=TRUE)
     rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
     vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
  } else {
    if (k==0) vec <- t(t(vec)*sqrt(1/ed$val[1])) else
    vec <- t(t(vec)/sqrt(ed$val[1:k]))
    if (k==1) rank <- 1
  }
 
  d <- t(vec)%*%(R%*%p)
  d <- sum(d^2) 

  rank1 <- rank ## rank for lower tail pval computation below

  ## note that for <1 edf then d is not weighted by EDF, and instead is 
  ## simply refered to a chi-squared 1

  if (nu>0) { ## mixture of chi^2 ref dist
     if (k1==1) rank1 <- val <- 1 else { 
       val <- rep(1,k1) ##ed$val[1:k1]
       rp <- nu+1
       val[k] <- (rp + sqrt(rp*(2-rp)))/2
       val[k1] <- (rp - val[k])
     }
   
     if (res.df <= 0) pval <- liu2(d,val) else ##  pval <- davies(d,val)$Qq else
     pval <- simf(d,val,res.df)
  } else { pval <- 2 }
  ## integer case still needs computing, also liu/pearson approx only good in 
  ## upper tail. In lower tail, 2 moment approximation is better (Can check this 
  ## by simply plotting the whole interesting range as a contour plot!)
  if (pval > .5) { 
    if (res.df <= 0) pval <- pchisq(d,df=rank1,lower.tail=FALSE) else
    pval <- pf(d/rank1,rank1,res.df,lower.tail=FALSE)
  }
  list(stat=d,pval=min(1,pval),rank=rank)
} ## end of testStat




####################################################
##### function to get all the summary information....
#############################################

model.matrix.scam <- function(object,...)
{ if (!inherits(object,"scam")) stop("`object' is not of class \"scam\"")
  predict(object,type="lpmatrix",...)
}



summary.scam <- function (object,dispersion = NULL,freq = FALSE,...) 
{
    pinv <- function(V, M, rank.tol = 1e-06) {
      ## a local pseudoinverse function
        D <- eigen(V,symmetric=TRUE)
        M1<-length(D$values[D$values>rank.tol*D$values[1]])
        if (M>M1) M<-M1 # avoid problems with zero eigen-values
  
        if (M+1<=length(D$values)) D$values[(M+1):length(D$values)]<-1
        D$values<- 1/D$values
        if (M+1<=length(D$values)) D$values[(M+1):length(D$values)]<-0
        res <- D$vectors%*%(D$values*t(D$vectors))  ##D$u%*%diag(D$d)%*%D$v
        attr(res,"rank") <- M
        res
    } ## end of pinv
    
     p.table <- pTerms.table <- s.table <- NULL
    if (freq) covmat <- object$Ve.t  else covmat <- object$Vp.t
    name <- names(object$coefficients.t)
   # name <- names(object$edf)
    dimnames(covmat) <- list(name, name)
    covmat.unscaled <- covmat/object$sig2
    est.disp <- object$scale.estimated
    
    if (!is.null(dispersion)) {
        covmat <- dispersion * covmat.unscaled
        object$Ve.t <- object$Ve.t*dispersion/object$sig2 ## freq
        object$Vp.t <- object$Vp.t*dispersion/object$sig2 ## Bayes
        est.disp <- FALSE
    }
    else dispersion <- object$sig2


    ## Now the individual parameteric coefficient p-values...

    se <- diag(covmat)^0.5
    residual.df <- length(object$y) - sum(object$edf)
    if (object$nsdf > 0) {
        p.coeff <- object$coefficients.t[1:object$nsdf]
        p.se <- se[1:object$nsdf]
        p.t <- p.coeff/p.se
        if (!est.disp) {
            p.pv <- 2 * pnorm(abs(p.t), lower.tail = FALSE)
            p.table <- cbind(p.coeff, p.se, p.t, p.pv)
            dimnames(p.table) <- list(names(p.coeff), c("Estimate", 
                "Std. Error", "z value", "Pr(>|z|)"))
        }
        else {
            p.pv <- 2 * pt(abs(p.t), df = residual.df, lower.tail = FALSE)
            p.table <- cbind(p.coeff, p.se, p.t, p.pv)
            dimnames(p.table) <- list(names(p.coeff), c("Estimate", 
                "Std. Error", "t value", "Pr(>|t|)"))
        }
    }
    else {
        p.coeff <- p.t <- p.pv <- array(0, 0)
    }

    ## Next the p-values for parametric terms, so that factors are treated whole... 
 
    term.labels <- attr(object$pterms, "term.labels")
    nt <- length(term.labels)
    if (nt > 0) {
        np <- length(object$assign)
        Vb <- matrix(covmat[1:np, 1:np],np,np) ## ,drop=FALSE)
        bp <- array(object$coefficients.t[1:np], np)
        pTerms.pv <- array(0, nt)
        attr(pTerms.pv, "names") <- term.labels
        pTerms.df <- pTerms.chi.sq <- pTerms.pv
        for (i in 1:nt) {
            ind <- object$assign == i
            b <- bp[ind]; V <- Vb[ind, ind]
            ## pseudo-inverse needed in case of truncation of parametric space
            if (length(b) == 1) {
                V <- 1/V
                pTerms.df[i] <- nb <- 1
                pTerms.chi.sq[i] <- V * b * b
            }
            else {
                V <- pinv(V, length(b), rank.tol = .Machine$double.eps^0.5)
                pTerms.df[i] <- nb <- attr(V, "rank")
                pTerms.chi.sq[i] <- t(b) %*% V %*% b
            }
            if (!est.disp) 
                pTerms.pv[i] <- pchisq(pTerms.chi.sq[i], df = nb, 
                  lower.tail = FALSE)
            else pTerms.pv[i] <- pf(pTerms.chi.sq[i]/nb, df1 = nb, 
                df2 = residual.df, lower.tail = FALSE)
        }
        if (!est.disp) {
            pTerms.table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)
            dimnames(pTerms.table) <- list(term.labels, c("df", 
                "Chi.sq", "p-value"))
        }
        else {
            pTerms.table <- cbind(pTerms.df, pTerms.chi.sq/pTerms.df, 
                pTerms.pv)
            dimnames(pTerms.table) <- list(term.labels, c("df", 
                "F", "p-value"))
        }
    }
    else {
        pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0, 0)
    }

    ## Now deal with the smooth terms....

    m <- length(object$smooth)  # number of smooth terms
    df <- edf1 <-edf <- s.pv <- chi.sq <- array(0, m)
    if (m > 0) { # form test statistics for each smooth
      if (!freq) { ## Bayesian p-values required 
        sub.samp <- max(1000,2*length(object$coefficients)) 
        if (nrow(object$model)>sub.samp) { ## subsample to get X for p-values calc.
          seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
          if (inherits(seed,"try-error")) {
            runif(1)
            seed <- get(".Random.seed",envir=.GlobalEnv)
          }
          kind <- RNGkind(NULL)
          RNGkind("default","default")
          set.seed(11) ## ensure repeatability
          ind <- sample(1:nrow(object$model),sub.samp,replace=FALSE)  ## sample these rows from X
          X <- predict(object,object$model[ind,],type="lpmatrix")
          RNGkind(kind[1],kind[2])
          assign(".Random.seed",seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
        } else { ## don't need to subsample 
          X <- model.matrix(object)
          }
        X <- X[!is.na(rowSums(X)),] ## exclude NA's (possible under na.exclude)
    
     } ## end if (!freq)

     for (i in 1:m) { ## loop through smooths
        start <- object$smooth[[i]]$first.para
        stop <- object$smooth[[i]]$last.para
            
        if (freq) { ## use frequentist cov matrix 
          V <- object$Ve.t[start:stop,start:stop,drop=FALSE] 
        } else V <- object$Vp.t[start:stop,start:stop,drop=FALSE] ## Bayesian
      
        p <- object$coefficients.t[start:stop] # transposed parameters of a smooth
            
        edf1[i] <- edf[i] <- sum(object$edf[start:stop]) # edf for this smooth
        ## extract alternative edf estimate for this smooth, if possible...
        ## edf1 is not done for scam output value...
        if (!is.null(object$edf1)) edf1[i] <-  sum(object$edf1[start:stop])   
       
        if (freq) {
           M1 <- object$smooth[[i]]$df
           M <- min(M1, ceiling(2 * sum(object$edf[start:stop]))) ## upper limit of 2*edf on rank
           V <- pinv(V, M)  # get rank M pseudoinverse of V
           chi.sq[i] <- t(p) %*% V %*% p
           df[i] <- attr(V, "rank")
        }  else { ## Better founded alternatives...
           Xt <- X[, start:stop,drop=FALSE]
           if (object$smooth[[i]]$null.space.dim==0&&!is.null(object$R)) { ## random effect or fully penalized term
             res <- reTest.scam(object,i)
           } else { ## Inverted Nychka interval statistics
             df[i] <- min(ncol(Xt),edf1[i])
             if (est.disp) rdf <- residual.df else rdf <- -1
             res <- testStat(p,Xt,V,df[i],type=0,res.df = rdf) ## was type=p.type
           }
           df[i] <- res$rank
           chi.sq[i] <- res$stat
           s.pv[i] <- res$pval 
        }
        names(chi.sq)[i]<- object$smooth[[i]]$label
      
        if (freq) {
          if (!est.disp)
            s.pv[i] <- pchisq(chi.sq[i], df = df[i], lower.tail = FALSE)
            else
              s.pv[i] <- pf(chi.sq[i]/df[i], df1 = df[i], df2 = residual.df, lower.tail = FALSE)
              ## p-values are meaningless for very small edf. Need to set to NA
              if (df[i] < 0.1) s.pv[i] <- NA
        }
      }
      if (!est.disp) {
        if (freq) {
          s.table <- cbind(edf, df, chi.sq, s.pv)      
          dimnames(s.table) <- list(names(chi.sq), c("edf", "Est.rank", "Chi.sq", "p-value"))
        } else {
          s.table <- cbind(edf, df, chi.sq, s.pv)      
          dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "Chi.sq", "p-value"))
        }
     } else {
       if (freq) {
         s.table <- cbind(edf, df, chi.sq/df, s.pv)      
         dimnames(s.table) <- list(names(chi.sq), c("edf", "Est.rank", "F", "p-value"))
       } else {
         s.table <- cbind(edf, df, chi.sq/df, s.pv)      
         dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "F", "p-value"))
       }
     }
   }
   w <- as.numeric(object$prior.weights)
   mean.y <- sum(w*object$y)/sum(w)
   w <- sqrt(w)
   nobs <- nrow(object$model)
   r.sq<- 1 - var(w*(as.numeric(object$y)-object$fitted.values))*(nobs-1)/(var(w*(as.numeric(object$y)-mean.y))*residual.df) 
   dev.expl<-(object$null.deviance-object$deviance)/object$null.deviance
   ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,
       s.pv=s.pv,scale=dispersion,r.sq=r.sq,family=object$family,formula=object$formula,n=nobs,
       dev.expl=dev.expl,edf=edf,dispersion=dispersion,pTerms.pv=pTerms.pv,pTerms.chi.sq=pTerms.chi.sq,
       pTerms.df = pTerms.df, cov.unscaled = covmat.unscaled, cov.scaled = covmat, p.table = p.table,
       pTerms.table = pTerms.table, s.table = s.table,method=object$method,sp.criterion=object$gcv.ubre,
       sp=object$sp,dgcv.ubre=object$dgcv.ubre,termcode=object$termcode,
       gcv.ubre=object$gcv.ubre,optimizer=object$optimizer)

    class(ret) <- "summary.scam"
    ret
}  ## end summary.scam

   
################## mgcv::: pinvXVX

pinvXVX <- function (X, V, rank = NULL) 
{
    k <- floor(rank)
    nu <- rank - k
    if (nu > 0) 
        k1 <- k + 1
    else k1 <- k
    qrx <- qr(X)
    R <- qr.R(qrx)
    V <- R %*% V[qrx$pivot, qrx$pivot] %*% t(R)
    V <- (V + t(V))/2
    ed <- eigen(V, symmetric = TRUE)
    vec <- qr.qy(qrx, rbind(ed$vectors, matrix(0, nrow(X) - ncol(X), 
        ncol(X))))
    if (k1 < ncol(vec)) 
        vec <- vec[, 1:k1, drop = FALSE]
    if (k == 0) {
        vec <- t(t(vec) * sqrt(nu/ed$val[1]))
        return(vec)
    }
    if (nu > 0) {
        if (k > 1) 
            vec[, 1:(k - 1)] <- t(t(vec[, 1:(k - 1)])/sqrt(ed$val[1:(k - 
                1)]))
        b12 <- 0.5 * nu * (1 - nu)
        if (b12 < 0) 
            b12 <- 0
        b12 <- sqrt(b12)
        B <- matrix(c(1, b12, b12, nu), 2, 2)
        ev <- diag(ed$values[k:k1]^-0.5)
        B <- ev %*% B %*% ev
        eb <- eigen(B, symmetric = TRUE)
        rB <- eb$vectors %*% diag(sqrt(eb$values)) %*% t(eb$vectors)
        vec[, k:k1] <- t(rB %*% t(vec[, k:k1]))
    }
    else {
        vec <- t(t(vec)/sqrt(ed$val[1:k]))
    }
    vec
}
 
################ mgcv::: eigXVX

eigXVX <- function (X, V, rank = NULL, tol = .Machine$double.eps^0.5) 
{
    qrx <- qr(X)
    R <- qr.R(qrx)
    V <- R %*% V[qrx$pivot, qrx$pivot] %*% t(R)
    V <- (V + t(V))/2
    ed <- eigen(V, symmetric = TRUE)
    ind <- abs(ed$values) > max(abs(ed$values)) * tol
    erank <- sum(ind)
    if (is.null(rank)) {
        rank <- erank
    }
    else {
        if (rank < erank) 
            ind <- 1:rank
        else rank <- erank
    }
    vec <- qr.qy(qrx, rbind(ed$vectors, matrix(0, nrow(X) - ncol(X), 
        ncol(X))))
    list(values = ed$values[ind], vectors = vec[, ind], rank = rank)
}


##### print.summary.scam .....
print.summary.scam <- function (x, digits = max(3, getOption("digits") - 3), 
         signif.stars = getOption("show.signif.stars"),...) 
## print method for scam, a clone of print.summary.gam of mgcv()
{
    print(x$family)
    cat("Formula:\n")
    print(x$formula)
    if (length(x$p.coeff) > 0) {
        cat("\nParametric coefficients:\n")
        printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    cat("\n")
    if (x$m > 0) {
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, 
            has.Pvalue = TRUE, na.print = "NA", cs.ind = 1, ...)
    }
    cat("\nR-sq.(adj) = ", formatC(x$r.sq, digits = 4, width = 5))
    if (length(x$dev.expl) > 0) 
        cat("   Deviance explained = ", formatC(x$dev.expl * 
            100, digits = 3, width = 4), "%\n", sep = "")
    cat( x$method," score = ", formatC(x$sp.criterion, digits = 5), 
            sep = "")
    cat("  Scale est. = ", formatC(x$scale, digits = 5, width = 8, 
        flag = "-"), "  n = ", x$n, "\n", sep = "")
    if (x$optimizer == "bfgs"){
               if (x$termcode!= 1) {
                   dgcv.ubre <- max(abs(x$dgcv.ubre)*max(abs(log(x$sp)),1)/max(abs(x$gcv.ubre),1))
                  cat("\nBFGS termination condition:\n", dgcv.ubre,"\n",sep = "")
               }   
    }
    cat("\n")
    invisible(x)
}



