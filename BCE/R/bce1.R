## based on BCE.R, version 1.2
## new in this release: 
## we use normal distributions for A and B, cut off for negative values
## stdev of A and B are absolute values and determined dynamically; a weighting is possible. 
## AX=B ipv X*rat=dat (compatibility with tlsce.r)
## Karel Van den Meersche
## 20081025
################################################################################
## work in progress
################################################################################

bce1 <- function(
                ## parameters
                A,
                B,
                Wa=NULL,
                Wb=NULL,
                jmpType  = "default",
                jmpA     = 0.1, 
                jmpX     = 0.1, 
                jmpCovar = NULL,
                initX    = NULL,
                initA   = NULL, 
                priorA  = "normal",
                minA    = NULL, 
                maxA    = NULL, 
                var0    = NULL, 
                wvar0   = 1e-6,    
                Xratios = TRUE, 
                verbose=TRUE,
                ...             
                )
  {

    ##=============================##
    ## warnings and error messages
    ##=============================##

    if (NROW(A)!=NROW(B)) stop("A and B must have same number of rows")

    
    ##===============================================##
    ## general initialisations
    ##===============================================##

    A0 <- as.matrix(A)
    B0 <- as.matrix(B)

    ## useful numbers & names
    nalg <- ncol(A0)
    nst  <- ncol(B0)
    npig <- nrow(A0)
    algnames   <- colnames(A0)       # number & names of taxonomic groups
    stnames    <- colnames(B0)       # number & names of stations or samples
    pignames   <- rownames(A0)       # number & names of biomarkers
    lr <- npig*nalg
    if (is.null(algnames)) algnames <- 1:nalg
    if (is.null(stnames)) stnames <- 1:nst
    if (is.null(pignames)) pignames <- 1:npig
    lc <- nalg*nst
    w <- A0>0
    lw <- length(which(w))
    
    if (is.null(Wa)) Wa <- matrix(1,npig,nalg) else
    if (length(Wa)==1) Wa <- matrix(Wa,npig,nalg) else Wa <- as.matrix(Wa)
    if (is.null(Wb)) Wb <- matrix(1,npig,nst) else
    if (length(Wb)==1) Wb  <- matrix(Wb,npig,nst) else Wb <- as.matrix(Wb)
    if (is.null(minA)) minA <- matrix(0,npig,nalg)
    if (is.null(maxA)) maxA <- matrix(Inf,npig,nalg)
    
    if (is.null(initX)|is.null(var0)|jmpType=="covar")
      tlsce0 <- tlsce(A0,B0,Wa,Wb,minA=minA,maxA=maxA,Xratios=Xratios)
    
    if (is.null(initX))  X0 <- tlsce0$X  else  X0 <- initX
    
    ##     for (i in 1:npig) for (j in 1:nst) 
    ##       if (B[i,j]==0 & Wb[i,j]==Inf) X0[A[i,]!=0,j] <- 0
    
    ## if Xratios: X = ZQ+P; EX=F (rowsums of X are 1); P=c(0,..,0,1)
    ## Q are the estimated model parameters if Xratios
    if (Xratios)
      {
        suppressWarnings(Z <- matrix(c(1,-1,rep(0,nalg-1)),nalg,nalg-1))
        P <- c(rep(0,nalg-1),1)
        Zinv <- matrix(1,nalg-1,nalg); Zinv[upper.tri(Zinv)] <- 0
        Q0 <- Zinv%*%X0
      } 
    
    
    ##==================================================##
    ## initialisation parameters and residuals function ##
    ##==================================================##

    if (is.null(initA)) initA <- A 
    if (Xratios) par <- c(initA[w],Q0) else par <- c(initA[w],X0)
    lp <- length(par)
    names(par)[1:lw] <-
      paste("A",
            abbreviate(rownames(A)[which(w,arr.ind=TRUE)[,1]]),
            abbreviate(colnames(A)[which(w,arr.ind=TRUE)[,2]]),
            sep="_")
    names(par)[-(1:lw)] <-
      paste("Q",
            1:(nalg-1),
            rep(abbreviate(colnames(B)),each=nalg-1),
            sep="_")
    if (is.null(var0)) var0 <- tlsce0$SS["total"]/lp
    
    ## residuals function
    residuals <- function(par,...)
      {
        A1 <- A0; A1[w] <- par[1:lw]
        Q <- matrix(par[-(1:lw)],ncol=nst)
        if (Xratios) X <- Z%*%Q+P else X <- Q
        AX <- A1%*%X
        if (priorA=="normal") resid <- c(Wa[w]*(A1[w]-A0[w]),Wb*(AX-B0)) else {
          if (priorA=="uniform") resid <- c(Wb*(AX-B0)) }
        return(resid)
      }

    
    ## prior information: all elements of A and X are positive; A is limited by minA, maxA
    lowerpar <- c(minA[w],rep(-Inf,lp-lw))
    upperpar <- c(maxA[w],rep(+Inf,lp-lw))
    if (Xratios)
      {
        prior <- function(par)
          {
            Q <- matrix(par[-(1:lw)],ncol=nst)
            X <- Z%*%Q+P
            ifelse(any(X<0),-Inf,0)
          }
      } else {
        prior <- function(par) ifelse(any(par[-(1:lw)]<0),-Inf,0)
      }


    ##=============================##
    ## initialisation jump lengths ##
    ##=============================##

    if (jmpType=="default")
      {
        if (!length(jmpA)%in%c(1,lr)) stop("The jump length of the ratio matrix should be either a single value, or a matrix with the same dimension as A")
        if (any(is.na(jmpA))) stop("missing values in jump ratio matrix; please specify a valid jump ratio matrix.")
        if (length(jmpA)==1) jmpA <- rep(jmpA,lw) else jmpA <- jmpA[w]
        if (length(jmpX)==1) jmpX <- matrix(jmpX,nalg,nst)
        if (Xratios) 
          {
            jmpQ <- (jmpX[-1,]+jmpX[-nalg,])*.5
          } else {
            jmpQ <- jmpX
          }
        jmp <- c(jmpA,jmpQ)
      }

    if (jmpType=="estimate")
      {
        if (length(jmpA)==1)
          jmpA <- jmpA*summary(tlsce0$fit)$cov.scaled*2.4^2/(lp-1)
        if (length(jmpX)==1)
          {
            jmpQ <- matrix(0,0,0)
            for (i in 1:nst)
              {
                BnotNA <- !is.na(B0[,i])  # remove NA from B
                if (Xratios) {
                  Qlseii <- lsei(A[BnotNA,]%*%Z,B0[BnotNA,i]-A[BnotNA,]%*%P,E=rep(0,nalg-1),F=0,G=Z,H=-P,Wa=Wb[BnotNA,i],fulloutput=TRUE)
                } else {
                  Qlseii <- lsei(A[BnotNA,],B0[BnotNA,i],E=rep(0,nalg),F=0,G=diag(nalg),H=rep(0,nalg),Wa=Wb[BnotNA,i],fulloutput=TRUE)
                }
                jmpQ <- bdiag(jmpQ,Qlseii$covar*2.4^2/(lp-1))
              }
            jmp <- as.matrix(bdiag(jmpA,jmpX*jmpQ))
          }
      }

    if (jmpType=="covar") jmp <- jmpCovar

        

    ##==========##
    ## mcmc
    ##==========##

    mcmc <- modMCMC(f=residuals,p=par,var0=var0,wvar0=wvar0,prior=prior,jump=jmp,lower=lowerpar,upper=upperpar,...)


    ##============##
    ## output     ##
    ##============##

    ##     outputlength <- nrow(mcmc$pars)
    ##     mcmc.A <- array(0,dim=c(npig,nalg,outputlength),dimnames=list(pignames,algnames,NULL))
    ##     mcmc.A[w] <- t(mcmc$pars[,1:lw]) ## check this!!
    ##     mcmc.X <- array(dim=c(nalg,nst,outputlength),dimnames=list(algnames,stnames,NULL))
    ##     if (Xratios)
    ##       {
    ##         mcmc.Q <- array(dim=c(nalg-1,nst,outputlength),dimnames=list(NULL,stnames,NULL))
    ##         mcmc.Q[] <- t(mcmc$pars[,-(1:lw)])
    ##         for (i in 1:outputlength) mcmc.X[,,i] <- Z%*%mcmc.Q[,,i]+P
    ##       } else {
    ##         mcmc.X[] <- t(mcmc$pars[,-(1:lw)])
    ##       }

    ## value <- list(A=mcmc.A,X=mcmc.X,mcmc=mcmc)
    value <- mcmc

    class(value) <- c("bce","modMCMC")
    attr(value,"A_not_null") <- w
    attr(value,"Xratios") <- Xratios
    attr(value,"pignames") <- pignames
    attr(value,"algnames") <- algnames
    attr(value,"stnames") <- stnames
    
    return(value)
  }
