################################################################
## Total Least Squares Composition Estimator
## use modFit
## this is an orthogonal alternative to chemtax
## Wb added (is called Wa in lsei...)
## is tested with all 4 examples (bce-tests.r) and performs
## as well or better than the previous tlsce.
## the use of modFit also allows for more flexibility in tuning
## the optimization algorithm, and more output details (hessian,...)
################################################################

tlsce <- function(A,
                  B,
                  Wa=NULL,
                  Wb=NULL,
                  minA=NULL,
                  maxA=NULL,
                  A_init=A,
                  Xratios=TRUE,
                  ## chemtax=FALSE,
                  ...)         
  {

    ##=================##
    ## initialisations ##
    ##=================##

    A <- as.matrix(A)
    B <- as.matrix(B)
    
    l <- nrow(A)                        # number of pigments
    m <- ncol(A)                        # number of species
    n <- NCOL(B)                        # number of samples
    w <- which(A>0)
    lw <- length(w)
    A_c <- A[w]                         # non-zero elements of A
    if (Xratios#|chemtax
        ) {
      E <- t(rep(1,m)); F <- t(rep(1,n))} else {
        E <- t(rep(0,m)); F <- t(rep(0,n))} # sum of species fractions is 1 or not
    G <- diag(1,m); H <- matrix(0,m,n)  # all elements positive
    if (is.null(Wa)) Wa <- 1            # weighting of elements of A and B
    if (length(Wa)==1) Wa <- matrix(Wa,l,m)
    if (length(Wa)==length(A)) Wa_c <- Wa[w]
    if (length(Wb)==1) Wb <- matrix(Wb,l,n)
    A_c_init <- A_init[w]
    if (is.null(minA)) minA_c <- rep(0,lw) else minA_c <- minA[w]
    if (is.null(maxA)) maxA_c <- rep(+Inf,lw) else maxA_c <- maxA[w]
    
    ##     if (chemtax)
    ##       {
    ##         A_rescaled <- rbind(A,1)
    ##         A_rescaled <- A_rescaled/rep(colSums(A_rescaled),each=l+1)
    ##         B_rescaled <- rbind(B,1)
    ##         B_rescaled <- B_rescaled/rep(colSums(B_rescaled),each=l+1)
    ##         if (!is.null(Wb)) Wb_rescaled <- rbind(Wb,1)
    ##         Wa_rescaled <- rbind(Wa,colMeans(Wa)) ## ad hoc solution...
    
    ##         residuals <- function(A_c_new)
    ##           {
    ##             A_new <- A
    ##             A_new[w] <- A_c_new
    ##             A_new_rescaled <- rbind(A_new,1)
    ##             A_new_rescaled <- A_new_rescaled/rep(colSums(A_new_rescaled),each=l+1)
    ##             X <- LSEI(A_new_rescaled,B_rescaled,E,F,G,H,Wa=Wb)$X
    ##             if (is.null(Wb)) return(c(Wa_rescaled*(A_new_rescaled-A_rescaled),A_new_rescaled%*%X-B_rescaled))
    ##             return(c(Wa_rescaled*(A_new_rescaled-A_rescaled),Wb_rescaled*(A_new_rescaled%*%X-B_rescaled)))
    ##           }
    ##       } else {
    residuals <- function(A_c_new)
      {
        A_new <- A
        A_new[w] <- A_c_new
        X <- LSEI(A_new,B,E,F,G,H,Wa=Wb)$X
        if (is.null(Wb)) return(c(Wa_c*(A_c-A_c_new),A_new%*%X-B))
        return(c(Wa_c*(A_c-A_c_new),Wb*(A_new%*%X-B)))
      }
    ##      }

    
    ##===========##
    ## model fit ##
    ##===========##
    
    tlsce_fit <- modFit(residuals,A_c,lower=minA_c,upper=maxA_c,...)

    
    ##========##
    ## output ##
    ##========##
    
    A_c_fit <- tlsce_fit$par
    A_fit <- A; A_fit[w] <- A_c_fit
##     if (chemtax)
##       {
##         A_fit_rescaled <- rbind(A_fit,1)
##         A_fit_rescaled <- A_fit_rescaled/rep(colSums(A_fit_rescaled),each=l+1)
##         LSEI_fit <- LSEI(A_fit_rescaled,B_rescaled,E,F,G,H,Wa=Wb)
##         X <- LSEI_fit$X; rownames(X) <- colnames(A); colnames(X) <- colnames(B)
##         B_fit_rescaled <- A_fit_rescaled%*%X
##         B_fit <- B_fit_rescaled[-(l+1),]/rep(B_fit_rescaled[l+1,],each=l)
##       }
##     else
##       {
        LSEI_fit <- LSEI(A_fit,B,E,F,G,H,Wa=Wb)
        X <- LSEI_fit$X; rownames(X) <- colnames(A); colnames(X) <- colnames(B)
        B_fit <- A_fit%*%X
##      }
    ssr <- tlsce_fit$ssr
    ssr_B <- LSEI_fit$solutionNorm
    ssr_A <- ssr-ssr_B
    solutionNorms <- c(ssr,ssr_A,ssr_B); names(solutionNorms) <- c("total","A","B")

    return(list(X=X,
                A_fit=A_fit,
                B_fit=B_fit,    # the fits
                SS=solutionNorms, # residual sums of squares
                fit=tlsce_fit)) # a modFit object
  }


##############################################################
## helper functions
##############################################################


LSEI <- function(A=NULL,B=NULL,E=NULL,F=NULL,G=NULL,H=NULL,Wa=NULL,...)
  {
    if (is.vector(B)) return(lsei(A,B,E,F,G,H,Wa=Wa,...))
    else
      {
        X <- matrix(NA,ncol(A),ncol(B))
        solutionNorm <- 0
        for (i in 1:ncol(B))
          {
            BnotNA <- !is.na(B[,i])  # remove NA from B
            ls <- lsei(A[BnotNA,],B[BnotNA,i],E,F[,i],G,H[,i],Wa=Wa[BnotNA,i],...)
            X[,i] <- ls$X
            solutionNorm <- solutionNorm + ls$solutionNorm
          }
        return(list(X=X,solutionNorm=solutionNorm))
      }
  } # LSEI


## LSEI <- function(A=NULL,B=NULL,E=NULL,F=NULL,G=NULL,H=NULL,Wa=NULL,...)
##   {
##     if (is.vector(B)) B <- as.matrix(B)

##     X <- matrix(NA,ncol(A),ncol(B))
##     solutionNorm <- 0
##     for (j in 1:ncol(B))
##       {
        
##         BnotNA <- !is.na(B[,j])  # remove NA from B (missing data)
##         Xpresent <- colSums(subset(A,B[,j]==0))==0  # remove missing groups from A (biomarker not found)
##         X[!Xpresent,j] <- 0
##         ls <- lsei(A[BnotNA,Xpresent],B[BnotNA,j],E[,Xpresent],F[,j],G[,Xpresent],H[,j],Wa=Wa[BnotNA,j],...)
##         X[Xpresent,j] <- ls$X
##         solutionNorm <- solutionNorm + ls$solutionNorm
##       }
##     return(list(X=X,solutionNorm=solutionNorm))
##   } # LSEI

