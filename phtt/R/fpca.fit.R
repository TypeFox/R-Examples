fsvd.pca <- function(Q,
                            allow.dual      = TRUE,
                            given.d         = NULL,
                            calcul.loadings = TRUE,
                            neglect.neg.ev  = FALSE, 
				    spar = NULL){
  ## extract data information
  nr      <- nrow(Q)
  nc      <- ncol(Q)

  ## save original Q-values 
  Q.non.smth <- Q 

  ## smoothing Q (small degree of undersmoothing) ===============================================#
  if(is.null(spar)) {													 #
  spar.low <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, method = 4       )$spar  * 0.75 }  #
  else spar.low = spar                                                                           #
  Q        <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, spar   = spar.low)$ysmth           #
  ##=============================================================================================#
  
  ## data informations
  beding = (nr>nc && calcul.loadings && allow.dual)
  dual <- ifelse(beding, TRUE, FALSE)
  if(dual){
    Q <- t(Q)
  }
  
  ## Compute spectral decompostion
  cov.mat         <- tcrossprod(Q)
  Spdec           <- eigen(cov.mat, symmetric= TRUE)
  Eval	          <- Spdec[[1]]
  Evec    	  <- Spdec[[2]]

  ## Compare rank and given.d
  
  nbr.pos.ev	   <- length(Eval[Eval > 0])
  max.rk 	   <- ifelse(neglect.neg.ev, nbr.pos.ev, length(Eval))
  if(is.null(given.d)){
    given.d <- max.rk
  }else{
    if(given.d > max.rk){
      warning(c("The given dimension 'given.d' is larger than the number of positve eigen values."))
    }
    given.d <- min(given.d, max.rk)
  }

  ## compute spectral variance decomposition

  ## Left side decomposion
  L                         <- Evec[,0:max.rk , drop= FALSE]   # dual==FALSE: (T x max.rk), dual==TRUE: (N x max.rk)

  ## sqrt of eigenvalues      
  if(!neglect.neg.ev){
    sqr.E                   <- c(sqrt(Eval[Eval > 0]), rep(0, (max.rk - nbr.pos.ev)))
  }else{
    sqr.E	            <- sqrt(Eval[1:max.rk ])
  }

  ## computation of the loadings-parameter (scores)
  if(calcul.loadings){
    S                       <- crossprod(Q, L)[, 0:max.rk , drop= FALSE] # crossprod: t(Q) %*% L
                                                                         # if dual: dim(S)= TxN, if non-dual: dim(S)=NxT

    R                       <- S %*% diag(diag(crossprod(S))^-{0.5}, max.rk) 


    ## no pca-fitting
    if(((given.d==max.rk)&&!neglect.neg.ev)){
      Q.fit                 <- Q
    ## pca-fitting  
    }else{
      Q.fit                 <- tcrossprod(L[, 0:given.d , drop= FALSE],
                                          L[, 0:given.d , drop= FALSE]) %*% Q
    }
  ## no computation of the loadings-parameter (scores)  
  }else{
    R                   <- NULL
    ## no pca-fitting
    if(((given.d==max.rk)&&!neglect.neg.ev)){
      Q.fit <- Q
    ## pca-fitting
    }else{
      Q.fit                 <- tcrossprod(L[, 0:given.d , drop= FALSE],
                                          L[, 0:given.d , drop= FALSE]) %*% Q
    }
  }

  ## re-convert dimension if dual covariance matrix was used
  if(dual){
    u      <- L
    L      <- R
    R      <- u
    Q.fit  <- t(Q.fit)
    Q      <- t(Q)          # (under)smoothed Q   
  }
  
  ## prepare return-values       
  d.seq <- seq.int(0, (max.rk-1)) # dimension-sequence
  E     <- Eval[1:max.rk]
  sum.e <- sum(E)
  cum.e <- cumsum(E)
  V.d   <- c(sum.e, sum.e-cum.e[-length(cum.e)])

  ## return
  structure(list(L              = L,
                 R              = R,
                 Q.orig         = Q.non.smth,
                 Q.orig.smth    = Q,
                 Q.fit          = Q.fit,
                 spar.low       = spar.low,
                 E              = E,
                 sqr.E          = sqr.E,
                 given.d        = given.d,
                 d.seq          = d.seq,
                 V.d            = V.d,
                 nr             = nr,
                 nc             = nc ,
                 cov.mat        = cov.mat,
                 dual           = dual),
            class   = "fsvd.pca")
}






############################
# Main function: fpac.fit  #
# Descrition:
# The function calculates functional spectral decomposition of a matrix in >>dat<<
# dimension: dim(dat) == T x N
# and, according to the argument >>given.d<<, does pca-fitting. Generally it is assumend
# that the column of the TxN-matrix >>dat<< are demaened by the "mean-column".
# Calls:                   #
# is.regular.panel()       #
# fsvc.pca()               #
# restrict.pca()           #
# Takes:===========================================================================================#
# dat                      (smoothed raw-data for fPCA)                                            #
# given.d = NULL           (user given dimension)                                                  #
# restrict.mode            ("restrict.factors" or "restrict.loadings")                             #
# allow.dual = TRUE        (possibility to switch of the dual-matrix calculations)                 #
# neglect.neg.ev = TRUE    (if TRUE:  max.rk = number of pos. eigenvalues                          #
#                           if FALSE: max.rk = number of pos. and evtl. numerical neg. eigenvalues)#
# Gives:===========================================================================================#
# factors
# loadings         
# fitted.values    
# orig.values.smth
# orig.values
# spar.low
# cov.matrix      
# eigen.values    
# Sd2              
# given.d        
# data.dim   
# dual       
# L          
####################################################################################################

fpca.fit <- function(dat,
                     given.d        = NULL,
                     restrict.mode  = c("restrict.factors","restrict.loadings"),
                     allow.dual     = TRUE,
                     neglect.neg.ev = TRUE, 
			   spar = NULL){


  ## Check input
  is.regular.panel(dat, stopper = TRUE)

  ## fPCA
  fpca.obj       <- fsvd.pca(Q              = dat,
                             given.d        = given.d,
                             allow.dual     = allow.dual,
                             neglect.neg.ev = neglect.neg.ev, 
				     spar = spar)
                           
  
  ## impose Restrictions (default: restrict.factors such that F'F/T = I )
  result        <- restrict.pca(fpca.obj)

  ## return
  structure(result, class = "fpca.fit")
}
