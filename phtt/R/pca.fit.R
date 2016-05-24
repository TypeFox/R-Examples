######## check the input data if is a regular panel ########################
# is.regular.panel check if the data is a matrix or dataframe 
# without NA NaN
# Input:
#	dat: input data
#	comment: if TRUE a appropriate message is given if the class of dat 
#		is not appropriate
#	stopper: if TRUE the function will be break 
# Output:
#	a lofical value TRUE or FALSE

#############################################################################
is.regular.panel <- function(dat,
                             comment = FALSE,
                             stopper = FALSE){
  
  if(!is.matrix(dat)&!is.data.frame(dat)){
    if(stopper) stop(expression("Each Variable of your model should be in a Txn-matrix 
				form or has a 'data.frame' class"),  call. = FALSE)
    if(comment) message(expression("Each Variable of your model should be in a Txn-matrix 
				form or has a 'data.frame' class"))
    FALSE
  }
  else{
    if(!is.numeric(dat)){
      if(stopper) stop(expression("The data is not numeric")
		,  call. = FALSE)
      if(comment) message(expression("The data is not numeric"))
      FALSE
    }
    else{
      if(any(is.na(dat))| any(is.nan(dat))){
        if(stopper) stop(expression("The data has NA or NaN values")
			,  call. = FALSE)
        if(comment) message(expression("The data has NA or NaN values"))
        FALSE
      }
      else TRUE
    }
  }
}

########### spectral variance decomposition and pca fitting #################
# Input:
#  Q 			 : an arbitrary n times p matrix (not ncessary squar)
#  allow.dual	 : if TRUE allow to calculate the eigenvalues and the 
#				eigenvectos from the dual matrix Q'Q instead of QQ'
#				if the dual has a smaller dimension				
#  given.d 		 : dimension of the orthogonal compnents
#  calcul.loadings : logical values specifies whether the loading parmetrs 
#				have to be calculated or not (default is TRUE);
#				an option to speed up computations in special cases
#  neglect.neg.ev  : logical values specifies whether the negative 
#				eigenvalues have to be negelcted of 
#				not  (default is FALSE)
# Output:
#  L : left side decomposion where L'L = I
#  R=: left side decomposion where R'R = I
#  Q.fit: the fitted values according to the given dimension 'given.d'
#  Q: the inpude data matrix (reproduced)
#  E: the eingenvalues of the non square matrix QQ'
#  sqr.E: the squared eigenvalues sqrt(E)
#  given.d: the given dimension to be used for pca
#  d.seq: a sequence of from 0 to maximal rank of QQ'
#  V.d : the cumulated sume of the remaining eingenvalues after the pca fit
#  nr: number of rows of Q
#  nc: number of colomns of Q
#  cov.mat: the non scaled covariance matrix QQ' 
#  dual: indicates whether the eigenvectors are from Q'Q or QQ' calculated
#  class = "svd.pca"
#############################################################################
svd.pca <- function(Q
		, given.d=NULL
		, calcul.loadings = TRUE
		, allow.dual = TRUE
		, neglect.neg.ev = FALSE){

  # data informations

	nr   <- nrow(Q)
	nc   <- ncol(Q)
	beding = (nr>nc && calcul.loadings && allow.dual)
	dual <- ifelse(beding, TRUE, FALSE)
	if(dual) Q <- t(Q)

  # Compute spectral decomposion 

	cov.mat	<- tcrossprod(Q)
	Spdec		<- eigen(cov.mat, symmetric= TRUE)
	Eval		<- Spdec[[1]]
	Evec	  	<- Spdec[[2]]

  # compare rank and given.d

	nbr.pos.ev	   <- length(Eval[Eval > 0])
	max.rk 	   <- ifelse(neglect.neg.ev, nbr.pos.ev, length(Eval))
	if(is.null(given.d)) given.d <- max.rk
	else {
		if(given.d > max.rk) warning(c("The given dimension 'given.d' 
			is larger than the number of positve eigen values."))
		given.d <- min(given.d, max.rk)
		}

  # compute spectral variance decomposition

	L      <- Evec[, 0:max.rk , drop= FALSE]
	if(!neglect.neg.ev) sqr.E <- c(sqrt(Eval[Eval > 0])
						,rep(0, (max.rk - nbr.pos.ev)))
	else sqr.E	 <- sqrt(Eval[0:max.rk ])

	if(calcul.loadings){
		S <- crossprod(Q, L)[, 0:max.rk , drop= FALSE]
		R <- S %*% diag(diag(crossprod(S))^-{0.5}, max.rk)
		if(((given.d==max.rk)&&!neglect.neg.ev)) Q.fit <- Q
		else Q.fit  <- tcrossprod(L[, 0:given.d , drop= FALSE]
				, L[, 0:given.d , drop= FALSE])%*%Q
		}

	else{
		R <- NULL
		if(((given.d==max.rk)&&!neglect.neg.ev)) Q.fit <- Q
		else Q.fit  <- tcrossprod(L[, 0:given.d , drop= FALSE]
				, L[, 0:given.d , drop= FALSE])%*%Q
		}

  # convert dimension if dual covariance matrix is  used

	if(dual){
		u <- L
		L <- R
		R <- u
		Q.fit <- t(Q.fit)
                # damit die daten bei der rückgabe die gleichen 
		    #	dimensionen haben:
                Q     <- t(Q)
		}

  # about dimension
	d.seq <- seq.int(0, (max.rk-1))
	E 	<- Eval[1:max.rk]
	sum.e <- sum(E)
	cum.e <- cumsum(E)
	V.d   <- c(sum.e, sum.e-cum.e[-length(cum.e)])

	structure(list(L=L, R=R, Q.fit=Q.fit, Q=Q, E=E
		, sqr.E=sqr.E, given.d=given.d
		, d.seq=d.seq, V.d=V.d, nr=nr, nc=nc 
		,cov.mat=cov.mat, dual=dual), class = "svd.pca")
}
###### test 
#test <- matrix(c(1,2,0, 3), 2,2 )
#svd.pca(test, given.d = 0)

########################## restrict.pca ####################################
#Input: 

# 	restrict mode: "restrict.factors" leads to F'F/T = I 
#	and "restrict.loadings" to Lamd'Lamd/N = I 
# 	 
#
#Output:
#	factors,	  : the factors after restrictions
#	loadings      : loadings parameters after restrictions,
#	fitted.values : the fitted values,
#	orig.values   : original data used as input svd.pca function,
#	cov.matrix    : the covariance matrix,
#	eigen.values  : Eigen values of the cov.natrix
#	Sd2           : the variance of the remaining residuals 
#			    corresponding to each dimension,
#	given.d       : the given or the used factor dimension 
#	data.dim      : c(nr, nc) where nr = number of rows 
#			    and nc = number of colomns
#	dual          : logical value informing us whether we used
#			     by the calculation the dual covariance matrix 
#			     or the classical covariance matrix
#     L             : the unrestricted factors 
############################################################################
restrict.pca <- function(obj
		, restrict.mode= c("restrict.factors","restrict.loadings")){

  # check the object class and it conformitiy

	if(class(obj) !="svd.pca"&&class(obj)!="fsvd.pca"){stop(c("Argument 
		has to be either a 'svd.pca' or a 'fsvd.pca' object."))}
	if(is.null(obj$R)){  stop(c("Loading-parameters are missing."))}

  # restrict variance and eigenvalues
	given.d	<- obj$given.d
	L		<- obj$L[, 0:given.d , drop= FALSE]
	R		<- obj$R[, 0:given.d , drop= FALSE]
	sqr.E		<- obj$sqr.E
	fitted.values	<- obj$Q.fit
	nr		<- nrow(fitted.values)
	nc		<- ncol(fitted.values)
	cov.matrix	<- obj$cov.mat/(nr*nc)
	Sd2		<- obj$V.d/(nr*nc)
	Eval		<- obj$E/(nr*nc)

 # restric factors and loadings

	re.mo <-match.arg(restrict.mode)
	switch(re.mo,
			restrict.factors={
				factors   <- L*sqrt(nr)
				loadings  <- R%*%diag(sqr.E[0:given.d]
						, given.d)/sqrt(nr)
			},
			restrict.loadings ={
				factors   <- L%*%diag(sqr.E[0:given.d]
						, given.d)/sqrt(nc)
				loadings  <- R*sqrt(nc)
			}
		)

  ## prepare return-list

	result <- list(cov.matrix    = cov.matrix,
			eigen.values  = Eval,
			Sd2           = Sd2,
			given.d       = given.d,
			L             = L,
			data.dim      = c(nr, nc),
			dual          = obj$dual,
			loadings      = loadings,
			factors       = factors,
			fitted.values = fitted.values)


  ## extend list according to the object class (svd.pca or fsvd.pca)

	if(class(obj) == "fsvd.pca"){
		result$orig.values = obj$Q.orig
		result$orig.values.smth = obj$Q.orig.smth
		result$spar.low		= obj$spar.low
	}
	else{
		result$orig.values = obj$Q
	}

  ## return list

	result
}


############################################################################
## The rapper function 
pca.fit <- function(dat, given.d=NULL 
		, restrict.mode= c("restrict.factors","restrict.loadings")
		, allow.dual = TRUE, neglect.neg.ev = TRUE){
	is.regular.panel(dat, stopper = TRUE)
	svd.pca.obj  <- svd.pca(dat, given.d = given.d
				, allow.dual= allow.dual
				, neglect.neg.ev = neglect.neg.ev)
	restrict.mode <-match.arg(restrict.mode)
	result <- restrict.pca(svd.pca.obj, restrict.mode= restrict.mode)
	structure(result, class = "pca.fit")
	}
############################################################################
