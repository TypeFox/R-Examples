#' @include getmf.R
NULL

#' generate fixed and random design matrices, X and Z
#'
#' Constructs fixed and random design matrices compromised of either or both growth curve (time-based)
#' components and user-defined nuisance fixed or random effects (input via a formula).  Used to conduct
#' Bayesian mixed effects modeling and to produce growth curve output.
#' 
#' @param time A vector of length{N} (number of subject-time cases) providing times for associated subject-measure observations.  
#'	Identical to \code{time} from \code{\link{dpgrowmm}}.
#' @param trt The treatment group membership vector of length \code{N}.  Assumed numeric and sequential. (e.g. \code{(0,0,1,1,1,2,2,...)}.
#' @param trt.lab Associated labels for the numeric treatment groups.  Each distinct treatment group assumed to have a unique label. 
#' @param subject Vector of length \code{N} providing subject-measure cases.  Identical to \code{subject} from \code{\link{dpgrowmm}}.
#'	\code{P} = number of subjects, \code{q} = number of random effect parameters, per subject.
#' @param n.random A scalar input providing the number of by-subject time-based random effect parameters.
#' @param n.fix_degree The highest polynomial degree to employ for constructing time-based fixed effects covariates. 
#' @param formula A formula of format \code{y ~ x_1 + x_2*x_3 | z_1*z_2 }
#'	where \code{|} separates fixed (to the left of \code{|}) and random effects.
#' @param random.only A boolean scalar used in the case that either fixed or random effects
#'	are entered in \code{formula}, but not both, which case the \code{|} is not entered
#'	(e.g. \code{y ~ x_1 + x_2*x_3 }.  Then, if \code{random.only == TRUE} the variables
#'	on the right-hand side are interpreted to be random effects; otherwise fixed
#'	for use in \code{\link{dpgrow}} and  \code{\link{dpgrowmm}}.
#' @param data Associated data.frame containing names variables in \code{formula}
#' @return A list object containing composed fixed and random effect design matrices, \code{X} and \code{Z}, with column names
#'	   	and subsets: \code{c(X.n, Z.n)} nuisance and \code{c(X.c, Z.c)} growth curve design matrices.  Also returns 
#'		data output, y, if included with formula; otherwise returns y = NULL.
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases XZcov XZ
#' @export
XZcov = function(time = NULL , trt = NULL, trt.lab = NULL, subject = NULL, n.random = NULL, n.fix_degree = NULL, formula = NULL, random.only = NULL, data = NULL)
{
  #################################################################
  ## compose growth curve part of Z and X; Z.c, X.c.
  #################################################################
  ## input check on subject - needed for both growth curve and nuisance covariates
  if( is.null(subject) ) stop("Must input 'subject' vector of length number equal to number of subject-time cases")

  if( !is.null(time) ) ## user wants growthCurve
  {
    ## data input check
    if( is.null(n.random)  ) stop("\nCannot create growth curve random effects design matrix, Z, without 'n.random' input\n")
    ## label unique treatments
    trt.lab 	= unique(trt.lab)
    Nlevel	= length(trt.lab)
    ## Z.c
    Ncase	= length(subject)
    Z.c		= matrix(0, Ncase, n.random)
    for(j in 1:n.random)
    {
	      Z.c[,j]		= time^(j-1)
    }
    colnames(Z.c) 	= paste("time^",0:(n.random-1),sep="")
    colnames(Z.c)[1]	= 1

    if( !is.null(n.fix_degree) ) ## user wants time-based fixed effects
    {
    	## X.c
	## data input check
    	if( is.null(trt)  ) stop("\nCannot create growth curve design matrix, X,  without 'trt' vector input\n")
    	Xc.trt_gen	= matrix(0, Ncase, (n.fix_degree+1))
    	for( k in 1:(n.fix_degree+1) )
    	{
    		Xc.trt_gen[,k]			= time^(k-1)
    	}
    	X.c			= Xc.trt_gen[,-1] ## non-treatment effect fixed covariates (excludes intercept)
    	names.x		= paste("time^",1:n.fix_degree,sep="")

    	## compose trt-based fixed effects
    	if( length(unique(trt)) > 1 ) ## trt is a treatment vector, so add trt*time covariates
    	{
    		for( i in 1:(Nlevel-1) ) ## trt defined to start at 0, so the first level is the baseline hold-out
    		{
			      this.trt	= as.integer( trt == i )
			      X.c		= cbind(X.c,this.trt*Xc.trt_gen)
			      names.x		= c(  names.x, paste("trt",trt.lab[i+1],"*time^",0:n.fix_degree,sep="_") )  ## trt.lab[2] corresponds to i = 1
    		}
    	} ## end conditional statement on whether trt is a treatment vector
    	colnames(X.c)		= names.x
    }else{ ## user wants time-based random effects, but not time-based fixed effects
    	X.c	<- NULL
    }## end conditional statement on whether user wants time-based fixed effects design matrix, X.c
  } ## end conditional statement on whether user wants growth curve
  
  #################################################################
  ## extract nuisance Z.n and X.n from formula and build X, Z
  #################################################################
  if ( !is.null(formula) ) ## user inputs nuisance covariates
  {
    ## check if random.only is not input
    if( is.null(random.only) )
    {
	random.only = FALSE
	warning("Didn't input 'random.only'; function sets it equal to FALSE")
    }
    ## extract model matrix and form nuisance fixed and random design matrices
    cov.res	<- getmf(formula, random.only = random.only, data = data, na.action = na.fail)
    X.n		= cov.res$X
    Z.n		= cov.res$Z
    y 		= cov.res$y ## use may or may not input y in the left-hand side of 'formula'. if not, is.null(y) = TRUE
    if( length(y) != length(subject) ) stop("length of response, y, and 'subject' vectors must both be in subject-time case format")

    ## check for redundant intercept in X.n (first checking that it is non-NULL)
    if(!is.null(X.n) )
    {
	if( all(X.n[,1] == 1) ) X.n = X.n[,-1]

    }

    if( !is.null(time)  ) ## the user wants a growth curve 
    {
    	## check for redundant Z.n, assuming it is non-NULL
    	if( !is.null(Z.n) )
	{
		## check to see if Z.n is identical to Z.c
		if( (nrow(Z.n) == nrow(Z.c)) & (ncol(Z.n) == ncol(Z.c))  )
			{ if( all( (Z.n - Z.c) == 0 ) ) Z.n = NULL }

		## check to see if any columns of Z.n are identical to Z.c
		remove = c()
		for( j in 1:ncol(Z.n) )
		{
			for( k in 1:ncol(Z.c) )
			if( all( (Z.n[,j] - Z.c[,k]) == 0 ) ) {remove = c(remove,j)}
		}
		if( !is.null(remove) ) 
		{
			zn.names		= colnames(Z.n)
			Z.n 			= as.matrix(Z.n[,-remove]) 
			zn.names		= zn.names[-remove]
			colnames(Z.n)	= zn.names
		}
	
	}

    	## check for redundancy between X.c and X.n and remove redundant columns from X.n
	if( !is.null(X.n) & !is.null(X.c) ) ## the user inputs nuisance fixed effects
	{
		remove = c()
    		for( i in 1:ncol(X.c) )
    		{
			for( j in 1:ncol(X.n) )
			{
				if(  all( (X.c[,i] - X.n[,j]) == 0 )  )
				{
					remove = c(remove, j)
				}
			}
    		} ## end loop i on checking for redundant columns between X.n and X.c
		if( !is.null(remove) ) X.n 	= X.n[,-remove] 
	} ## end conditional statement on whether is.null(X.n) 

        ## since one of *.e and *.c exist, combine them into a single X and Z
	X 	= as.matrix(cbind(X.c, X.n))
	Z	= as.matrix(cbind(Z.c, Z.n))

     }else{ ## X.c and Z.c don't exist; user doesn't want a growth curve
	X		= X.n
	Z		= Z.n
	X.c		= NULL
	Z.c		= NULL
     } ## end conditional statment on whether is.null(Z.c)

  }else{ ## is.null(formula) = TRUE, so nuisance covariates created
  	X	= X.c
	Z	= Z.c
	X.n	= NULL
	Z.n	= NULL
	y	= NULL ## set y to NULL since user separately inputs not part of 'data'
  } ## end conditional statement on inclusion of nuisance covariates (formula is not NULL)

  ## reorder objects to subject to ensure that subject is contiguous
  o	<- order(subject)
  X	<- X[o,]
  Z	<- Z[o,]
  if( !is.null(X.c) & !is.null(Z.c) ) 
  {
	X.c	<- X.c[o,]
  	Z.c	<- Z.c[o,]
  }

  if( !is.null(X.n) & !is.null(Z.n) ) 
  {
  	X.n	<- X.n[o,] ## reordering a NULL object still leaves it with value = NULL
  	Z.n	<- Z.n[o,]
  }
  if(!is.null(y)){y <- y[o]}

  ## collinearity check, drop collinear vars if needed
  m <- lm.fit(X, rep(1,nrow(X)), singular.ok=TRUE)
  stopifnot(all(names(m$coef) == colnames(X)))
  v <- which(is.na(m$coef))
  if(length(v) != 0)  ## non-singular
  {        
      	## cat(paste("\nWarning: dropping following collinear fixed covariates:",paste(names(v), collapse=" "), "\n\n"))
	warning(paste("\nDropping following collinear fixed covariates:",paste(names(v), collapse=" "), "\n\n"),call. = FALSE)
      	X <- X[,-as.vector(v)]
  }

  out 	= list(X = X, X.c = X.c, X.n = X.n, Z.n = Z.n, Z = Z, Z.c = Z.c, y = y)
  return(out)

} ## end function XZcov.R