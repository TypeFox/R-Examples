Mort2Dsmooth_checker <-
function(x, y, Z, offset, W,
                                 overdispersion,
                                 ndx, deg, pord,
                                 lambdas, df, method,
                                 coefstart,
                                 control){
  ## Input:
  ## x: abcissa of data
  ## y: ordinate of data
  ## Z: matrix of count response
  ## offset: an a priori known component (optional)
  ## overdispersion: logical on the
  ##                 presence of overdispersion
  ## W: a matrix of weights to be used
  ##    in the fitting process 
  ## ndx: a vector for the numbers
  ##      of internal knots -1 for both axes. 
  ## deg: a vector for the degrees of the B-splines for
  ##      the x-axis and y-axis. 
  ## pord: a vector for the order of
  ##       differences for both axes. 
  
  ## lambdas: a vector of smoothing parameters
  ##          for both axes (optional)
  ## df: degree of freedom for both axes (optional)
  ## method: the method for controlling
  ##         the amount of smoothing.
  ## coefstart: eventual initial coefficients
  ## control: a list of control parameters
  
  ## Output: a list containing CHECKED arguments
  ##         for the Mort2Dsmooth function
  
  ## about x and y:
  if(missing(x)){
    x <- 1:nrow(Z)
  }
  if(missing(y)){
    y <- 1:ncol(Z)
  }
  ## about the offset
  offsetINIT <- offset
  m <- length(x)
  n <- length(y)
  ## about infinitive offset, i.e. zero exposures 
  whioff <- which(is.infinite(offset))
  whiwei <- which(W==0)
  if(any(!whioff%in%whiwei)){
    stop("weights different from zero associated with infinitive offset values")
  }
  offset[c(whioff, whiwei)] <- 100

  ## about lengths and wrong values
  if (length(x)!=nrow(Z)) 
    stop("length of x must be equal to number of rows in Z")
  if (length(y)!=ncol(Z)) 
    stop("length of y must be equal to number of columns in Z")
  if (dim(Z)[1] != m | dim(Z)[2] != n)
    stop("Argument arrays of wrong length")
  if (dim(offset)[1] != m | dim(offset)[2] != n)
    stop("Argument arrays of wrong length")
  if (deg[1]<1 | deg[1]>=10 | deg[2]<1 | deg[2]>=10) 
    stop("Wrong values for deg")    
  if (pord[1]<=0 | pord[1]>=5 | pord[2]<=0 | pord[2]>=5) 
    stop("Wrong value for pord")
  if (ndx[1]<2 | ndx[1]>=floor(m*.9) |
      ndx[2]<2 | ndx[2]>=floor(n*.8))
    stop("Wrong value for ndx")
  coefstart.check <- is.null(coefstart)
  if(!coefstart.check){
    if(nrow(coefstart)!=(ndx[1]+deg[1]) |
       ncol(coefstart)!=(ndx[2]+deg[2])){
      stop("coefstart must be a ndx[1]+deg[1] times ndx[2]+deg[2] matrix")
      }
  }
  ## about method
  if (method != 1 & method != 2 & method != 3 & method != 4) 
    stop("Wrong value for method")
  ## method = 1 adjusts lambda so that the BIC is minimized.
  ## method = 2 adjusts lambda so that the AIC is minimized.
  ## method = 3 uses the value supplied for lambdas. 
  ## method = 4 adjusts lambdas so that the degrees of freedom is
  ##          equal to df (isotopic smoothing).

  ## check-point methods
  lambdas.check <- is.null(lambdas)
  df.check <- is.null(df)
  MET <- NULL
  ## both lambdas and df NULL
  if(lambdas.check & df.check & method==1){MET=1}
  if(lambdas.check & df.check & method==2){MET=2}
  if(lambdas.check & df.check & method==3){
    stop("with method 3, provide lambdas")
  }
  if(lambdas.check & df.check & method==4){
    stop("with method 4, provide df")
  }
  ## lambdas NULL and df GIVEN
  if(lambdas.check & !df.check & method==1){
    stop("df and method 1 cannot be chosen together")
  }
  if(lambdas.check & !df.check & method==2){
    stop("df and method 2 cannot be chosen together")
  }
  if(lambdas.check & !df.check & method==3){
    stop("df and method 3 cannot be chosen together")
  }
  if(lambdas.check & !df.check & method==4){
    MET=4
    warning("Isotropic smoothing is applied", call.=FALSE)
  }
  ## lambdas GIVEN and df NULL
  if(!lambdas.check & df.check & method==1){
    stop("lambdas and method 1 cannot be chosen together")
  }
  if(!lambdas.check & df.check & method==2){
    stop("lambdas and method 2 cannot be chosen together")
  }
  if(!lambdas.check & df.check & method==3){MET=3}
  if(!lambdas.check & df.check & method==4){
    stop("lambdas and method 4 cannot be chosen together")
  }
  ## both lambdas and df GIVEN, regardless method
  if(!lambdas.check & !df.check){
    stop("lambdas and df cannot be chosen together")
  }  
  ## impossible values for lambdas and df
  if(!lambdas.check && lambdas<0)
    stop("lambdas must be positive")
  if(!df.check && df < (pord[1] + pord[2]))
    stop("df must be larger than the sum of pord")
  if(!df.check && df > ((ndx[1]+deg[1])*(ndx[2]+deg[2])))
    stop("df must be smaller than the product of (ndx+deg) values")
  if(!df.check & length(df)!=1)
    stop("df must be length 1")
  if(!lambdas.check &length(lambdas)!=1 &length(lambdas)!=2)
    stop("lambda must be length 1 or 2")
  if (!lambdas.check & length(lambdas)==1){
    lambdas <- rep(lambdas, 2)
    warning("Isotropic smoothing is applied", call.=FALSE)
  }
  ## setting control-parameters
  con <- list(MON=FALSE, TOL1=1e-06, TOL2=0.5,
              RANGEx=c(10^-4, 10^6), RANGEy=c(10^-4, 10^6), 
              MAX.IT=50)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  ## stop about weights
  if(nrow(W)!=m | ncol(W)!=n){
    stop("dimensions of W and Z must be equal")
  }
  ## warning about interpolation/extrapolation
  if(any(W==0)){
    warning("Interpolation and/or extrapolation is taking place",
            call. = FALSE)
  }
  ## about the overdispersion parameter
  if(overdispersion & method==3)
    warning("given method 3, overdispersion is computed a posteriori")
  if(overdispersion & method==4)
    warning("given method 4, overdispersion is computed a posteriori")
  ## warning about weights
  if(min(W) < 0) {
    warning(paste("At least one weight entry is negative"))
  }
  ## returning
  llist <- list(x=x, y=y, Z=Z, offset=offset, W=W, m=m, n=n,
                offsetINIT=offsetINIT,
                overdispersion=overdispersion,
                ndx=ndx, deg=deg, pord=pord,
                lambdas=lambdas, df=df, method=method,
                coefstart=coefstart,
                control=con)
  llist
}
