Mort1Dsmooth_checker <-
function(x, y, offset, w,
                                 overdispersion,
                                 ndx, deg, pord, 
                                 lambda, df, method,
                                 coefstart,
                                 control){
  ## Input:
  ## x: abcissae of data
  ## y: count response
  ## offset: an a priori known component (optional)
  ## w: weights
  ## overdispersion: logical on the presence of
  ##                  overdispersion
  ## ndx: number of internal knots -1.
  ##      Default: floor(length(x)/5)
  ## deg: degree of the B-splines. Default: 3
  ## pord: order of differences. Default: 2
  ## lambda: smoothing parameter. Default: NULL (optional)
  ## df: a number which specifies the degrees of freedom.
  ##     Default: NULL (optional)
  ## method: the method for controlling the amount of
  ##         smoothing. Default: 4
  ## coefstart: eventual initial coefficients
  ## control: a list of control parameters
  
  ## Output: a list containing CHECKED arguments
  ##         for the Mort1Dsmooth function
  m <- length(y)
  
  offsetINIT <- offset
  ## about infinitive or NA offset
  whioff <- which(is.infinite(offset) | is.na(offset))
  whiwei <- which(w==0)
  if(any(!whioff%in%whiwei)){
    stop("weights different from zero associated with infinitive or NA offset values")
  }
  offset[c(whioff, whiwei)] <- 100
  ## about lengths and wrong values
  if (length(x)!=length(y)) 
    stop("Arguments must have same length")
  if (length(y) != m | length(offset) != m) 
    stop("Argument arrays of wrong length")
  if (deg < 1 | deg >= 10) 
    stop("Wrong value for deg")
  if (pord <= 0 | pord >= 5) 
    stop("Wrong value for pord")
  if (ndx < 2 | ndx >= floor(m*.9)) 
    stop("Wrong value for ndx")
  coefstart.check <- is.null(coefstart)
  if(!coefstart.check){
    if(length(coefstart)!=(ndx+deg)){
        stop("coefstart must have length equal to ndx+deg")
      }
  }
  ## about method
  if (method != 1 & method != 2 &
      method != 3 & method != 4) 
    stop("Wrong value for method")
  ## method = 1 adjusts lambda so that the BIC is minimized
  ## method = 2 adjusts lambda so that the AIC is minimized
  ## method = 3 uses the value supplied for lambda
  ## method = 4 adjusts lambda so that the degrees of
  ##          freedom is equal to df
  ## check-point methods
  lambda.check <- is.null(lambda)
  df.check <- is.null(df)
  MET <- NULL
  ## both lambda and df NULL
  if(lambda.check & df.check & method==1){MET=1}
  if(lambda.check & df.check & method==2){MET=2}
  if(lambda.check & df.check & method==3){
    stop("with method 3, provide lambda")
  }
  if(lambda.check & df.check & method==4){
    stop("with method 4, provide df")
  }
  ## lambda NULL and df GIVEN
  if(lambda.check & !df.check & method==1){
    stop("df and method 1 cannot be chosen together")
  }
  if(lambda.check & !df.check & method==2){
    stop("df and method 2 cannot be chosen together")
  }
  if(lambda.check & !df.check & method==3){
    stop("df and method 3 cannot be chosen together")
  }
  if(lambda.check & !df.check & method==4){MET=4}
  ## lambda GIVEN and df NULL
  if(!lambda.check & df.check & method==1){
    stop("lambda and method 1 cannot be chosen together")
  }
  if(!lambda.check & df.check & method==2){
    stop("lambda and method 2 cannot be chosen together")
  }
  if(!lambda.check & df.check & method==3){MET=3}
  if(!lambda.check & df.check & method==4){
    stop("lambda and method 4 cannot be chosen together")
  }
  ## both lambda and df GIVEN, regardless method
  if(!lambda.check & !df.check){
    stop("lambda and df cannot be chosen together")
  }  
  ## impossible values for lambda and df
  if(!lambda.check && lambda<0)
    stop("lambda must be positive")
  if(!df.check && df<pord)
    stop("df must be larger than pord")
  if(!df.check && df>c(ndx+deg))
    stop("df must be smaller than ndx+deg")
  if (!df.check & length(df)!=1)
    stop("df must be length 1")
  if (!lambda.check & length(lambda)!=1)
    stop("lambda must be length 1")
  ## setting control-parameters
  con <- list(MON=FALSE, TOL1=1e-06, TOL2=0.5,
              RANGE=c(10^-4, 10^6), MAX.IT=50)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in control: ",
            paste(noNms, collapse = ", "))
  ## stop about weights
  if(length(w)!=m){ 
    stop("length of w and y must be equal")
  }
  ## warning about interpolation/extrapolation
  if(any(w==0)){
    warning("Interpolation and/or extrapolation is taking place", call. = FALSE)
  }
  ## about the overdispersion parameter
  if(overdispersion & method==3)
    warning("given method 3, overdispersion is computed a posteriori")
  if(overdispersion & method==4)
    warning("given method 4, overdispersion is computed a posteriori")
  ## warning about weights
  if(min(w) < 0) {
    warning(paste("At least one weight entry is negative"))
  }
  ## returning
  llist <- list(x=x, y=y, offset=offset, w=w,
                offsetINIT=offsetINIT,
                overdispersion=overdispersion, m=m,
                ndx=ndx, deg=deg, pord=pord, 
                lambda=lambda, df=df, method=method,
                coefstart=coefstart,
                control=con)
  llist
}
