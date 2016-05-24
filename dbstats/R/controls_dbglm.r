
 ########################
 #### controls_dbglm #####
 ########################

 ## description: Internal function. Check if the attributes are consistents.
 ##
 ##        Inputs:  G, weights, eff.rank, rel.gvar, method and y.
 ##        Outputs: eff.rank
 ##


controls_dbglm <- function(distance,weights,offset,rel.gvar,maxiter,eps1,
                            eps2,y, method)
{
  # program controls: y
   if (missing(y)||is.null(y))
    stop("the response variable must be defined")
   if (!is.numeric(y)&&!is.data.frame(y)&&!is.matrix(y)&&!is.factor(y)&&!is.logical(y))
    stop("the response 'y' must be numeric, data.frame, matrix or a factor")

  nobs <- nrow(as.matrix(y))
  # program controls: weights
   if (missing(weights)||is.null(weights))
     weights<-rep(1,nobs)

   if (!is.null(weights)) {
      if (!is.numeric(weights))
        stop(" 'weights' must be a numeric vector")
      if (any(weights < 0))
        stop(" negative weights not allowed")
      if (sum(weights)==0)
        stop(" Weights array weights is not valid: sum(weights)=0")
      if (length(weights) != nobs)
         stop(gettextf("number of weights is %d should equal %d (number of observations)",
              length(weights),nobs), domain = NA)
   }

   # program controls: offset
   if (missing(offset)||is.null(offset))
     offset<-rep(0,nobs)

   if (!is.null(offset)) {
     if (length(offset) != nobs)
         stop(gettextf("number of offsets is %d should equal %d (number of observations)",
              length(offset),nobs), domain = NA)
   }
   
   # program controls: rel.gvar
   if (method=="rel.gvar"&&is.numeric(rel.gvar)&&((rel.gvar<0)||(rel.gvar>1)))
      stop("'rel.gvar must be between 0 and 1")
   if (method=="rel.gvar"&&!is.numeric(rel.gvar))
      stop("'rel.gvar' must be a float number between 0 and 1")
      
   # program controls: maxiter
   if (is.null(maxiter))
      maxiter <- 100
   if (length(maxiter)>1)
      stop("'maxiter' must be of length one")
   if (maxiter<=0)
      stop("'maxiter' must be > 0")

    # program controls: epsilons
   if (is.null(eps1))
    eps1=1e-10
   if (is.null(eps2))
    eps2=1e-10
   if (eps1>=1||eps1<=0)
    stop("'eps1' must be between the interval (0,1)")
   if (eps2>=1||eps2<=0)
    stop("'eps2' must be between the interval (0,1)")

    return(list(weights=weights,offset=offset,rel.gvar=rel.gvar,maxiter=maxiter,
                  eps1=eps1,eps2=eps2))

}
