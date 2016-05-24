
 ########################
 #### controls_dbpls ####
 ########################

 ## description: Internal function. Check if the attributes are consistents.
 ##
 ##        Inputs:  distance, weights, eff.rank, y.
 ##        Outputs: weights
 ##


controls_dbplsr <- function(distance,weights,ncomp,y){

   # program controls: weights
   if (missing(weights)||is.null(weights))
    weights<-rep(1,nrow(as.matrix(y)))
   if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
   if (sum(weights<0)>0)
    stop("Weights array weights is not valid: some weights are negative")
   if (sum(weights)==0)
    stop("Weights array weights is not valid: sum(weights)=0")

   # program controls: response y
   if (missing(y)||is.null(y))
    stop("the response variable must be defined")
   if (!is.numeric(y)&&!is.data.frame(y)&&!is.matrix(y))
    stop("the response 'y' must be numeric, data.frame or matrix")
   n <- nrow(as.matrix(y)) # number of observations

   # program controls: ncomp
   if (missing(ncomp))
    stop("'ncomp' is not defined. You must assign a number to ncomp")
   if (!missing(ncomp) && !is.numeric(ncomp))
    stop("'ncomp' must be numeric")
   if (!missing(ncomp)&&(ncomp<0||ncomp>(n-1)))
    stop("'ncomp' must be between 1 and number of observations - 1.")

   # control that the distance dimension is the same tha y dimension.
   if (nrow(as.matrix(y))!=nrow(as.matrix(distance))) {
    stop(gettextf("The dimensions of the distance matrix: nrow=ncol= %d, should equal %d (number of observations)",
                nrow(as.matrix(distance)),nrow(as.matrix(y))))
   }
   return(weights)
 }