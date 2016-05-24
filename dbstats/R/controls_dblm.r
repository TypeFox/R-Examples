
 ########################
 #### controls_dblm #####
 ########################

 ## description: Internal function. Check if the attributes are consistents.
 ##
 ##        Inputs:  G, weights, eff.rank, rel.gvar, method and y.
 ##        Outputs: eff.rank
 ##


controls_dblm <- function(G,weights,eff.rank,rel.gvar,method,y){
                    
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

   # program controls: eff.rank
   if (missing(eff.rank)&&(method=="eff.rank"))
    stop("eff.rank is not defined. You must assign a number to eff.rank")
   else if (missing(eff.rank)||method!="eff.rank")
    eff.rank<-0
   if (!missing(eff.rank) && !is.numeric(eff.rank))
    stop("'eff.rank' must be numeric")
   if (!(missing(eff.rank))&&(method=="eff.rank")){
    if (eff.rank<0 || eff.rank>(n-1))
     stop("'eff.rank' must be in the interval (0,n-1]")
    else
     ini_eff.rank<-eff.rank
   }
   
   # program controls: rel.gvar
   if (method=="rel.gvar"&&is.numeric(rel.gvar)&&((rel.gvar<0)||(rel.gvar>1)))
      stop("'rel.gvar must be between 0 and 1")
   if (method=="rel.gvar"&&!is.numeric(rel.gvar))
      stop("'rel.gvar' must be a float number between 0 and 1")

   # control that the G dimension is the same tha y dimension.
   if (nrow(as.matrix(y))!=nrow(as.matrix(G))) {
    stop(gettextf("The dimensions of the G matrix: nrow=ncol= %d, should equal %d (number of observations)",
                nrow(as.matrix(G)),nrow(as.matrix(y))))
   }
   return(list(eff.rank=eff.rank,weights=weights))
 }