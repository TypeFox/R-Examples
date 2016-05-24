DiffSim <- function(n,start,Delta,driftpar,Sigma,seed,thin=1,Model){
    start <- matrix(start,ncol=2)
    driftpar <- matrix(driftpar,nrow=1)
   if(!(Model %in% c("OUa","OU","CIR","FHN","FHN5"))){
     stop("Model must be either 'OU', 'CIR', 'FHN' or 'FHN5'.")
   }
    ## if(is.null(seed)){
    ##   seed <- sample(x=1:100000,size=1)
    ## }
    A <- .Call( "DiffSim",
                 n,
                 start,
                 Delta,
                 driftpar,
                 Sigma,
                 seed,
                 thin,
                 Model,
#                 DUP=FALSE,
                 PACKAGE = "BIPOD" )
      return(A)
}

