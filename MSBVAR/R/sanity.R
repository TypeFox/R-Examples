"sanity.check.prior" <- function(args){

    if(args$lambda0<=0 || args$lambda0>1){
        stop("\n\n\t-- Overall prior tightness (lambda0) must be between 0 and 1\n")
    }
    else if(args$lambda1<0){
        stop("\n\n\t--
        Prior tightness around AR(1) parameters (lambda1) must be greater than 0\n")
    }
    else if(args$lambda3<0){
        stop("\n\n\t-- Order of lag decay (lambda3) must be greater than 0\n")
    }
    else if(args$lambda4<0){
        stop("\n\n\t-- Prior tightness around intercept (lambda4) must be greater than 0\n")
    }
    else if(args$lambda5<0){
        stop("\n\n\t--
        Prior tightness around exogenous parameters (lambda5) must be greater than 0\n")
    }
    else if(args$mu5<0){
        stop("\n\n\t-- Prior weight of sum of coefficients must be non-negative\n")
    }
    else if(args$mu6<0){
        stop("\n\n\t-- Prior values of initial dummy observations or of cointegrating prior must be non-negative\n")
    }
    else{ return(1) }

}


"sanity.check.var" <- function(args){

    arg.names <- names(args)

    if(length(which(arg.names=="dat"))){
        Y <- args$dat
    } else if (length(which(arg.names=="Y"))){
        Y <- args$Y
    }

    # Check time series properties of input Y
    if(is.null(tsp(Y))==TRUE){
        stop("\n\n\t-- Input data [Y] must be of class ts() or mts()\n")
    } else {
        T <- nrow(Y)
        m <- ncol(Y)
    }

##     # Check time series properties of input z
##     if(is.null(args$z)==FALSE || is.ts(args$z)==FALSE){
##         stop("\n\n\t-- Input data [z] must be NULL, a ts() object, or an mts() object\n")
##     }

    # Check if num lags is possible given the amount of data supplied
    if(args$p<=0 || args$p>=T){
        stop("\n\n\t-- Number of lags (p): 0 < p < number of observations (T) \n")
    }

    return(1)

}

"sanity.check.bvar" <- function(args){

    arg.names <- names(args)

    if(length(which(arg.names=="dat"))){
        Y <- args$dat
    } else if (length(which(arg.names=="Y"))){
        Y <- args$Y
    }

    sanity.check.var(list(Y=Y, p=args$p, z=args$z))

    # Make sure prior specifications are at least within the bounds
    # outlined in the documentation

    prior.vals <- list(lambda0=args$lambda0, lambda1=args$lambda1,
                       lambda3=args$lambda3, lambda4=args$lambda4,
                       lambda5=args$lambda5, mu5=args$mu5, mu6=args$mu6)
    sanity.check.prior(prior.vals)

    if(args$qm!=4 && args$qm!=12){
        warning("\n\n\t-- qm should be set to 4 (quarterly) or 12 (monthly): qm=4 selected as default\n")
    }

    # BVAR-specific argument checks

    if(length(which(arg.names=="prior")) && (args$prior<0 || args$prior>2)){
        stop("\n\n\t-- Please specify valid prior form: (0) Normal-Wishart; (1) Normal-flat; or (2) Flat-flat\n")
    }

    if(length(which(arg.names=="posterior.fit")) && !(is.logical(args$posterior.fit))){
        stop("\n\n\t-- For BVAR models, posterior.fit argument must be logical (TRUE/FALSE)\n")
    }

    return(1)
}


"sanity.check.bsvar" <- function(args){

    arg.names <- names(args)

    sanity.check.bvar(args)

    m <- ncol(args[[1]])

    # Check identification MATRIX

    if(is.matrix(args$ident)==FALSE) {
        stop("\n\n\t-- Check Identification: ident argument should be of the class matrix, not list or dataframe\n")
    }

    if(length(which(args$ident==1))>(m*(m+1))/2){
        stop("\n\n\t-- Check Identification: no more than m*(m+1)/2 free parameters allowed\n")
    }

    if(length(which(args$ident==1))<m){
        stop("\n\n\t-- Check Identification: too few free parameters\n")
    }

    return(1)
}


"sanity.check.gibbs" <- function(args){

    if(args$N1<1000) warning("\n\n\t-- Burnin iteration count may be too low for sensible results (<1000)\n")
    if(args$N2<1000) warning("\n\n\t-- Gibbs sampler iteration count may be too low for sensible results (<1000)\n")
    if(args$thin<1) stop("\n\n\t-- Thinning parameter must be greater than or equal to 1\n")
    methodlist <- c("DistanceMLA", "DistanceMLAhat", "Euclidean", "PositiveDiagA", "PositiveDiagAinv")

    if(length(which(methodlist==args$normalization))==0){
        warning("\n\n\t-- Specified normalization method does not exist: default set to 'DistanceMLA'\n")
        return("DistanceMLA")
    }
    else { return(1) }

}


"sanity.check.irf" <- function(args){

    if(is.null(args$nsteps)){
        stop("\n\n\t-- 'nsteps' undefined: number of steps for IRF must be user defined\n")
    }
    else{ return(1) }

}


"sanity.check.mc.irf" <- function(args){

    arg.names <- names(args)

    if(args$nsteps<=0) stop("\n\n\t-- 'nsteps' must be greater than 0\n")

    if(attr(args$varobj,"class")=="VAR" && !(args$nsteps>0)){
        stop("\n\n\t--
        For VAR models, argument 'nsteps' must be defined and non-negative integer\n")
    }

    if(attr(args$varobj,"class")=="BVAR" && !(args$draws>0)){
        stop("\n\n\t--
        For BVAR models, argument 'draws' must be defined and non-negative integer\n")
    }

    if(attr(args$varobj,"class")=="BSVAR" && is.null(args$A0.posterior)){
        stop("\n\n\t-- mc.irf.BSVAR requires A0.posterior object from gibbs.A0()\n")
    }

    return(1)
}
