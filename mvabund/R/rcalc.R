################################################################################
# RCALC Sample correlation matrix with 0 variance correction and 2 methods to
# ensure non-singularity
#
#    Finds the correlation matrix R and covariance matrix S for DATA. If any
#    columns of DATA have zero variance, then the nondiagonal elements of R for
#    the corresponding rows and columns are set to 0.
################################################################################

rcalc <- function( dat, cor.type="I", param=0 )  {

    dat       <- as.matrix(dat)
    nVars     <- ncol( dat )
    nRows     <- nrow( dat )
    eyeP      <- matrix(0, nrow=nVars,ncol=nVars) 
    eyeP[1+0:(nVars-1)*(nVars+1)] <- 1

    if (nRows ==1) stop ("More than one observation is necessary.")
    c         <- var( dat )
    # Obtain the diagonal of c.
    d         <- c(c)[1+0:(nVars-1)*(nVars+1)]
    isVar0    <- ( d==0 )
    # Correction for zero variance (so that log-lik. is zero for this variable)
    d[isVar0] <- ( nRows-1 ) /( nRows^2 ) 	
    if ( cor.type=="R"  & (nRows <= nVars) ){
        stop("An unstructured correlation matrix should only be used if N>>p.")
    } else if ( cor.type=="augvar" )  {
        if (nRows <= nVars  & param==0)  {
          stop("An unstructured correlation matrix should only be used if N>>p.")
        }
        param <- param / (nRows-1)
        c     <- c + param*eyeP
        d     <- d + param
    } 

    R                <- c / sqrt( d%*%t(d) )
    # Substitutes NA values that occured because of division with 0:
    R[isVar0,isVar0] <- 0
    R[isVar0,isVar0][1+0:(sum(isVar0)-1)*(sum(isVar0)+1)] <- 1
    
    if ( cor.type=="I" )  {
        R <- eyeP
    } else if ( cor.type=="shrink" ) {
    # Reduce corr matrix by a factor param.
        R <- param*R + (1 - param)*eyeP
   
    }  else  if (! (cor.type=="augvar" | cor.type=="R") )  {
        stop("rcalc: Invalid 'cor.type'!")
    }


    if(length(d)>1) {
        sd <- matrix(0, nVars,nVars)
        sd[1+0:(nVars-1)*(nVars+1)] <- sqrt( d )
    } else sd <- as.matrix(sqrt( d ))
    S  <- sd%*%R%*%sd

    list( R=R, S=S )

}

