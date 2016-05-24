cmcr.cournot <- function(shares,mktElast){

    if(!is.vector(shares) || length(shares)!=2) { stop("'shares'  must be a vector of length 2")}
    if(!is.vector(mktElast) || length(mktElast)!=1) { stop("'mktElast'  must be a vector of length 1")}
    if(any(shares < 0,na.rm=TRUE) ||
       any(shares >1,na.rm=TRUE) ) {
        stop("'shares'  must be between 0 and 1")}
    if(any(mktElast < 0,na.rm=TRUE)) { stop("'mktElast'  must be positive")}

    mcDelta <- 2*prod(shares)/(mktElast*sum(shares) - sum(shares^2))


    return(mcDelta*100)
}


cmcr.bertrand <- function(prices, margins, diversions, ownerPre, ownerPost=matrix(1,ncol=length(prices), nrow=length(prices)), labels=paste("Prod",1:length(prices),sep=""))
{




    if(!(is.vector(prices) & is.vector(margins))){
        stop("'prices' and 'margins' must be vectors")}

    nprod = length(prices)


    if(nprod != length(margins)){
        stop("'prices'and 'margins' vectors must be the same length")}

    if(any(prices < 0,na.rm=TRUE)){ stop("'prices'  must be non-negative")}
    if(any(margins < 0 || margins > 1,na.rm=TRUE) ){ stop("'margins' vector elements  must be between 0 and 1")}

    if(!is.matrix(diversions)){ stop("'diversions'  must be a matrix")}
    if(!all(diag(diversions) == -1,na.rm=TRUE)){ stop("'diversions' diagonal elements must all equal -1")}
    if(any( abs(diversions) > 1,na.rm=TRUE)){ stop("'diversions' elements must be between -1 and 1")}
    if(ncol(diversions)!=nrow(diversions) ||
       ncol(diversions)!= nprod){
        stop("'diversions' must be a square matrix whose dimension equals the length of 'prices'")
    }

    if(!is.matrix(ownerPost)){ stop("'ownerPost'  must be a matrix")}
    if(any(ownerPost < 0 || ownerPost > 1,na.rm=TRUE)){ stop("'ownerPost' elements must be between 0 and 1")}
    if(
       !isTRUE(all.equal(colSums(unique(ownerPost)),rep(1,nprod)))){
        stop("The columns of the matrix formed from the unique rows of 'ownerPost' must sum to 1")
    }

     ## transform pre-merger ownership vector into matrix##
    if(is.vector(ownerPre) ||
       is.factor(ownerPre)){

        if(nprod != length(ownerPre)){
            stop("'prices'and 'ownerPre' vectors must be the same length")}

        owners <- as.numeric(factor(ownerPre))
        ownerPre <- matrix(0,ncol=nprod,nrow=nprod)

        for( o in unique(owners)){
            ownerPre[owners == o, owners == o] = 1
        }

        rm(owners)
    }




    else if(!is.matrix(ownerPre) ||
            ncol(ownerPre) != nrow(ownerPre) ||
            any(ownerPre < 0 || ownerPre > 1,na.rm=TRUE) ||
            ncol(ownerPre) != nprod ||
            !isTRUE(all.equal(colSums(unique(ownerPre)),rep(1,nprod)))
            ){
            stop("'ownerPre' must be a square matrix whose dimension equals the length of 'prices' and whose elements are between 0 and 1. Also, the columns of the matrix formed from the unique rows of 'ownerPre' must sum to 1")}




    ## weight diversion ratios by price ratios and ownership matrices ##
    priceRatio = tcrossprod(1/prices, prices)
    Bpre =  -1 * diversions * priceRatio * ownerPre;
    Bpost = -1 * diversions * priceRatio * ownerPost;



    ##calculate post-merger margin ##
    marginPost = as.vector( solve(Bpost) %*%
    ((diag(ownerPost)/diag(ownerPre)) * as.vector(Bpre %*% margins))
    )

    ## calculate proportional change in marginal cost ##
    mcDelta= (marginPost - margins)/(1 - margins)

    names(mcDelta) <- labels

    return(mcDelta*100)
}
