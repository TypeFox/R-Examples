HHI <- function(shares,
    owner=diag(length(shares)),
    control){

    nprod <- length(shares)

    if(!(is.vector(shares))){
        stop("'shares' must be a vector")}

    if(any(shares < 0,na.rm=TRUE) ||
       any(shares >1,na.rm=TRUE) ) {
        stop("'shares'  must be between 0 and 1")}


    ## transform pre-merger ownership vector into matrix##
    if(is.vector(owner) ||
       is.factor(owner)){

        if(nprod != length(owner)){
            stop("'shares' and 'owner' vectors must be the same length")}

        owners <- as.numeric(factor(owner))
        owner <- matrix(0,ncol=nprod,nrow=nprod)

        for( o in unique(owners)){
            owner[owners == o, owners == o] = 1
        }

        rm(owners)

    }

    else if(!is.matrix(owner) ||
            ncol(owner) != nrow(owner) ||
            ncol(owner) != nprod ||
            any(owner < 0 | owner > 1,na.rm=TRUE)){

               stop("'owner' must be a square matrix whose dimensions equal the length of 'shares' and whose elements are between 0 and 1")
    }


    if(missing(control)){
        control <- owner>0
        }

    else if(is.vector(control) ||
       is.factor(control)){

        if(nprod != length(control)){
            stop("'shares' and 'control' vectors must be the same length")}

        controls <- as.numeric(factor(control))
        control <- matrix(0,ncol=nprod,nrow=nprod)

        for( c in unique(controls)){
            control[controls == c, controls == c] = 1
        }

        rm(controls)
    }

    else if (!is.matrix(control) ||
             ncol(control) != nrow(control) ||
             ncol(control) != nprod ||
             any(control < 0 | control > 1,na.rm=TRUE)
             ){
        stop("'control' must be a square matrix whose dimensions equal the length of 'shares' and whose elements are between 0 and 1")

        }

    weights <- crossprod(control,owner)
    weights <- t(t(weights)/diag(weights)) # divide each element by its corresponding diagonal

    shares <- shares*100
    result <- as.vector(shares %*% weights %*% shares)



    return(result)




    }
