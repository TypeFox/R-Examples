LKrig.make.a.wght<- function( a.wght, nlevel, mx, my){
    # if a.wght is not a list with nlevel components then
    # assume this is a scalar or vector of values for center value.
    if (!is.list(a.wght)) {
        # some checks on a.wght
        # coerce a.wght to list if it is passed as something else (most likely a vector)
        if (nlevel == 1) {
            a.wght <- list(a.wght)
        }
        else {
            # repeat a.wght to fill out for all levels.
            if (length(a.wght) == 1) {
                a.wght <- rep(a.wght, nlevel)
            }
            a.wght <- as.list(c(a.wght))
        }
    }
    #################################
    # check length of a.wght list
    if (length(a.wght) != nlevel) {
        stop("length of a.wght list differs than of nlevel")
    }
    #
    # now figure out if the model is stationary
    # i.e. a.wght pattern is to be  repeated for each node
    # this is the usual case
    # if not stationary a.wght should lists of arrays that
    # give values for each node separately
    stationary <- is.null(dim(a.wght[[1]]))
    first.order<- rep( NA, nlevel)
    # simple check on sizes of arrays
    if (stationary) {      
        for (k in 1:length(a.wght)) {
             N.a.wght <- length(a.wght[[1]])
        # allowed lengths for a.wght are just the center 1 values
        # or 9 values for center,  first, and second order neighbors
         if (is.na(match(N.a.wght, c(1, 9)))) {
               stop("a.wght needs to be of length 1 or 9")
        }
           first.order[k]<- ifelse( N.a.wght == 1, TRUE, FALSE)
        }
    
    }
    else {
        for (k in 1:length(a.wght)) {
            dim.a.wght <- dim(a.wght[[k]])
            if( (dim.a.wght[1] != mx[k]) | (dim.a.wght[2] != my[k]) ) {
                stop(paste("a.wght lattice at level", k,
                    " has wrong first two dimensions") )
            }
        first.order[k] <- length( dim.a.wght) == 2
        }
     }
    
    fastNormalization<- stationary & all( first.order )
    if( fastNormalization & any( is.na( unlist( a.wght))) ) {
         fastNormalization<- FALSE}
    
    return( list( fastNormalization= fastNormalization, stationary= stationary, first.order=first.order,
                          a.wght=a.wght))
 }          
