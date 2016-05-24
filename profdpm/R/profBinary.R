profBinary <- function(formula, data, clust, param, method="agglomerative", 
                       maxiter=1000, crit=1e-6, verbose=FALSE, sampler=FALSE) {
    ###################################################
    if(missing(data))
        data <- sys.frame(sys.parent(0))
    mf <- model.frame(formula, data)
    mm <- model.matrix(formula, data)

    if(missing(clust)) {
        clust <- FALSE 
    } else { 
        clust <- as.factor(clust)
        if(nrow(mm) != length(clust))
          stop("nrow(data) must equal length(clust)")  
    }
    if(missing(param)) {
        param <- list() 
    } else if(!is.list(param)) {
        warning("param must be a list, using defaults")
        param <- list() 
    } else {
        if(length(names(param)) == 0) {
            warning("param argument does not include any named items, using defaults")
            param <- list()
        } else if(length(param) > length(names(param))) {
            warning("param contains unnamed items")
        }
    }

    if(is.null(param$alpha)) param$alpha <- 1/1000
    if(is.null(param$a0)) param$a0 <- 1.00
    if(is.null(param$b0)) param$b0 <- 1.00

    if(!is.character(method)) {
        warning("method must be a character string, using default")
        method <- "stochastic"
    }
    if(!is.numeric(maxiter) | maxiter < 0) {
        warning("maxiter must be numeric and non-negative, using default")
        maxiter <- 1000
    }
    if(!is.numeric(crit) | crit < 0) { 
        warning("crit must be numeric and non-negative, using default") 
        crit <- 1e-5
    }
    if(!is.logical(verbose)) {
        warning("verbose must be a logical, using default")
        verbose <- FALSE
    }
    if(!is.logical(sampler)) {
        warning("sampler must be a logical, using default")
        sampler <- FALSE 
    }
    if(sampler && method != "gibbs")
        warning("'sampler' has no effect for methods other than 'gibbs'")

    ###################################################
    #remove missing observations, issue warning
    miss <- apply( is.na( mm ), 1, any ) 
    rmm <- mm[!miss,]
    rmf <- mf[!miss,]
    if( is.logical(clust) ) { rc <- FALSE }
    else { rc <- as.integer(unclass(clust[!miss])-1) }
    if( any( miss ) ) {
        warning( "removed observations with missing values: ", 
        paste(" ", which(miss), sep="") )
    }

    ###################################################
    #convert method to integer
    if(      method == "none" )          { method <- 0 }
    else if( method == "stochastic" )    { method <- 1 }
    else if( method == "agglomerative" ) { method <- 2 }
    else if( method == "gibbs" )         { method <- 3 }
    else if( method == "fast" )          { method <- 4 }
    else {
        method <- 2 #default is "agglomerative"
        warning("method must be \'stochastic\', \'agglomerative\',
            \'gibbs\', \'fast\' or \'none\'", )
    }

    ###################################################
    #convert mm to integer storage, 
    #check that all are binary
    storage.mode(rmm) <- "integer"
    if( any( (rmm != 1L) & (rmm != 0L) ) )
        stop("formula must specify only 0s and 1s")

    ###################################################
    #call the C function
    ret <- .Call("profBinary", rmm, rc, as.list(param), as.integer(method),
                  as.integer(maxiter), as.double(crit), verbose, sampler, 
                  PACKAGE="profdpm")  
    ret$model <- rmf
    return(ret)
}
