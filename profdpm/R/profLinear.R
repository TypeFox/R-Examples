profLinear <- function(formula, data, group, clust, param, method="agglomerative", 
                       maxiter=1000, crit=1e-6, verbose=FALSE, sampler=FALSE) {
    ###################################################
    #argument checking
    if(missing(data))
        data <- sys.frame(sys.parent(0))
    mf <- model.frame(formula, data, drop.unused.levels=TRUE)
    mm <- model.matrix(formula, mf)
    mr <- model.response(mf, "numeric")
    if(missing(group)) {
        data$group <- seq(1, length(mr))
    } else {
        data$group <- as.factor(eval(substitute(group), data, environment(formula)))
    }
    mf <- model.frame(formula, data, group=group, drop.unused.levels=TRUE)

    if(missing(clust)) { 
        clust <- FALSE 
    } else {
        data$clust <- as.factor(clust)
        mf <- model.frame(formula, data, group=group, clust=clust,
            drop.unused.levels=TRUE)
        mg <- model.extract(mf, "group")
        mc <- model.extract(mf, "clust")
        for(g in unique(mg))
            if(length(unique(mc[mg==g])) > 1)
                stop("clust and group are conflicting") 
    }

    if(missing(param)) { 
        param <- list()
    } else if(!is.list(param)) {
        warning("param must be a list, using defaults")
        param <- list()
    } else {
        if(length(names(param)) == 0) {
            warning("param argument has no named items, using defaults")
            param <- list()
        } else if(length(param) > length(names(param))) {
            warning("param has unnamed items")
        }
    }

    if(is.null(param$alpha)) param$alpha <- 1/1000
    if(is.null(param$a0)) param$a0 <- 0.001
    if(is.null(param$b0)) param$b0 <- 0.001
    if(is.null(param$m0)) param$m0 <- rep(0,ncol(mm))
    if(is.null(param$s0)) param$s0 <- 0.001 

    if(!is.character(method))
        stop("method must be a character string")
    if(!is.numeric(maxiter) | maxiter < 0)
        stop("maxiter must be numeric and non-negative")
    if(!is.numeric(crit) | crit < 0)
        stop("crit must be numeric and non-negative") 
    if(!is.logical(verbose))
        stop("verbose must be a logical")
    if(!is.logical(sampler))
        stop("sampler must be a logical")
    if(sampler && method != "gibbs")
        warning("'sampler' has no effect for methods other than 'gibbs'")

    ###################################################
    #order the data according to group
    #convert ordered group to integers from 0,1,...
    #convert ordered clust to integers from 0,1,...
    ord <- order(mf[["(group)"]])
    mfo <- mf
    mf  <- mf[ord,]
    mr  <- model.response(mf, "double")
    mm  <- model.matrix(formula, mf)
    mg  <- model.extract(mf, "group")
    mg  <- as.integer(unclass(mg)-1)
    if(!is.logical(clust)) {
        mc  <- model.extract(mf, "clust")
        mc  <- as.integer(unclass(mc)-1)
    } else {
        mc  <- FALSE
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
    #call the C function
    ret <- .Call("profLinear", mr, mm, mg, mc, as.list(param), as.integer(method),
                  as.integer(maxiter), as.double(crit), verbose, sampler, PACKAGE="profdpm")

    ###################################################
    #undo ordering
    ret$y[ord] <- ret$y
    ret$x[ord,] <- ret$x
    ret$group[ord] <- ret$group
    ret$clust[ord] <- ret$clust
    #ret$clust <- unclass(as.factor(ret$clust))
    #attributes(ret$clust) <- NULL
    ret$model <- mfo
    return(ret)  
}
