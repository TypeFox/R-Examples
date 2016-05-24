dreFitCond <- function(object, omodel = TRUE, rootFinder = findRoots, ...){

    ## Fit exposure nuisance model
    e.fit <- geeFitCond(y = object$a,
                        x = object$z,
                        link = object$elink,
                        id = object$id,
                        rootFinder = rootFinder,
                        ...)
    alpha.hat <- e.fit$coefficients
    ## names(alpha.hat) <- colnames(object$z)
    x.res.e <- apply(object$x, 2, '*', e.fit$res)

    if (omodel) {

        ## ## Center variables
        ## v.cent <- object$v - apply(object$v,2,function(z) ave(z,object$id))
        ## Fit outcome nuisance model
        o.fit <- geeFitCond(y = object$y,
                            x = cbind(object$ax,object$v),
                            link = object$olink,
                            id = object$id,
                            rootFinder = rootFinder, ...)
        beta1.hat <- o.fit$coefficients[1:ncol(object$ax)]
        ## names(beta1.hat) <- colnames(object$ax)
        gamma.hat <- o.fit$coefficients[-(1:ncol(object$ax))]
        ## names(gamma.hat) <- colnames(object$v)

        v.cent <- .Call("center", object$v, object$id, PACKAGE = "drgee")

    }

    y.cent <- .Call("center", object$y, object$id, PACKAGE = "drgee")
    ax.cent <- .Call("center", object$ax, object$id, PACKAGE = "drgee")

    if (object$olink == "identity") {

        ## y.cent <- object$y - ave(object$y, object$id)
        
        if (omodel) {
            y.star.cent <- y.cent -  v.cent %*% gamma.hat
        } else {
            y.star.cent <- y.cent
        }


        ## lhs %*% beta = rhs
        lhs <- crossprod( x.res.e, ax.cent)
        rhs <- crossprod( x.res.e, y.star.cent)
        beta.hat <- as.vector( solve(lhs) %*% rhs )
        names(beta.hat) <- colnames(object$ax)

        ## Estimating equations and derivatives of the doubly robust estimating equations
        ## The rows are the columns in u1
        ## and the columns are the partial derivatives

        s.o.cent <- as.vector( y.star.cent - ax.cent %*% beta.hat )
        x.res.o <- apply(object$x, 2, '*', s.o.cent)

        u1 <- apply( x.res.e, 2, '*', s.o.cent)

        d.u1.beta <- crossprod( x.res.e, -ax.cent )

        d.u1.alpha <- crossprod(x.res.o, e.fit$d.res)

        if (omodel) {
            d.u1.beta1 <- matrix( rep(0, ncol(object$ax)^2), ncol=ncol(object$ax) )
            d.u1.gamma <- crossprod( x.res.e, -v.cent )
        }

        optim.object <- NULL

    } else if (object$olink == "log") {

        if (omodel) {
            y.star <- object$y * exp(-object$v %*% gamma.hat)
        } else {
            y.star <- object$y
            beta1.hat <- coef(geeFitCond(y = object$y,
                                         x = object$ax,
                                         link = object$olink,
                                         id = object$id,
                                         rootFinder = rootFinder,
                                         ...))
        }

        u.func <- function(beta, arg.list) {
            s.o <- arg.list$y.star * exp(-arg.list$ax %*% beta)
            s.o.cent <- .Call("center", as.matrix(s.o), arg.list$id, PACKAGE = "drgee")
            ## s.o.cent <- s.o - ave(s.o, arg.list$id)
            return( as.vector( crossprod( arg.list$x.res.e, s.o.cent) ) )
        }

        all.args <- c(list(beta.init = beta1.hat, eq.func = u.func, d.eq.func = NULL,
                           arg.list = list(y.star = y.star, ax = object$ax,
                               x.res.e = x.res.e, id = object$id)
                           ),
                      list(...)
                      )

        ## Call equation solver with beta.init as initial guess and eq.func as estimation function
        root.object <- do.call(rootFinder, all.args)
        beta.hat <- root.object$roots
        optim.object <- root.object$optim.object

        ## Estimating equations and derivatives of the doubly robust estimating equations
        ## The rows of d.u1 are the columns in u.dr
        ## and the columns are the partial derivatives
        s.o <- as.vector( y.star * exp(- object$ax %*% beta.hat) )
        s.o.cent <- .Call("center", as.matrix(s.o), object$id, PACKAGE = "drgee")
        ## s.o.cent <- s.o - ave(s.o, object$id)
        x.s.o.cent <- apply(object$x, 2, '*', s.o.cent)

        u1 <- apply(x.res.e, 2, '*', s.o.cent)

        d.s.o.beta <- apply(-object$ax, 2, '*', s.o)
        ## d.s.o.beta.cent <- d.s.o.beta - apply(d.s.o.beta, 2, function(z) ave(z,
        ## object$id))
        d.s.o.beta.cent <- .Call("center", as.matrix(d.s.o.beta), object$id, PACKAGE = "drgee")
        d.u1.beta <- crossprod( x.res.e, d.s.o.beta.cent  )

        d.u1.alpha <- crossprod(x.s.o.cent, e.fit$d.res)

        d.u1 <- cbind(d.u1.beta, d.u1.alpha)

        if (omodel) {
            d.u1.beta1 <- matrix( rep(0, ncol(object$ax)^2), ncol = ncol(object$ax))

            d.s.o.gamma <- apply(-object$v, 2, '*', s.o)
            ## d.s.o.gamma.cent <- d.s.o.gamma - apply(d.s.o.gamma, 2, function(z)
            ## ave(z,object$id))
            d.s.o.gamma.cent <- .Call("center", as.matrix(d.s.o.gamma), object$id, PACKAGE = "drgee")
            d.u1.gamma <- crossprod( x.res.e, d.s.o.gamma.cent  )
        }


    } else {
        
        stop("\nOutcome link needs to be identity or log")
        
    }

    ## Estimating equations and derivatives of outcome and exposure estimating equations
    ## The rows are the columns in u
    ## and the columns are the partial derivatives

    z.cent <- .Call("center", object$z, object$id, PACKAGE = "drgee")
    ## Exposure nuisance model estimating equations
    u2 <- apply(z.cent, 2, '*', e.fit$res)

    d.u2.beta <-  matrix( rep(0, ncol(object$z) * ncol(object$ax) ),
                         nrow = ncol(object$z) )

    d.u2.alpha <- crossprod(z.cent, e.fit$d.res)

    ## Outcome nuisance model estimating equations
    if (omodel) {

        u3 <- apply(cbind(ax.cent,v.cent), 2, '*', o.fit$res)

        u <- cbind(u1, u2, u3)

        d.u2.beta1.gamma <-  matrix( rep(0, ncol(object$z) * (ncol(object$ax) +
                                                              ncol(object$v))),
                                    nrow = ncol(object$z))

        d.u3.beta.alpha <- matrix( rep(0, (ncol(object$ax) + ncol(object$v)) *
                                       (ncol(object$ax) + ncol(object$z))), ncol
                                  = ncol(object$ax) + ncol(object$z))

        d.u3.beta1.gamma <- crossprod( cbind(ax.cent, v.cent), o.fit$d.res )

        d.u <- rbind(cbind(d.u1.beta, d.u1.alpha, d.u1.beta1, d.u1.gamma),
                     cbind(d.u2.beta, d.u2.alpha, d.u2.beta1.gamma),
                     cbind(d.u3.beta.alpha, d.u3.beta1.gamma) ) / nrow(u)

        coefficients <- c(beta.hat, alpha.hat, beta1.hat, gamma.hat)
        coef.names <- c(object$ax.names, object$z.names, object$ax.names, object$v.names)
        
    } else {

        u <- cbind(u1, u2)

        d.u <- rbind(cbind(d.u1.beta, d.u1.alpha),
                     cbind(d.u2.beta, d.u2.alpha)) / nrow(u)

        coefficients <- c(beta.hat, alpha.hat)
        coef.names <- c(object$ax.names, object$z.names)
    }

    names(coefficients) <- coef.names
    
    ## Calculate variance of all estimates

    vcov <- robVcov(u, d.u, object$id)

    dimnames(vcov) <- list(coef.names, coef.names)

    result <- list(coefficients = coefficients, vcov = vcov,
                   optim.object = optim.object)

    return(result)

}
