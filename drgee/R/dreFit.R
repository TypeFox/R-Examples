dr.logit.eq.func <- function(beta, arg.list) {

    exp.beta.x <- exp( arg.list$x.f %*% beta )
    exp.beta.ax <- exp( arg.list$ax.f %*% beta )
    
    p <- exp.beta.x * arg.list$exp.alpha.z * ( 1 + arg.list$exp.gamma.v )
    q <- p + 1 + exp.beta.x * arg.list$exp.gamma.v
    e.star <- as.vector( p / q )

    res.a.e.star <- as.vector(arg.list$a.f - e.star)
    res.y <- as.vector( arg.list$y.f - 1 + 1/(1 + exp.beta.ax * arg.list$exp.gamma.v) )
    
    return( as.vector( crossprod( arg.list$x.f, res.a.e.star * res.y ) ) )

}


d.dr.logit.eq.func <- function(beta, arg.list) {

    exp.beta.x <- as.vector(exp( arg.list$x.f %*% beta ))
    exp.beta.ax <- as.vector(exp( arg.list$ax.f %*% beta ))

    p <- exp.beta.x * arg.list$exp.alpha.z * (1 + arg.list$exp.gamma.v)
    q <- p + 1 + exp.beta.x * arg.list$exp.gamma.v
    e.star <- p / q

    res.a.e.star <- as.vector(arg.list$a.f - e.star)
    res.y <- as.vector( arg.list$y.f - 1 +
                        1/(1 + exp.beta.ax * arg.list$exp.gamma.v) )

    d.res.a.e.star <- -arg.list$x.f * (e.star / q) * ( q - p - exp.beta.x * arg.list$exp.gamma.v)
    d.res.y <- -arg.list$ax.f * exp.beta.ax * arg.list$exp.gamma.v / (1 + exp.beta.ax * arg.list$exp.gamma.v)^2
    
    return( crossprod( arg.list$x.f , d.res.a.e.star * res.y +
                                      d.res.y * res.a.e.star ) )
}

dreFit <-
    function(object, omodel = TRUE, rootFinder = findRoots, ...) {

        if (object$olink == "logit") {

            o.fit <- geeFit(object$y, cbind(object$ax, object$v), "logit")

            e.fit <- geeFit(object$a, cbind(object$yx, object$z), "logit")

            beta1.hat <- o.fit$coefficients[object$ax.names]
            gamma.hat <- o.fit$coefficients[object$v.names]
            beta2.hat <- e.fit$coefficients[object$yx.names]
            alpha.hat <- e.fit$coefficients[object$z.names]

            exp.gamma.v <- as.vector(exp( object$v %*% gamma.hat ))
            exp.alpha.z <- as.vector(exp( object$z %*% alpha.hat ))

            ## all.args <- c(list(beta.init = beta1.hat,
            ##                    eq.func = dr.logit.eq.func,
            ##                    d.eq.func = d.dr.logit.eq.func,
            ##                    arg.list = list(y.f = object$y,
            ##                                    ax.f = object$ax,
            ##                                    a.f = object$a,
            ##                                    x.f = object$x,
            ##                                    exp.gamma.v = exp.gamma.v,
            ##                                    exp.alpha.z = exp.alpha.z)),
            ##               list(...))

            all.args <- c(list(beta.init = beta1.hat,
                               eq.func = dr.logit.eq.func,
                               d.eq.func = NULL,
                               arg.list = list(y.f = object$y,
                                               ax.f = object$ax,
                                               a.f = object$a,
                                               x.f = object$x,
                                               exp.gamma.v = exp.gamma.v,
                                               exp.alpha.z = exp.alpha.z)),
                          list(...))

            ## Call equation solver with beta.init as initial guess and eq.func as estimation function
            root.object <- do.call(rootFinder, all.args)
            beta.hat <- root.object$roots
            optim.object <- root.object$optim.object

            exp.beta.ax <- as.vector(exp( object$ax %*% beta.hat ) )
            exp.beta.x <- as.vector(exp(  object$x %*% beta.hat ) )
            exp.beta1.ax <- as.vector(exp( object$ax %*% beta1.hat ) )
            exp.beta2.yx <- as.vector(exp( object$yx %*% beta2.hat ) )
            exp.gamma.v <- as.vector(exp( object$v %*% gamma.hat ) )
            exp.alpha.z <- as.vector(exp( object$z %*% alpha.hat ) )

            p <- as.vector(exp.beta.x * exp.alpha.z * (1 + exp.gamma.v) )
            q <- as.vector(p + 1 + exp.beta.x * exp.gamma.v)
            e.star <- as.vector(p / q)

            res.a.e.star <- as.vector(object$a - e.star)
            res.y <- as.vector(object$y - 1 + 1/(1 + exp.beta.ax * exp.gamma.v))

            d.res.a.e.star.gamma <- -object$v * exp.gamma.v * (exp.beta.x / q) *
                (exp.alpha.z - e.star * (exp.alpha.z + 1))
            d.res.y.gamma <- -object$v * exp.beta.ax * exp.gamma.v / (1 +
                                                                      exp.beta.ax * exp.gamma.v)^2

            d.res.a.e.star.alpha <- -object$z * e.star * (1 - e.star)

            U <- cbind(object$x * res.a.e.star * res.y,
                       cbind(object$yx, object$z) * e.fit$res,
                       cbind(object$ax, object$v) * o.fit$res)

            ## Derivative of the doubly robust estimating equations
            d.U1.beta <- d.dr.logit.eq.func(beta.hat, all.args$arg.list)
            d.U1.beta2 <-
                matrix(rep(0, ncol(object$x) * ncol(object$yx)), nrow = ncol(object$x))
            d.U1.alpha <- crossprod( object$x , d.res.a.e.star.alpha * res.y )
            d.U1.beta1 <-
                matrix(rep(0, ncol(object$x) * ncol(object$ax)), nrow = ncol(object$x))
            d.U1.gamma <- crossprod( object$x , d.res.a.e.star.gamma * res.y +
                                                d.res.y.gamma * res.a.e.star )
            d.U1 <- cbind(d.U1.beta, d.U1.beta2, d.U1.alpha, d.U1.beta1, d.U1.gamma)

            ## Derivative of the estimating equations for
            ## the exposure model
            d.U2.beta <- matrix(rep(0, (ncol(object$yx) + ncol(object$z)) *
                                    ncol(object$ax)), ncol = ncol(object$ax))
            d.U2.beta2 <- crossprod( cbind(object$yx, object$z), e.fit$d.res[,colnames(object$yx)])
            d.U2.alpha <- crossprod( cbind(object$yx, object$z), e.fit$d.res[,colnames(object$z)])
            d.U2.beta1 <- matrix(rep(0, (ncol(object$yx) + ncol(object$z)) *
                                     ncol(object$ax)),ncol = ncol(object$ax))
            d.U2.gamma <- matrix(rep(0, (ncol(object$yx) + ncol(object$z)) *
                                     ncol(object$v)), ncol = ncol(object$v))
            d.U2 <- cbind(d.U2.beta, d.U2.beta2, d.U2.alpha, d.U2.beta1, d.U2.gamma)
            
            
            ## Derivative of the estimating equations for
            ## the outcome model
            d.U3.beta <- matrix(rep(0, (ncol(object$ax) + ncol(object$v)) *
                                    ncol(object$ax)), ncol = ncol(object$ax))
            d.U3.beta2 <- matrix(rep(0, (ncol(object$ax) + ncol(object$v)) *
                                     ncol(object$yx)), ncol = ncol(object$yx))
            d.U3.alpha <- matrix(rep(0, (ncol(object$ax) + ncol(object$v)) * ncol(object$z)),ncol=ncol(object$z))
            d.U3.beta1 <- crossprod( cbind(object$ax, object$v), o.fit$d.res[,colnames(object$ax)])
            d.U3.gamma <- crossprod( cbind(object$ax, object$v), o.fit$d.res[,colnames(object$v)])

            d.U3 <- cbind(d.U3.beta, d.U3.beta2, d.U3.alpha, d.U3.beta1, d.U3.gamma)
                
            d.U <- rbind( cbind(d.U1.beta, d.U1.beta2, d.U1.alpha, d.U1.beta1, d.U1.gamma),
                         cbind(d.U2.beta, d.U2.beta2, d.U2.alpha, d.U2.beta1, d.U2.gamma),
                         cbind(d.U3.beta, d.U3.beta2, d.U3.alpha, d.U3.beta1, d.U3.gamma)
                         ) / nrow(U)

            coefficients <- c(beta.hat, beta2.hat, alpha.hat, beta1.hat, gamma.hat)
            coef.names <- c(object$ax.names, object$ax.names, object$z.names, object$ax.names, object$v.names)

            ## If outcome link is identity or log
        } else {

            if(omodel){
                
                o.fit <- geeFit(object$y, cbind(object$ax, object$v), object$olink)
                beta1.hat <- o.fit$coefficients[colnames(object$ax)]
                gamma.hat <- o.fit$coefficients[colnames(object$v)]
                v.gamma.hat <- object$v %*%o.fit$coefficients[colnames(object$v)]

            }else{
                
                v.gamma.hat <- rep(0, nrow(object$y))
                
            }
            
            e.fit <- geeFit(object$a, object$z, object$elink)
            alpha.hat <- e.fit$coefficients

            if(object$olink=="identity"){

                w <- t(object$x * e.fit$res)

                beta.hat <- as.vector(solve(w %*% object$ax) %*% w %*%
                                      (object$y - v.gamma.hat))
                
                res.y <- as.vector(object$y - object$ax %*% beta.hat - v.gamma.hat)

                optim.object <- NULL

                d.res.y.beta <- -object$ax

                if(omodel){
                    d.res.y.gamma <- -object$v
                }

            } else if(object$olink == "log") {

                eq.func <- function(beta, arg.list) {
                    crossprod(arg.list$w, arg.list$y * exp(-arg.list$ax %*% beta) - arg.list$exp.v.gamma)
                }

                d.eq.func <- function(beta, arg.list) {
                    crossprod(arg.list$w * as.vector(arg.list$y *
                                                     exp(-arg.list$ax %*% beta)), -arg.list$ax)
                }

                w <- object$x * e.fit$res

                if(omodel) {
                    exp.v.gamma <- exp(object$v %*% gamma.hat)
                } else {
                    ## Even if we don't have an outcome nuisance model
                    ## we need a start value for beta
                    o.fit <- geeFit(object$y, cbind(1,object$ax), object$olink)
                    beta1.hat <- o.fit$coefficients[-1]
                    exp.v.gamma <- rep(0, nrow(object$y))
                }

                all.args <- c(list(beta.init = beta1.hat, eq.func = eq.func,
                                   d.eq.func = d.eq.func,
                                   arg.list = list(w = w,
                                       exp.v.gamma = exp.v.gamma,
                                       ax = object$ax, y = object$y)),
                              list(...))

                ## call equation solver with beta.init as initial guess and eq.func as estimation function
                root.object <- do.call(rootFinder, all.args)
                beta.hat <- root.object$roots
                optim.object <- root.object$optim.object
                res.y <- as.vector(object$y * exp(-object$ax %*% beta.hat) - exp.v.gamma)
                d.res.y.beta <- -object$ax * as.vector((object$y) * exp(-object$ax %*% beta.hat))

                if(omodel){
                    d.res.y.gamma <- -object$v * as.vector(exp(object$v %*% gamma.hat))
                }
            }

            U1 <- object$x * e.fit$res * res.y

            d.U1.beta <- crossprod(object$x * e.fit$res, d.res.y.beta)
            d.U1.alpha <- crossprod(object$x * res.y, e.fit$d.res)

            U2 <- e.fit$eq.x * e.fit$res
            d.U2.beta <- matrix(rep(0, ncol(object$z) * ncol(object$ax)), nrow = ncol(object$z))
            d.U2.alpha <- crossprod(e.fit$eq.x, e.fit$d.res)

            if(omodel) {

                U3 <- o.fit$eq.x * o.fit$res
                U <- cbind( U1, U2, U3 )

                d.U1.beta1 <- matrix(rep(0, ncol(object$x) * ncol(object$ax)),
                                     nrow = ncol(object$x))
                d.U1.gamma <- crossprod(object$x * e.fit$res, d.res.y.gamma) / nrow(U)

                d.U2.beta1.gamma <-
                    matrix(rep(0, ncol(object$z) * (ncol(object$ax) +
                                                    ncol(object$v))), nrow = ncol(object$z))

                d.U3.beta.alpha <-
                    matrix(rep(0, (ncol(object$ax) + ncol(object$v)) *
                               (ncol(object$ax) + ncol(object$z) ) ), ncol = ncol(object$ax) + ncol(object$z))
                d.U3.beta1.gamma <- crossprod(o.fit$eq.x, o.fit$d.res) / nrow(U)

                d.U <- rbind( cbind(d.U1.beta, d.U1.alpha, d.U1.beta1, d.U1.gamma),
                             cbind(d.U2.beta, d.U2.alpha, d.U2.beta1.gamma),
                             cbind(d.U3.beta.alpha, d.U3.beta1.gamma)) / nrow(U)

                coefficients <- c(beta.hat, alpha.hat, beta1.hat, gamma.hat)
                coef.names <- c(object$ax.names, object$z.names, object$ax.names, object$v.names)

            } else {
                
                U <- cbind(U1, U2)

                d.U <- rbind( cbind(d.U1.beta, d.U1.alpha),
                             cbind(d.U2.beta, d.U2.alpha) ) / nrow(U)

                coefficients <- c(beta.hat, alpha.hat)
                coef.names <- c(object$ax.names, object$z.names)

            }

        }

        names(coefficients) <- coef.names
        
        ## Calculate variance of all estimates
        
        vcov <- robVcov(U, d.U, object$id)

        dimnames(vcov) <- list(coef.names, coef.names)

        result <- list(coefficients = coefficients, vcov = vcov,
                       optim.object = optim.object,
                       optim.object.o = NULL, optim.object.e = NULL)

        return(result)

    }
