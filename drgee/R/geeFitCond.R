geeFitCond <-
    function(y,
             x,
             y.names = colnames(y),
             x.names = colnames(x), 
             link = c("identity", "log", "logit"),
             id,
             rootFinder = findRoots,
             ...) {
        
        link <- match.arg(link)

        if( is.null(x) ) {
            
            return( list(coefficients = NULL,
                         res = NULL, 
                         d.res = NULL, 
                         eq.x = NULL,
                         optim.object = NULL,
                         naive.var = NULL))

        } else if (link == "identity") {

            y.cent <- .Call("center", y, id, PACKAGE = "drgee")
            x.cent <- .Call("center", x, id, PACKAGE = "drgee")

            ## Solve the equation     t(instr) %*% (  y.cent - x.cent %*% beta ) = 0
            ## t(instr) %*% y.cent = t(instr) %*% x.cent %*% beta
            ## beta.hat = solve( crossprod(instr, x.cent) ) %*% t(instr) %*% y.cent

            beta.hat <- as.vector( solve(crossprod(x.cent, x.cent)) %*% crossprod(x.cent, y.cent) )

            return( list(coefficients = beta.hat,
                         res = as.vector(y.cent - x.cent %*% beta.hat),
                         d.res = -x.cent,
                         eq.x = x.cent,
                         optim.object = NULL,
                         naive.var = NULL))

        } else if (link == "log") {

            x.cent <- .Call("center", x, id, PACKAGE = "drgee")
            
            ## Create an initial estimate of beta
            ## with intercept being the mean of the outcome
            intercept.init <- log( colMeans(y) )
            beta.init <- geeFit(y = y, x = cbind(rep(intercept.init, nrow(x)), x), link = link)$coefficients[-1]

            u.func <- function(beta, arg.list) {
                y.star <- arg.list$y * exp(-arg.list$x %*% beta)
                ## y.star.cent <- y.star - ave(y.star, arg.list$id)
                y.star.cent <- .Call("center", y.star, arg.list$id, PACKAGE = "drgee")
                as.vector( crossprod( arg.list$x.cent , y.star.cent))
            }

            all.args <- c(list(beta.init = beta.init, eq.func = u.func, d.eq.func = NULL,
                               arg.list = list(y = y, x = x, x.cent = x.cent, id = id)),
                          list(...))

            ## Call equation solver with beta.init as initial guess and eq.func as estimation function
            root.object <- do.call(rootFinder, all.args)
            beta.hat <- root.object$roots

            y.star.hat <- as.vector(y * exp(-x %*% beta.hat))
            
            d.y.star.hat.beta <- -x * y.star.hat
            
            res  <- .Call("center", as.matrix(y.star.hat), id, PACKAGE = "drgee")
            d.res <- .Call("center", d.y.star.hat.beta, id, PACKAGE = "drgee")

            return( list(coefficients = beta.hat,
                         res = as.vector(res),
                         d.res = d.res,
                         eq.x = x.cent,
                         optim.object = root.object$optim.object,
                         naive.var = NULL) )

        } else if (link == "logit") {

            return( condit(y = y,
                           x = x,
                           y.names = y.names,
                           x.names = x.names, 
                           id = id)
                   )
            
        }

    }
