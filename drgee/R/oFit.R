oFit <-
    function(object, inv = FALSE, ...) {

        if(class(object) != "drgeeData") {
            stop("An object of class \"drgeeData\" is expected")
        }

        if(inv){

            ## For e-estimation with logit link, 
            ## let y and a switch place and replace v with z
            ## and run retrospective logistic regression

            if (object$olink != "logit" | object$elink != "logit")
                stop("\nReverse regression only possible in the logit-logit case\n")
            
            if (object$cond) {
                fit <- geeFitCond(y = object$a, x = cbind(object$yx, object$z),
                link = object$olink, id = object$id, ...)
            } else {
                fit <- geeFit(y = object$a, x = cbind(object$yx, object$z), link = object$olink)
            }

            coef.names <- c(object$yx.names, object$z.names)
            
        } else {
            
            if (object$cond) {
                fit <- geeFitCond(y = object$y, x = cbind(object$ax, object$v),
                link = object$olink, id = object$id, ...)
            } else {
                fit <- geeFit(y = object$y, x = cbind(object$ax, object$v), link = object$olink)
            }

            coef.names <- c(object$ax.names, object$v.names)

        }

        u <- apply(fit$eq.x, 2, '*', fit$res)
        d.u <- crossprod( fit$eq.x , fit$d.res ) / nrow(u)

        coefficients <- fit$coefficients
        names(coefficients) <- coef.names

        vcov <- robVcov(u, d.u, object$id)
        dimnames(vcov) <- list(coef.names, coef.names)


        result <- list(coefficients = coefficients, vcov = vcov, optim.object = fit$optim.object,
                       optim.object.o = fit$optim.object, optim.object.e = NULL)

        return(result)
    }
