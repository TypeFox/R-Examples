predict.gv <-
function(object, newdata, ...) {

    if (missing(newdata) || is.null(newdata)) {
        newdata <- object$lag[!is.na(object$lag)]
    }

    mtest <- mtest.gv(object)
    if (!mtest) stop("Object has no model built!")

    model  <- object$model$type
    S <- object$model$sill
    N <- object$model$nugget
    R <- object$model$range    

    if (model == 1) {
    # Gaussian model
        vFUN <- function(x) N + (S-N) * (1-exp(-3*(x**2)/(R**2)))

    } else if (model ==2 ) {
    # Exponential mode
        vFUN <- function(x) N + (S-N) * (1-exp(-x/(R)))

    } else if (model == 3) {
    # Spherical model
        vFUN <- function(x) {
            y <- x * 0 + S
            y[x<R] <- N+(S-N)*(((3*(x[x<R]))/(2*R))-((x[x<R])**3/(2*R**3)))
            return(y)
        }
    } else if (model == 4) {
        Sl <- object$model$slope
        if (is.infinite(S)) {
            # Linear model
            vFUN <- function(x) N + Sl*x
        } else {
            # Linear model with sill
            vFUN <- function(x) {
                y <- x * 0 + S
                y[x<R] <- N + Sl * x[x<R]
                return(y)
            }
        }            
    }

    y <- vFUN(newdata)
    y
}
