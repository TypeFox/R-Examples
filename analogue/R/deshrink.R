`deshrink` <- function(env, wa.env,
                       type = c("inverse", "classical",
                       "expanded", "none","monotonic")) {
### Inline Functions ############################################
    ## inverse deshrinking
    `inverse` <- function(env, wa.env) {
        X <- cbind(rep(1, length(wa.env)), wa.env)
        QR <- qr(X)
        coef <- qr.coef(QR, env)
        pred <- qr.fitted(QR, env)
        return(list(coefficients = coef, env = pred))
    }
    ## classical deshrinking
    `classical` <- function(env, wa.env) {
        X <- cbind(rep(1, length(env)), env)
        QR <- qr(X)
        coef <- drop(qr.coef(QR, wa.env))
        ##coef <- c(-coef[1], 1)/coef[2]
        ##pred <- deshrink.pred(wa.env, coef)
        pred <- (wa.env - coef[1]) / coef[2]
        return(list(coefficients = coef, env = pred))
    }
    ## deshrink to equal SD
    ## A bit like in vegan's wascores, but wascores uses
    ## weighted sd which would need row and column sums in the
    ## function call, and this would make the function API
    ## incompatible with other *.deshrink functions.
    `expanded` <- function(env, wa.env) {
        b1 <- sd(env)/sd(wa.env)
        b0 <- mean(env) - b1 * mean(wa.env)
        pred <- b0 + b1 * wa.env
        return(list(coefficients = c(b0, b1), env = pred))
    }
    ## No deshrinking
    ## Do not deshrink: for those who think they know  what they
    ## are doing
    `none` <- function(env, wa.env) {
        return(list(coefficients = c(0, 1), env = wa.env))
    }
    ## Monotonic deshrinking via pcls() in mgcv
    ## Use a spline fit to deshrink instead of the linear inverse
    ## or classical deshrinking methods
    ## Needs to be constrained to be monotonic so ?mono.con
    `monotonic` <- function(env, wa.env) {
        df <- data.frame(env = env, wa.env = drop(wa.env))
        mod <- gam(env ~ s(wa.env, k = 10, bs = "cr"), data = df)
        ## grab bits for setting up a monotonic spline, see ?pcls
        sm <- smoothCon(s(wa.env, k = 10, bs = "cr"), data = df,
                        knots = NULL)[[1]]
        ## Fm are the constraints to enforce monotonicity
        Fm <- mono.con(sm$xp)
        G <- list(X = sm$X, C = matrix(0,0,0), sp = mod$sp, p = sm$xp,
                  y = env, w = env*0+1, Ain = Fm$A, bin = Fm$b, S = sm$S,
                  off = 0)
        p <- pcls(G)
        ## predict for the current data
        pred <- Predict.matrix(sm, data = data.frame(wa.env = wa.env)) %*% p
        pred <- drop(pred)
        ## return
        list(coefficients = list(sm = sm, p = p), env = pred)
    }
### End Inline Functions #########################################
    if(missing(type))
        type <- "inverse"
    type <- match.arg(type)
    res <- switch(type,
                  inverse = inverse(env, wa.env),
                  classical = classical(env, wa.env),
                  expanded = expanded(env, wa.env),
                  none = none(env, wa.env),
                  monotonic = monotonic(env, wa.env))
    class(res) <- c("deshrink","list")
    attr(res, "type") <- type
    return(res)
}
