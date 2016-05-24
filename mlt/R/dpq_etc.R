
.frmt <- function(q) {
    if (is.factor(q)) 
        return(as.character(q))
    return(formatC(q, digits = 3, width = 5))
}

### transformation function
tmlt <- function(object, newdata = object$data, q = NULL, ...) {

    vn <- unlist(variable.names(object$model$model))
    vnx <- vn[!(vn %in% object$response)]
    y <- object$response
    model <- object$model$model

    stopifnot(!is.null(newdata[[y]]) || !is.null(q))

    if (is.data.frame(newdata)) {

        ### unconditional
        if (length(vnx) == 0 & !is.null(q)) {
            newdata <- data.frame(q)
            names(newdata) <- y
            q <- NULL
        }

        ### in sample predictions
        ### this will _not_ work for censored responses
        if (!is.null(newdata[[y]]) & is.null(q)) {
            stopifnot(is.atomic(newdata[[y]]))
            ret <- c(predict(model, newdata = newdata, 
                             coef = coef(object), ...))
            names(ret) <- rownames(newdata)
            ### P(Y \le y_K) = 1 but trafo can be < Inf
            ### depending on parameterisation
            if (is.factor(f <- newdata[[y]])) {
                i <- f == levels(f)[nlevels(f)]
                if (any(i)) 
                    ret[i] <- Inf
            } else {
                ### Y in (b[1], b[2]) => P(Y \le b[1]) = 0, P(Y \le b[2]) = 1
                ### <FIXME> what happens with deriv in ...? </FIXME>
                b <- object$bounds[[y]]
                ret[f < (b[1] - .Machine$double.eps)] <- -Inf
                ret[f > (b[2] - .Machine$double.eps)] <- Inf
            }
            return(ret)
        }

        ### extra quantiles, compute transformation
        ### for each q and each row of newdata
        stopifnot(is.atomic(q))
        stopifnot(length(unique(q)) == length(q))
        dim <- c(length(q), nrow(newdata))

        ### <FIXME> this triggers a trick in 
        ### basefun:::predict.basis; better checks needed </FIXME>
        names(dim) <- c(y, vnx[1])
        newdata <- as.list(newdata)
        newdata[[y]] <- q
        ret <- predict(object$model$model, newdata = newdata, 
                       coef = coef(object), dim = dim, ...)
        dn <- vector(mode = "list", length = 2)
        names(dn) <- c(y, "newdata") ### deparse(substitute(newdata))) ?
        dn[[y]] <- .frmt(q)
        dn[["newdata"]] <- rownames(newdata)
        dimnames(ret) <- dn

        ### trafo of last level is always Inf, see above
        if (is.factor(f <- newdata[[y]])) {
            i <- f == levels(f)[nlevels(f)]
            if (any(i))
                ret[i,] <- Inf
        } else {
            ### Y in (b[1], b[2]) => P(Y \le b[1]) = 0, P(Y \le b[2]) = 1
                ### <FIXME> what happens with deriv in ...? </FIXME>
            b <- object$bounds[[y]]
            ret[f < (b[1] - .Machine$double.eps)] <- -Inf
            ret[f > (b[2] - .Machine$double.eps)] <- Inf
        }
        return(ret)
    }

    ### need to generate newdata outside tmlt such that
    ### the rows of expand.grid(newdata) match the elements of
    ### the return value
    stopifnot(is.atomic(newdata[[y]]))
    stopifnot(is.null(q))
    stopifnot(y %in% names(newdata))
    ret <- predict(object$model$model, newdata = newdata, 
                   coef = coef(object), dim = TRUE, ...)
    dn <- lapply(newdata, .frmt)
    dimnames(ret) <- dn

    ### trafo of last level is always Inf, see above
    if (is.factor(f <- newdata[[y]])) {
        i <- f == levels(f)[nlevels(f)]
        if (any(i)) {
            args <- lapply(names(dn), function(d) {
                if (d == y)
                    return(i)
                return(1:(dim(ret)[which(names(dn) == d)]))
            })
            ret <- do.call("[<-", c(list(i = ret), args, 
                                    list(value = Inf)))
        }
    } else {
        ### Y in (b[1], b[2]) => P(Y \le b[1]) = 0, P(Y \le b[2]) = 1
                ### <FIXME> what happens with deriv in ...? </FIXME>
        b <- object$bounds[[y]]
        i <- (f < (b[1] - .Machine$double.eps))
        if (any(i)) {
            args <- lapply(names(dn), function(d) {
                if (d == y)
                    return(i)
                return(1:(dim(ret)[which(names(dn) == d)]))
            })
            ret <- do.call("[<-", c(list(i = ret), args, 
                                    list(value = -Inf)))
        }
        i <- (f > (b[2] - .Machine$double.eps))
        if (any(i)) {
            args <- lapply(names(dn), function(d) {
                if (d == y)
                    return(i)
                return(1:(dim(ret)[which(names(dn) == d)]))
            })
            ret <- do.call("[<-", c(list(i = ret), args, 
                                    list(value = Inf)))
        }
    }
    return(ret)
}

### distribution function
pmlt <- function(object, newdata = object$data, q = NULL)
    object$model$todistr$p(tmlt(object = object, newdata = newdata,
                                q = q))

### survivor function
smlt <- function(object, newdata = object$data, q = NULL)
    1 - pmlt(object = object, newdata = newdata, q = q)

### cumulative hazard function
Hmlt <- function(object, newdata = object$data, q = NULL)
    -log(smlt(object = object, newdata = newdata, q = q))

### numerical inversion of distribution function
### to get quantile function
.p2q <- function(prob, q, p, interpolate = FALSE, 
                 discrete = FALSE, bounds = c(-Inf, Inf)) {

    prob <- cbind(0, prob, 1)
    i <- rowSums(prob < p)
    if (discrete)
        return(q[i])

    ### Note: mkgrid for bounded variables already contains the bounds
    q <- c(bounds[1], q, bounds[2])

    ### return interval censored quantiles
    if (!interpolate)
        return(R(cleft = q[i], cright = q[i + 1]))

    FIN <- is.finite(q[i]) & is.finite(q[i + 1])

    ptmp <- cbind(prob[cbind(1:nrow(prob), i)],    
                  prob[cbind(1:nrow(prob), i + 1)])
    beta <- (ptmp[,2] - ptmp[,1]) / (q[i + 1] - q[i])
    alpha <- ptmp[,1] - beta * q[i]
    ret <- (p - alpha) / beta
    ret[!is.finite(beta)] <- q[i][!is.finite(beta)]

    if (all(FIN)) return(ret)
    ret[!FIN] <- NA
    left <- q[i]
    left[FIN] <- NA
    right <- q[i + 1]
    right[FIN] <- NA
    return(R(ret, cleft = left, cright = right))
}

### quantile function
qmlt <- function(object, newdata = object$data, p = .5, n = 50, 
                 interpolate = TRUE) {

    y <- object$response
    ### don't accept user-generated quantiles
    q <- mkgrid(object, n = n)[[y]]
    if (!is.null(newdata) & !is.data.frame(newdata)) {
        newdata[[y]] <- NULL
        nm <- names(newdata)
        newdata[[y]] <- q
        newdata <- newdata[c(y, nm)]
        prob <- pmlt(object, newdata)
    } else {
        prob <- pmlt(object, newdata = newdata, q = q)
    } 

    ### convert potential array-valued distribution function
    ### to matrix where rows correspond to observations newdata 
    ### and columns to quantiles q
    ptmp <- t(matrix(prob, nrow = length(q)))
    nr <- nrow(ptmp)
    ptmp <- ptmp[rep(1:nr, each = length(p)),,drop = FALSE]
    pp <- rep(p, nr) ### p varies fastest
    discrete <- !inherits(as.vars(object)[[object$response]],
                          "continuous_var")
    bounds <- bounds(object)[[y]]
    ret <- .p2q(ptmp, q, pp, interpolate = interpolate, bounds = bounds,
                discrete = discrete)

    ### arrays of factors are not allowed
    if (is.factor(q)) return(ret)

    ### return "response" object
    if (!interpolate || inherits(ret, "response")) 
        return(ret)

    dim <- dim(prob)
    dim[1] <- length(p)
    dn <- c(list(p = .frmt(p)), dimnames(prob)[-1])
    return(array(ret, dim = dim, dimnames = dn))
}

### density
dmlt <- function(object, newdata = object$data, q = NULL, log = FALSE) {

    response <- object$data[[y <- object$response]]

    ### Lebesgue density only for double
    if (.type_of_response(response) %in% c("double", "survival")) {
        trafo <- tmlt(object, newdata = newdata, q = q)
        deriv <- 1
        names(deriv) <- y
        trafoprime <- tmlt(object, newdata = newdata, q = q, 
                           deriv = deriv)
        ### <FIXME> trafoprime is +/-Inf at boundaries, so use 0 density
        trafoprime[!is.finite(trafoprime)] <- .Machine$double.eps
        trafoprime <- pmax(.Machine$double.eps, trafoprime)
        if (log)
            return(object$model$todistr$d(trafo, log = TRUE) + log(trafoprime))
        return(object$model$todistr$d(trafo) * trafoprime)
    }

    stopifnot(!is.null(newdata[[y]]) || !is.null(q))

    ### for factors and integers compute density as F(y) - F(y - 1)
    lev <- levels(response)
    if (is.data.frame(newdata)) {

        ### in sample density
        if (!is.null(newdata[[y]]) & is.null(q)) {
            stopifnot(is.atomic(newdata[[y]]))
            q <- newdata[[y]]
            first <- q == lev[1]
            qwoK <- factor(lev[pmax(unclass(q) - 1, 1)], 
                           levels = lev, labels = lev, ordered = is.ordered(q))
            p <- pmlt(object, newdata = newdata)
            newdata[[y]] <- qwoK
            pwoK <- pmlt(object, newdata = newdata)
            pwoK[first] <- 0
            ret <- p - pwoK
        } else {
            ### extra quantiles, compute density
            ### for each q and each row of newdata 
            stopifnot(is.atomic(q))
            first <- q == lev[1]
            qfirst <- q[first]
            qwoK <- q[q != lev[length(lev)]]
            qwo1 <- q[q != lev[1]]

            pfirst <- pmlt(object, newdata = newdata, q = qfirst)
            pwo1 <- pmlt(object, newdata = newdata, q = qwo1)
            pwoK <- pmlt(object, newdata = newdata, q = qwoK)
            ret <- matrix(0, nrow = length(first), ncol = NCOL(pfirst))
            ret[!first,] <- pwo1 - pwoK
            ret[first,] <- pfirst
            rownames(ret) <- as.character(q)
       }
    } else {

        ### need to generate newdata outside tmlt such that
        ### the rows of expand.grid(newdata) match the elements of
        ### the return value
        stopifnot(is.atomic(newdata[[y]]))
        stopifnot(is.null(q))
        stopifnot(y %in% names(newdata))
        dim <- sapply(newdata, NROW)
        q <- newdata[[y]]

        first <- q == lev[1]
        qfirst <- q[first]
        qwoK <- q[q != lev[length(lev)]]
        qwo1 <- q[q != lev[1]]

        newdata[[y]] <- qfirst
        pfirst <- pmlt(object, newdata = newdata)
        newdata[[y]] <- qwo1
        pwo1 <- pmlt(object, newdata = newdata)
        newdata[[y]] <- qwoK
        pwoK <- pmlt(object, newdata = newdata)

        dn <- dim(pfirst)
        names(dn) <- names(dimnames(pfirst))

        frst <- lapply(names(dn), function(d) {
            if (d != y) return(1:dn[d])
            return(first)
        })
        ntfrst <- lapply(names(dn), function(d) {
            if (d != y) return(1:dn[d])
            return(!first)
        })

        dn <- lapply(names(dn), function(d) {
            if (d != y) return(dimnames(pfirst)[[d]])
            return(as.character(q))
        })

        ret <- array(0, dim = dim, dimnames = dn)
        ret <- do.call("[<-", 
            c(list(i = ret), frst, list(value = pfirst)))
        ret <- do.call("[<-", 
            c(list(i = ret), ntfrst, list(value = pwo1 - pwoK)))
    }

    if (log) return(log(ret))
    return(ret)
}

### hazard function
hmlt <- function(object, newdata = object$data, q = NULL, log = FALSE) {
    if (log) 
        return(dmlt(object, newdata = newdata, q = q, log = TRUE) -
               log(smlt(object, newdata = newdata, q = q)))
    return(dmlt(object, newdata = newdata, q = q, log = FALSE) /
           smlt(object, newdata = newdata, q = q))
}
