
R <- function(object, ...)
    UseMethod("R")

R.Surv <- function(object, ...) {

    type <- attr(object, "type")
    stopifnot(type %in% c("left", "right", "interval", 
                          "interval2", "counting"))
    status <- object[, "status"]

    ret <- switch(type,
        "right" = R(object = ifelse(status == 1, object[, "time"], NA),
                    cleft = ifelse(status != 1, object[, "time"], NA),
                    cright = ifelse(status != 1, Inf, NA)),
        "left" =  R(object = ifelse(status == 1, object[, "time"], NA),
                    cleft = ifelse(status != 1, -Inf, NA),
                    cright = ifelse(status != 1, object[, "time"], NA)),

        "interval2" = {
            ret <- cbind(left = object[, "time1"], 
                         right = object[, "time2"])
            ret$left[is.na(ret$left)] <- -Inf
            ret$right[is.na(ret$right)] <- Inf
            R(cleft = ret$left, cright = ret$right)
        },
        "interval" = {
            status <- factor(status, levels = 0:3, 
                             labels = c("right", "exact", "left", "interval"))
            tmp <- matrix(NA, nrow = nrow(object), ncol = 2)
            colnames(tmp) <- c("left", "right")
            for (s in levels(status)) {
                idx <- which(status == s)
                tmp[idx, ] <- switch(s, 
                    "right" = cbind(object[idx, "time1"], Inf),
                    "exact" = cbind(object[idx, "time1"], NA),
                    "left" = cbind(-Inf, object[idx, "time1"]),
                    "interval" = object[idx, c("time1", "time2")])
            }
            R(object = ifelse(is.na(tmp[, "right"]), tmp[, "left"], NA), 
              cleft = ifelse(is.na(tmp[, "right"]), NA, tmp[, "left"]), 
              cright = tmp[, "right"])
        },
        ### left truncation, right censoring
        "counting" = R(object = ifelse(status == 1, object[, "stop"], NA),
                       cleft = ifelse(status != 1, object[, "stop"], NA),
                       cright = ifelse(status != 1, Inf, NA),
                       tleft = object[, "start"])
    )
    attr(ret, "prob") <- function(weights) {
        w <- weights[weights > 0]
        sf <- survival::survfit(object ~ 1, subset = weights > 0, weights = w)
        function(y) {
            s <- summary(sf, times = y)$surv[order(y)]
            s[is.na(s)] <- 0
            1 - s
        }
    }
    ret
}


R.factor <- function(object, ...) {

    warning("response is unordered factor;
             results may depend on order of levels")
    return(R(as.ordered(object), ...))
}

R.ordered <- function(object, cleft = NA, cright = NA, ...) {

    lev <- levels(object)
    ret <- .mkR(exact = object, cleft = cleft, cright = cright, ...)
    ret[is.na(ret$cright), "cright"] <- ret$exact[is.na(ret$cright)]
    ret[is.na(ret$cleft), "cleft"] <- factor(unclass(object)[is.na(ret$cleft)] - 1,
        levels = 1:length(lev), labels = lev, exclude = 0, ordered = TRUE)
    ret$exact <- NA
    ret[ret$cright == lev[nlevels(object)], "cright"] <- NA
    attr(ret, "prob") <- function(weights) {
        prt <- prop.table(xtabs(weights ~ object))
        function(y) prt[y]
    }
    ret
}

R.integer <- function(object, cleft = NA, cright = NA, bounds = c(0L, Inf), ...) {

    ret <- .mkR(exact = object, cleft = cleft, cright = cright, ...)
    ret$cright[is.na(ret$cright)] <- ret$exact[is.na(ret$cright)]
    ret$cright[ret$cright == bounds[2]] <- NA
    ret$cleft[is.na(ret$cleft)] <- ret$exact[is.na(ret$cleft)] - 1
    ret$cleft[ret$cleft < bounds[1]] <- NA
    ret$exact <- NA
    attr(ret, "prob") <- function(weights)
        .wecdf(object, weights)
    ret
}

### handle exact integer / factor as interval censored
R.numeric <- function(object = NA, cleft = NA, cright = NA, 
                      tleft = NA, tright = NA, tol = sqrt(.Machine$double.eps), 
                      ...) {

    ### treat extremely small intervals as `exact' observations
    d <- cright - cleft
    if (any(!is.na(d) | is.finite(d))) {
        if (any(d < 0, na.rm = TRUE)) stop("cleft > cright")
        i <- (d < tol)
        if (any(i, na.rm = TRUE)) {
            object[i] <- cleft[i]
            cleft[i] <- cright[i] <- NA
        }
    }
    ret <- .mkR(exact = object, cleft = cleft, cright = cright,
                tleft = tleft, tright = tright)
    attr(ret, "prob") <- function(weights)
         .wecdf(object, weights)
    ret
}

### for object = NA, ie censored observations only
R.logical <- R.numeric

R.response <- function(object, ...)
    return(object)

R.default <- function(object, ...)
    stop("cannot deal with response class", class(object))

.mkR <- function(...) {

    args <- list(...)
    cls <- unique(sapply(args, function(a) class(a)[1]))
    cls <- cls[cls != "logical"]
    stopifnot(length(cls) <= 1)

    n <- unique(sapply(args, length))
    stopifnot(length(n) <= 2)
    if (length(n) == 2) stopifnot(min(n) == 1)

    ret <- do.call("as.data.frame", list(x = args))
    if (is.null(ret$exact)) ret$exact <- NA
    if (is.null(ret$cleft) || all(is.na(ret$cleft))) {
        ret$cleft <- NA
        if (is.ordered(ret$exact)) 
            ret$cleft <- factor(ret$cleft, levels = 1:nlevels(ret$exact), 
                                labels = levels(ret$exact), ordered = TRUE)
    }
    if (is.null(ret$cright) || all(is.na(ret$cright))) {
        ret$cright <- NA
        if (is.ordered(ret$exact)) 
            ret$cright <- factor(ret$cright, levels = 1:nlevels(ret$exact), 
                                 labels = levels(ret$exact), ordered = TRUE)
    }
    if (is.null(ret$tleft)) ret$tleft <- NA
    if (is.null(ret$tright)) ret$tright <- NA

    if (all(is.finite(ret$exact))) {
#        ret$rank <- rank(ret$exact, ties.method = "max")
        ret$approxy <- ret$exact
    } else {
        ### some meaningful ordering of observations
        tmpexact <- as.numeric(ret$exact)
        tmpleft <- as.numeric(ret$cleft)
        tmpright <- as.numeric(ret$cright)
        tmpler <- c(tmpleft, tmpexact, tmpright)
        tmpleft[!is.finite(tmpleft)] <- min(tmpler[is.finite(tmpler)])
        tmpright[!is.finite(tmpright)] <- max(tmpler[is.finite(tmpler)])
        tmpexact[is.na(tmpexact)] <-
            (tmpleft + ((tmpright - tmpleft) / 2))[is.na(tmpexact)]
#        ret$rank <- rank(tmpexact, ties.method = "max")   
        ret$approxy <- tmpexact
    }
    class(ret) <- c("response", class(ret))
    ret[, c("tleft", "cleft", "exact", "cright", "tright", "approxy"), drop = FALSE]
}

.exact <- function(object)
    !is.na(object$exact)

.cleft <- function(object)
    is.finite(object$cleft) 

.cright <- function(object)
    is.finite(object$cright)

.cinterval <- function(object)
    .cleft(object) | .cright(object)

.tleft <- function(object)
    is.finite(object$tleft) 

.tright <- function(object)
   is.finite(object$tright)

.tinterval <- function(object)
    .tleft(object) | .tright(object)

.mm_exact <- function(model, data, response, object) {

    e <- .exact(object)
    if (!any(e)) return(NULL)
    tmp <- data[e,,drop = FALSE]
    tmp[[response]] <- object$exact[e]
    Y <- model.matrix(model, data = tmp)
    deriv <- 1
    names(deriv) <- response
    Yprime <- model.matrix(model, data = tmp, deriv = deriv)

    trunc <- NULL
    if (any(.tinterval(object) & e)) {
        Ytleft <- matrix(-Inf, nrow = nrow(Y), ncol = ncol(Y))
        Ytright <- matrix(Inf, nrow = nrow(Y), ncol = ncol(Y))
        if (any(il <- (.tleft(object) & e))) {
            tmp <- data[il,]
            tmp[[response]] <- object$tleft[il]
            Ytleft[il,] <- model.matrix(model, data = tmp)
        }
        if (any(ir <- (.tright(object) & e))) {
            tmp <- data[ir,,drop = FALSE]
            tmp[[response]] <- object$tright[ir]
            Ytright[ir,] <- model.matrix(model, data = tmp)
        }
        trunc <- list(left = Ytleft, right = Ytright)
    }

    list(Y = Y, Yprime = Yprime, trunc = trunc, which = which(e))
}

.mm_interval <- function(model, data, response, object) {

    i <- .cinterval(object)
    if (!any(i)) return(NULL)
    tmpdata <- data[i,,drop = FALSE]
    object <- object[i,, drop = FALSE]

    Yleft <- NULL
    if (any(il <- .cleft(object))) {
        tmp <- tmpdata[il,,drop = FALSE]
        tmp[[response]] <- object$cleft[il]
        Ytmp <- model.matrix(model, data = tmp)
        Yleft <- matrix(-Inf, nrow = length(il), ncol = ncol(Ytmp))
        colnames(Yleft) <- colnames(Ytmp)
        attr(Yleft, "constraint") <- attr(Ytmp, "constraint")
        attr(Yleft, "Assign") <- attr(Ytmp, "Assign")
        Yleft[il,] <- Ytmp
    }

    Yright <- NULL
    if (any(ir <- .cright(object))) {
        tmp <- tmpdata[ir,, drop = FALSE]
        tmp[[response]] <- object$cright[ir]
        Ytmp <- model.matrix(model, data = tmp)
        Yright <- matrix(Inf, nrow = length(ir), ncol = ncol(Ytmp))
        colnames(Yright) <- colnames(Ytmp)
        attr(Yright, "constraint") <- attr(Ytmp, "constraint")
        attr(Yright, "Assign") <- attr(Ytmp, "Assign")
        Yright[ir,] <- Ytmp
    }

    if (is.null(Yright)) { 
        Yright <- matrix(Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        colnames(Yright) <- colnames(Yleft)
        attr(Yright, "constraint") <- attr(Yleft, "constraint")
        attr(Yright, "Assign") <- attr(Yleft, "Assign")
    }
    if (is.null(Yleft)) {
        Yleft <- matrix(-Inf, nrow = nrow(Yright), ncol = ncol(Yright))
        colnames(Yleft) <- colnames(Yright)
        attr(Yleft, "constraint") <- attr(Yright, "constraint")
        attr(Yleft, "Assign") <- attr(Yright, "Assign")
    }

    trunc <- NULL
    if (any(.tinterval(object))) {
        Ytleft <- matrix(-Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        Ytright <- matrix(Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        colnames(Ytleft) <- colnames(Ytright) <- colnames(Yleft)
        if (any(il <- (.tleft(object)))) {
            tmp <- tmpdata[il,,drop = FALSE]
            tmp[[response]] <- object$cleft[il]
            Ytleft[il,] <- model.matrix(model, data = tmp)
        }
        if (any(ir <- (.tright(object)))) {
            tmp <- tmpdata[ir,,drop = FALSE]
            tmp[[response]] <- object$cleft[ir]
            Ytright[ir,] <- model.matrix(model, data = tmp)
        }
        trunc <- list(left = Ytleft, right = Ytright)
    }

    list(Yleft = Yleft, Yright = Yright, trunc = trunc, which = which(i))
}

.wecdf <- function(x, weights) {
    ### from: spatstat::ewcdf
    ox <- order(x)
    x <- x[ox]
    w <- weights[ox]
    vals <- sort(unique(x))
    xmatch <- factor(match(x, vals), levels = seq_along(vals))
    wmatch <- tapply(w, xmatch, sum)
    wmatch[is.na(wmatch)] <- 0
    cumwt <- cumsum(wmatch) / sum(wmatch)
    approxfun(vals, cumwt, method = "constant", yleft = 0, 
              yright = sum(wmatch), f = 0, ties = "ordered")
}
