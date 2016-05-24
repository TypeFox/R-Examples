
starting <- function(data) {
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    fam.clear()
    for (i in seq(along = data$families))
        fam.set(data$families[[i]])
    out <- .C("aster_starting_theta",
        nnode = length(data$regroup),
        group = as.integer(data$regroup),
        code = as.integer(data$recode),
        theta = double(length(data$regroup)),
        PACKAGE = "aster2")
    fam.clear()
    return(out$theta)
}

inner.loop <- function(x, gradfun, hessfun, typ.grad.size,
    fudge, tol1, tol2, lambda, old.search.dir, old.gradient,
    old.inverse.hessian, old.point,
    type = c("steep", "polak-ribere", "newton", "quasi-newton")) {

    type <- match.arg(type)

    stopifnot(is.atomic(x))
    stopifnot(is.numeric(x))
    stopifnot(is.finite(x))
    stopifnot(is.function(gradfun))
    stopifnot(is.function(hessfun))

    # done in outer loop, omit here
    # gradient <- gradfun(x)
    # if (! is.atomic(gradient))
    #     stop("gradfun does not return atomic vector")
    # if (! is.numeric(gradient))
    #     stop("gradfun does not return numeric vector")
    # if (length(gradient) != length(x))
    #     stop("gradfun does not return vector of same length as x")
    # if (! all(is.finite(gradient)))
    #     stop("gradfun does not return vector with all components finite")
    # 
    # hessian <- hessfun(x)
    # if (! is.atomic(hessian))
    #     stop("hessfun does not return atomic matrix")
    # if (! is.numeric(hessian))
    #     stop("hessfun does not return numeric matrix")
    # if (! is.matrix(hessian))
    #     stop("hessfun does not return numeric matrix")
    # if (! all.equal(dim(hessian), rep(length(x), 2)))
    #     stop("hessing does not return d by d matrix where d == length(x)")
    # if (! all(is.finite(hessian)))
    #     stop("gradfun does not return matrix with all components finite")

    if (missing(typ.grad.size)) {
        base <- rep(1, length(x))
    } else {
        stopifnot(is.atomic(typ.grad.size))
        stopifnot(is.numeric(typ.grad.size))
        stopifnot(is.finite(typ.grad.size))
        base <- abs(typ.grad.size)
        base <- pmax(base, max(base) * sqrt(.Machine$double.eps))
    }
    gradient <- gradfun(x)
    rel.norm.grad <- sqrt(sum((gradient / base)^2))

    if (rel.norm.grad < .Machine$double.eps * 8) {
        #### already converged ####
        return(list(x = x, gradient = gradient, type = "none", iter = 0,
            rel = rel.norm.grad, lambda = as.double(NA)))
    }

    if (rel.norm.grad < sqrt(.Machine$double.eps)) {
        #### do Newton without line search regardless of type ####
        type <- "pure-newton"
    }

    if (type == "quasi-newton") {
        if (missing(old.point) || missing(old.gradient)
            || missing(old.inverse.hessian)) {
            type <- "steep"
        } else {
            stopifnot(is.atomic(old.point))
            stopifnot(is.numeric(old.point))
            stopifnot(is.finite(old.point))
            stopifnot(is.atomic(old.gradient))
            stopifnot(is.numeric(old.gradient))
            stopifnot(is.finite(old.gradient))
            stopifnot(is.atomic(old.inverse.hessian))
            stopifnot(is.numeric(old.inverse.hessian))
            stopifnot(is.finite(old.inverse.hessian))
            stopifnot(is.matrix(old.inverse.hessian))
            stopifnot(length(old.point) == length(old.gradient))
            stopifnot(length(old.point) == nrow(old.inverse.hessian))
            stopifnot(length(old.point) == ncol(old.inverse.hessian))
            delta.point <- x - old.point
            delta.gradient <- gradient - old.gradient
            delta.dot <- sum(delta.point * delta.gradient)
            foo <- as.numeric(old.inverse.hessian %*% delta.gradient)
            bar <- sum(delta.gradient * foo)
            baz <- outer(delta.point, delta.gradient)
            qux <- outer(delta.point, delta.point)
            inverse.hessian <- old.inverse.hessian -
                (baz + t(baz)) / delta.dot +
                (1 + bar / delta.dot) * (qux / delta.dot)
            search.dir <- as.numeric(inverse.hessian %*% gradient)
            if (sum(search.dir * gradient) < 0) {
                type <- "steep"
                search.dir <- gradient
            }
        }
    }
    if (type == "polak-ribere") {
        if (missing(old.search.dir) || missing(old.gradient)) {
            type <- "steep"
        } else {
            stopifnot(is.atomic(old.gradient))
            stopifnot(is.numeric(old.gradient))
            stopifnot(is.finite(old.gradient))
            stopifnot(is.atomic(old.search.dir))
            stopifnot(is.numeric(old.search.dir))
            stopifnot(is.finite(old.search.dir))
            gammaPR <- sum((gradient - old.gradient) * gradient) /
                sum(old.gradient^2)
            if (gammaPR <= 0) {
                gammaPR <- 0
                type <- "steep"
            }
            search.dir <- gradient + gammaPR * old.search.dir
        }
    }
    if (type == "steep") {
        search.dir <- gradient
    }
    if (type == "newton" || type == "pure-newton") {
        hessian <- hessfun(x)
        search.dir <- solve(hessian, - gradient)
    }

    if (type == "pure-newton") {
        x.lambda <- x + search.dir
        gradient.lambda <- gradfun(x.lambda)
        rel.norm.grad <- sqrt(sum((gradient.lambda / base)^2))
        return(list(x = x.lambda, gradient = gradient.lambda, type = type,
            iter = 1, rel = rel.norm.grad, lambda = 1,
            search.dir = search.dir, old.gradient = gradient, old.point = x))
    }

    ##### line search #####

    criterion0 <- sum(gradient * search.dir)

    lambda.low <- 0
    lambda.hig <- Inf
    criterion.low <- criterion0
    criterion.hig <- (- Inf)
    if (missing(lambda)) {
        lambda <- 1
    }

    iter <- 1
    repeat {

        x.lambda <- x + lambda * search.dir
        gradient.lambda <- try(gradfun(x.lambda), silent = TRUE)

        if (inherits(gradient.lambda, "try-error")) {
            criterion.lambda <- (- Inf)
        } else {
            criterion.lambda <- sum(gradient.lambda * search.dir)
        }

        if (criterion.lambda < (- tol1 * criterion0)) {
            lambda.hig <- lambda
            criterion.hig <- criterion.lambda
            foo <- criterion.low / (criterion.low - criterion.hig)
            lambda <- ((foo + fudge) * lambda.hig +
                (1 - foo + fudge) * lambda.low) / (1 + 2 * fudge)
        } else if (criterion.lambda > tol2 * criterion0) {
            lambda.low <- lambda
            criterion.low <- criterion.lambda
            if (lambda.hig == Inf) {
                lambda <- 2 * lambda
            } else {
                foo <- criterion.low / (criterion.low - criterion.hig)
                lambda <- ((foo + fudge) * lambda.hig +
                    (1 - foo + fudge) * lambda.low) / (1 + 2 * fudge)
            }
        } else {
            break
        }

        iter <- iter + 1
    }

    rel.norm.grad <- sqrt(sum((gradient.lambda / base)^2))
    return(list(x = x.lambda, gradient = gradient.lambda, type = type,
        iter = iter, rel = rel.norm.grad, lambda = lambda,
        search.dir = search.dir, old.gradient = gradient, old.point = x,
        inverse.hessian =
        if(exists("inverse.hessian")) inverse.hessian else diag(length(x)),
        hessian = if(exists("hessian")) hessian else NULL))
}

outer.loop <- function(x, gradfun, hessfun, typ.grad.size, fudge, tol1, tol2,
    lambda, blather = FALSE) {

    stopifnot(is.atomic(x))
    stopifnot(is.numeric(x))
    stopifnot(is.finite(x))
    stopifnot(is.function(gradfun))
    stopifnot(is.function(hessfun))

    gradient <- gradfun(x)
    if (! is.atomic(gradient))
        stop("gradfun does not return atomic vector")
    if (! is.numeric(gradient))
        stop("gradfun does not return numeric vector")
    if (length(gradient) != length(x))
        stop("gradfun does not return vector of same length as x")
    if (! all(is.finite(gradient)))
        stop("gradfun does not return vector with all components finite")

    hessian <- hessfun(x)
    if (! is.atomic(hessian))
        stop("hessfun does not return atomic matrix")
    if (! is.numeric(hessian))
        stop("hessfun does not return numeric matrix")
    if (! is.matrix(hessian))
        stop("hessfun does not return numeric matrix")
    if (! all.equal(dim(hessian), rep(length(x), 2)))
        stop("hessing does not return d by d matrix where d == length(x)")
    if (! all(is.finite(hessian)))
        stop("gradfun does not return matrix with all components finite")

    if (missing(typ.grad.size)) {
        typ.grad.size <- rep(1, length(x))
    } else {
        stopifnot(is.atomic(typ.grad.size))
        stopifnot(is.numeric(typ.grad.size))
        stopifnot(is.finite(typ.grad.size))
    }

    stopifnot(is.atomic(fudge))
    stopifnot(is.numeric(fudge))
    stopifnot(is.finite(fudge))
    stopifnot(length(fudge) == 1)
    stopifnot(fudge > 0)
    stopifnot(is.atomic(tol1))
    stopifnot(is.numeric(tol1))
    stopifnot(is.finite(tol1))
    stopifnot(length(tol1) == 1)
    stopifnot(tol1 > 0)
    stopifnot(is.atomic(tol2))
    stopifnot(is.numeric(tol2))
    stopifnot(is.finite(tol2))
    stopifnot(length(tol2) == 1)
    stopifnot(tol2 > 0)

    if (missing(lambda)) {
        lambda <- 1
    } else {
        stopifnot(is.atomic(lambda))
        stopifnot(is.numeric(lambda))
        stopifnot(is.finite(lambda))
        stopifnot(length(lambda) == 1)
        stopifnot(lambda > 0)
    }

    stopifnot(is.atomic(blather))
    stopifnot(is.logical(blather))
    stopifnot(length(blather) == 1)

    d <- length(x)
    out <- list(x = x, lambda = lambda)

    if (blather) {
        x.save <- x
        gradient <- gradfun(x)
        gradient.save <- gradient
        lambda.save <- as.double(NA)
        base <- abs(typ.grad.size)
        base <- pmax(base, max(base) * sqrt(.Machine$double.eps))
        rel.norm.grad <- sqrt(sum((gradient / base)^2))
        rel.save <- rel.norm.grad
        iter.save <- 0
        type.save <- "start"

    }

    for (i in seq(1, d / 2)) {

        out <- inner.loop(out$x, gradfun, hessfun, typ.grad.size, fudge,
            tol1, tol2, lambda = out$lambda, type = "steep")

        if (blather) {
            x.save <- rbind(x.save, out$x, deparse.level = 0)
            gradient.save <- rbind(gradient.save, out$gradient,
                deparse.level = 0)
            lambda.save <- c(lambda.save, out$lambda)
            rel.save <- c(rel.save, out$rel)
            iter.save <- c(iter.save, out$iter)
            type.save <- c(type.save, out$type)
        }
    }

    for (i in seq(1, 2 * d)) {

        out <- inner.loop(out$x, gradfun, hessfun, typ.grad.size, fudge,
            tol1, tol2, lambda = 1, type = "newton")

        if (blather) {
            x.save <- rbind(x.save, out$x, deparse.level = 0)
            gradient.save <- rbind(gradient.save, out$gradient,
                deparse.level = 0)
            lambda.save <- c(lambda.save, out$lambda)
            rel.save <- c(rel.save, out$rel)
            iter.save <- c(iter.save, out$iter)
            type.save <- c(type.save, out$type)
        }

        if (out$type %in% c("pure-newton", "none")) {
             break
        }
    }

    if (blather) {
        return(list(x = out$x, gradient = out$gradient, hessian = out$hessian,
            blather = list(x = x.save, gradient = gradient.save,
            lambda = lambda.save, relative.error = rel.save,
            iter = iter.save, type = type.save),
            converged = out$type %in% c("pure-newton", "none")))
    } else {
        return(list(x = out$x, gradient = out$gradient, hessian = out$hessian,
            converged = out$type %in% c("pure-newton", "none")))
    }
}

