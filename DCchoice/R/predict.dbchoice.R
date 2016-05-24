predict.dbchoice <- function(object, newdata = NULL, 
    type = c("utility", "probability"), bid = NULL, ...)
{
    type <- match.arg(type)
    COEF <- matrix(object$coefficients, ncol = 1)
    dist <- object$distribution

    if (!is.null(newdata)) {
        vars <- attr(object$terms, "term.labels")
        nvars <- length(vars)
        vars <- vars[-c(1, 2, nvars)]
        new.formula <- as.formula(paste("~ ", paste(vars, collapse= "+")))
        mf.newX <- model.frame(new.formula, newdata, xlev = object$xlevels)
        mm.newX <- model.matrix(new.formula, mf.newX, contrasts.arg = object$contrasts)

        V <- mm.newX %*% COEF
        V <- as.vector(V)

        if(type == "probability") {
            if(dist == "logistic" | dist == "log-logistic") {
                P <- plogis(-V, lower.tail = FALSE, log.p = FALSE)
            } else if (dist == "normal" | dist == "log-normal") {
                P <- pnorm(-V, lower.tail = FALSE, log.p = FALSE)
            } else if (dist == "weibull") {
                P <- pweibull(exp(-V), shape = 1, lower.tail = FALSE, log.p = FALSE)
            }
            P <- as.vector(P)
        }
    } else {
        bid.names <- colnames(bid)
        key.var <- bid.names[1]
        mf.X <- object$data.name
        mm.X <- model.matrix(object$formula, data = mf.X, rhs = 1:2)
        mm.X <- mm.X[, -c(ncol(mm.X))]
        covariate.names <- colnames(mm.X)[-c(ncol(mm.X))]

        tmp.sort.id <- c(1:nrow(mm.X))
        mm.X <- cbind(mm.X, tmp.sort.id)
        mm.X <- merge(bid, mm.X, by = key.var)
        mm.X <- mm.X[order(mm.X$tmp.sort.id), ]
        mm.X <- subset(mm.X, select = -tmp.sort.id)

        ivar.names.mm.XF <- c(covariate.names, bid.names[1])
        ivar.names.mm.XU <- c(covariate.names, bid.names[2])
        ivar.names.mm.XL <- c(covariate.names, bid.names[3])
        mm.XF <- as.matrix(subset(mm.X, select = ivar.names.mm.XF))
        mm.XU <- as.matrix(subset(mm.X, select = ivar.names.mm.XU))
        mm.XL <- as.matrix(subset(mm.X, select = ivar.names.mm.XL))

        VF <- mm.XF %*% COEF
        VU <- mm.XU %*% COEF
        VL <- mm.XL %*% COEF
        V <- cbind(VF, VU, VL)
        colnames(V) <- c("f", "u", "l")
        rownames(V) <- NULL

        if(dist == "logistic" | dist == "log-logistic") {
            P.yy = plogis(-VU, lower.tail = FALSE, log.p = FALSE)
            P.nn = plogis(-VL, lower.tail = TRUE,  log.p = FALSE)
            P.yn = plogis(-VU, lower.tail = TRUE,  log.p = FALSE) -
                   plogis(-VF, lower.tail = TRUE,  log.p = FALSE)
            P.ny = plogis(-VF, lower.tail = TRUE,  log.p = FALSE) -
                   plogis(-VL, lower.tail = TRUE,  log.p = FALSE)
        } else if (dist == "normal" | dist == "log-normal") {
            P.yy = pnorm(-VU, lower.tail = FALSE, log.p = FALSE)
            P.nn = pnorm(-VL, lower.tail = TRUE,  log.p = FALSE)
            P.yn = pnorm(-VU, lower.tail = TRUE,  log.p = FALSE) -
                   pnorm(-VF, lower.tail = TRUE,  log.p = FALSE)
            P.ny = pnorm(-VF, lower.tail = TRUE,  log.p = FALSE) -
                   pnorm(-VL, lower.tail = TRUE,  log.p = FALSE)
        } else if (dist == "weibull") {
            P.yy = pweibull(exp(-VU), shape = 1, lower.tail = FALSE, log.p = FALSE)
            P.nn = pweibull(exp(-VL), shape = 1, lower.tail = TRUE,  log.p = FALSE)
            P.yn = pweibull(exp(-VU), shape = 1, lower.tail = TRUE,  log.p = FALSE) -
                   pweibull(exp(-VF), shape = 1, lower.tail = TRUE,  log.p = FALSE)
            P.ny = pweibull(exp(-VF), shape = 1, lower.tail = TRUE,  log.p = FALSE) -
                   pweibull(exp(-VL), shape = 1, lower.tail = TRUE,  log.p = FALSE)
        }

        P <- cbind(P.yy, P.nn, P.yn, P.ny)
        colnames(P) <- c("yy", "nn", "yn", "ny")
        rownames(P) <- NULL
    }

    if(type == "utility") {
        return(V)
    } else {
        return(P)
    }
}
