fixef.jointModel <-
function (object, process = c("Longitudinal", "Event"), 
        include.splineCoefs = FALSE, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    process <- match.arg(process)
    if (process == "Longitudinal") {
        object$coefficients$betas
    } else {
        gammas <- object$coefficients$gammas
        if (object$method == "ch-Laplace" && !include.splineCoefs) {
            ng <- length(gammas)
            nw <- ncol(object$x$W)
            gammas <- if (is.null(nw)) NULL else gammas[seq(ng - nw + 1, ng)]
        }
        out <- c(gammas, "Assoct" = as.vector(object$coefficients$alpha), 
            "Assoct.s" = as.vector(object$coefficients$Dalpha))
        if (object$method == "weibull-AFT-GH")
            out <- - out
        jj <- grep("Assoct[!^\\.s]", names(out))
        ii <- setdiff(grep("Assoct", names(out)), jj)
        if (length(ii) > 1) {
            nn <- names(object$coefficients$alpha)
            names(out)[ii] <- if (length(nn) == 1) "Assoct" else {
                if (nn[1] == "") 
                    c("Assoct", paste("Assoct", nn[-1], sep = ":"))
                else
                    paste("Assoct", nn, sep = ":")
            }
        }
        if (length(jj) > 1) {
            nn <- names(object$coefficients$Dalpha)
            names(out)[jj] <- if (length(nn) == 1) "Assoct.s" else {
                if (nn[1] == "") 
                    c("Assoct.s", paste("Assoct.s", nn[-1], sep = ":"))
                else
                    paste("Assoct.s", nn, sep = ":")
            }
        }
        if ((lag <- object$y$lag) > 0) {
            kk <- grep("Assoct", names(out), fixed = TRUE)
            names(out)[kk] <- paste(names(out)[kk], "(lag=", lag, ")", sep = "")
        }
        out
    }
}
