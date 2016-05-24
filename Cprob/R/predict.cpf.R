predict.cpf <- function(object, timepoints, ...) {
    if (!inherits(object, "cpf")) {
        stop("'object' must be of class 'cpf'")
    }
    if (sum(timepoints < 0) >= 1) {
        stop("negative timepoints' may be problematic")
    }
    ns <- length(object$size.strata)
    size <- c(0, cumsum(object$size.strata))
    pred <-lapply(seq_len(ns), function(i) {
        times <- object$time[(size[i] + 1):size[i + 1]]
        ind <- findInterval(timepoints, times)
        cp <- object$cp[(size[i] + 1):size[i + 1]][ind]
        var <- object$var[(size[i] + 1):size[i + 1]][ind]
        lower <- object$lower[(size[i] + 1):size[i + 1]][ind]
        upper <- object$upper[(size[i] + 1):size[i + 1]][ind]
        n.risk <- object$n.risk[(size[i] + 1):size[i + 1]][ind]
        if (!is.null(object$X)) {
            group <- rep(object$X[i, ], length(timepoints))
            data.frame(time = timepoints, cp = cp, var = var,
                       lower = lower, upper = upper, n.risk = n.risk,
                       group = group)
        }
        else {
            data.frame(time = timepoints, cp = cp, var = var,
                       lower = lower, upper = upper, n.risk = n.risk)
        }
    })
    pred <- do.call(rbind, pred)
    pred
}
