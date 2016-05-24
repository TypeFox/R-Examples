summary.cpf <- function(object, ...) {
    if (!inherits(object, "cpf"))
        stop("'object' must be of class 'cpf'")
    if (is.null(object$X)) {
         temp <- list(data.frame(time = object$time,
                                 n.risk = object$n.risk,
                                 n.event = object$n.event[, 1],
                                 n.event.other = object$n.event[, 2],
                                 cp = object$cp,
                                 var.cp = object$var,
                                 lower = object$lower,
                                 upper = object$upper))
     }
    else {
        size <- c(0, cumsum(object$size.strata))
        cova <- sapply(1:length(object$size.strata), function(i) {
            paste(colnames(object$X), object$X[i, ], sep = " = ")
        })
        temp <- lapply(seq_along(cova), function(i) {
            tmp <- data.frame(time = object$time[(size[i] + 1):size[i + 1]],
                              n.risk = object$n.risk[(size[i] + 1):size[i + 1]],
                              n.event = object$n.event[(size[i] + 1):size[i + 1], 1],
                              n.event.other = object$n.event[(size[i] + 1):size[i + 1], 2],
                              cp = object$cp[(size[i] + 1):size[i + 1]],
                              var = object$var[(size[i] + 1):size[i + 1]],
                              lower = object$lower[(size[i] + 1):size[i + 1]],
                              upper = object$upper[(size[i] + 1):size[i + 1]])
            tmp
        })
        names(temp) <- cova
    }
    res <- list(est = temp, call = object$call, X = object$X, z = object$z, p = object$p)
    class(res) <- "summary.cpf"
    res
}
