ranef.jointModel <-
function (object, type = c("mean", "mode"), postVar = FALSE, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    type <- match.arg(type)
    if (type == "mean") {
        out <- as.matrix(object$EB$post.b)
        rownames(out) <-  if (!object$CompRisk && !object$LongFormat)
            names(object$y$logT)
        else 
            seq_len(nrow(out))
        if (postVar) {
            n <- nrow(out)
            ncz <- ncol(out)
            vars <- vector("list", n)
            for (i in 1:n) {
                vars[[i]] <- matrix(as.matrix(object$EB$post.vb)[i, ], ncz, ncz)
                dimnames(vars[[i]]) <- list(colnames(out), colnames(out))
            }
            names(vars) <- rownames(out)
            attr(out, "postVar") <- vars
        }
    } else {
        mv <- log.posterior.b2(object)
        out <- mv$modes
        if (postVar)
            attr(out, "postVar") <- mv$vars
    }
    out
}
