coef.ltm <-
function (object, standardized = FALSE, prob = FALSE, order = FALSE, ...) {
    if (!inherits(object, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    cof <- object$coef
    cof <- if (standardized) {
        factors <- object$ltst$factors
        if (any(object$ltst$inter, object$ltst$quad.z1, object$ltst$quad.z2)) {
            warning("standardized loadings are returned for the simple one- and two-factor models.\n")
            colnames(cof) <- object$ltst$nams
            cof
        }
        std.z <- cof[, -1] / sqrt(rowSums(cof[, -1, drop = FALSE] * cof[, -1, drop = FALSE]) + factors)
        if (factors == 1) {
            cof <- cbind(cof, std.z)
            colnames(cof) <- c("(Intercept)", "z1", "std.z1")
            cof
        } else {
            cof <- cbind(cof[, 1:2], std.z[, 1], cof[, 3], std.z[, 2])
            colnames(cof) <- c("(Intercept)", "z1", "std.z1", "z2", "std.z2")
            cof
        }
    } else
        cof
    if (object$IRT.param){
        cofIRT <- IRT.parm(object)$parms
        cof[, "(Intercept)"] <- cofIRT[, 1]
        cof[, "z1"] <- cofIRT[, 2]
        colnames(cof)[1:2] <- c("Dffclt", "Dscrmn")
    }
    if (prob)
        cof <- cbind(cof, "P(x=1|z=0)" = plogis(object$coef[, 1]))
    if (order)
        cof <- cof[order(cof[, 1]), ]
    cof
}
