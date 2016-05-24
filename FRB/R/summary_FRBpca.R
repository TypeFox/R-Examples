
summary.FRBpca <- function(object, confmethod = c("BCA","basic","both"), digits=3, ...) {

confmethod <- match.arg(confmethod)
q <- length(object$eigval)

if (!is.null(dimnames(object$est$Mu)[[2]]))
    dimnames(object$eigvec) <- list(dimnames(object$est$Mu)[[2]], paste("PC",1:q,sep=""))
else
    dimnames(object$eigvec) <- list(paste("V",1:q,sep=""), paste("PC",1:q,sep=""))

eigvalstoprint <- t(as.matrix(object$eigval))
if (confmethod == "BCA") {
  eigvalstoprint <- rbind(eigvalstoprint,t(object$eigval.CI.bca))
  dimnames(eigvalstoprint) <- list(c("      estimates", paste("  BCa ", object$conf*100, "% lower",sep=""),paste("  BCa ", object$conf*100, "% upper",sep="")), paste("PC", 1:q, sep=""))
}
else if (confmethod == "basic") {
  eigvalstoprint <- rbind(eigvalstoprint,t(object$eigval.CI.basic))
  dimnames(eigvalstoprint) <- list(c("      estimates", paste("basic ", object$conf*100, "% lower",sep=""),paste("basic ", object$conf*100, "% upper",sep="")), paste("PC", 1:q, sep=""))
}
else if (confmethod == "both") {
  eigvalstoprint <- rbind(eigvalstoprint,t(object$eigval.CI.bca),t(object$eigval.CI.basic))
  dimnames(eigvalstoprint) <- list(c("      estimates", paste("  BCa ", object$conf*100, "% lower",sep=""),paste("  BCa ", object$conf*100, "% upper",sep=""),paste("basic ", object$conf*100, "% lower",sep=""),paste("basic ", object$conf*100, "% upper",sep="")), paste("PC", 1:q, sep=""))
}
avgangletoprint <- t(as.matrix(object$avgangle))
dimnames(avgangletoprint) <- list("", paste("PC", 1:q, sep=""))

if (confmethod == "BCA") {
  pvartable <- cbind(object$pvar*100, object$pvar.CI.bca*100)
  pvartable <- rbind(pvartable, rep(100,3))
  dimnames(pvartable) <- list(paste("k=", 1:q, sep=""),c("Est.", paste("(BCa ",object$conf*100,"% lower",sep=""), "  upper)"))
}
else if (confmethod == "basic") {
  pvartable <- cbind(object$pvar*100, object$pvar.CI.basic*100)
  pvartable <- rbind(pvartable, rep(100,3))
  dimnames(pvartable) <- list(paste("k=", 1:q, sep=""),c("Est.", paste("(basic ",object$conf*100,"% lower",sep=""), "  upper)"))
}
else if (confmethod == "both") {
  pvartable <- cbind(object$pvar*100, object$pvar.CI.bca*100, object$pvar.CI.basic*100)
  pvartable <- rbind(pvartable, rep(100,5))
  dimnames(pvartable) <- list(paste("k=", 1:q, sep=""),c("Est.", paste("(BCa ",object$conf*100,"% lower",sep=""), "  upper)", paste("(basic ",object$conf*100,"% lower",sep=""), "  upper)"))
}

res <- list(eigvals=eigvalstoprint, eigvecs=object$eigvec, avgangle=avgangletoprint, pvars=pvartable, method=object$method, digits=digits)
class(res) <- "summary.FRBpca"

res

}

###############################################################################

print.summary.FRBpca <- function(x, ...) {

cat(x$method, "\n\n")
cat("Eigenvalues:\n")
print(x$eigvals, digits = x$digits)
cat("\n")
cat("Principal components loadings:\n")
print(x$eigvecs, digits = x$digits)
cat("\n")
cat("Average angle between PC and its bootstrapped versions:\n")
print(x$avgangle, digits = x$digits)
cat("(in [0 - pi/2], cf. aligned - perpendicular)\n")
cat("\n")
cat("Percentage of variance explained by first k components:\n")
print(x$pvars, digits = x$digits)

}

###############################################################################

print.FRBpca <- function(x, digits=3, ...) {

q <- length(x$eigval)

eigvalstoprint <- t(as.matrix(x$eigval))
eigvalstoprint <- rbind(eigvalstoprint,t(x$eigval.CI.bca))
dimnames(eigvalstoprint) <- list(c("      estimates", paste("  BCa ", x$conf*100, "% lower",sep=""),paste("  BCa ", x$conf*100, "% upper",sep="")), paste("PC", 1:q, sep=""))

cat(x$method, "\n\n")
cat("Standard deviations:\n")
print(sqrt(pmax(eigvalstoprint,0)), digits = digits)
cat("\n")

}