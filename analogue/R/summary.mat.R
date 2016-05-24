###########################################################################
##                                                                       ##
## mat - function to perform the modern analogue technique for           ##
##       environmental reconstruction                                    ##
##                                                                       ##
## Created       : 17-Apr-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.0-1                                                 ##
## Last modified : 17-Apr-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
##                                                                       ##
###########################################################################
summary.mat <- function(object, k = 10,
                        digits = min(2, getOption("digits") - 4),
                        ...)
  {
    tbl <- cbind(object$standard$rmsep[1:k], object$standard$r.squared[1:k],
                 object$standard$avg.bias[1:k], object$standard$max.bias[1:k])
    tbl.w <- cbind(object$weighted$rmsep[1:k], object$weighted$r.squared[1:k],
                   object$weighted$avg.bias[1:k], object$weighted$max.bias[1:k])
    tbl <- as.matrix(format(tbl, digits = digits))
    tbl.w <- as.matrix(format(tbl.w, digits = digits))
    tbl <- cbind(as.integer(1:k), tbl)
    tbl.w <- cbind(as.integer(1:k), tbl.w)
    rownames(tbl) <- rownames(tbl.w) <- rep("", nrow(tbl))
    colnames(tbl) <- colnames(tbl.w) <- c("k",
                                          "RMSEP","R2","Avg Bias","Max Bias")
    W.Est <- object$weighted$est[k, ]
    Est <- object$standard$est[k, ]
    Obs <- object$orig.y
    W.Resi <- resid(object, k = k, weighted = TRUE)$residuals#[, k]
    Resi <- resid(object, k = k)$residuals#[, k]
    minDC <- apply(object$Dij, 2, minDij)
    minW.Resi <- apply(object$weighted$resid[1:k, , drop = FALSE], 2,
                       function(x) {min(abs(x))})
    minResi <- apply(object$standard$resid[1:k, , drop = FALSE], 2,
                       function(x) {min(abs(x))})
    W.n.closest <- apply(object$weighted$resid[1:k, , drop = FALSE], 2,
                         function(x) {which.min(abs(x))})
    n.closest <- apply(object$standard$resid[1:k, , drop = FALSE], 2,
                       function(x) {which.min(abs(x))})
    structure(list(summ = data.frame(Obs = Obs, Est = Est, Resi = Resi,
                     W.Est = W.Est, W.Resi = W.Resi,
                     minDC = minDC,
                     minResi = minResi, k = n.closest,
                     minW.Resi = minW.Resi, k.W = W.n.closest),
                   tbl = tbl, tbl.w = tbl.w, call = object$call,
                   quantiles = quantile(object$Dij[lower.tri(object$Dij)],
                     probs = c(0.01, 0.02, 0.05, 0.1, 0.2))),
              class = "summary.mat",
              k = k)
  }

print.summary.mat <- function(x,
                              digits = min(3, getOption("digits") - 4),
                              ...)
  {
    cat("\n")
    writeLines(strwrap("Modern Analogue Technique", prefix = "\t"))
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    cat("\nPercentiles of the dissimilarities for the training set:\n\n")
    print(x$quantiles, digits = digits)
    cat("\nInferences based on the mean of k-closest analogues:\n\n")
    print(x$tbl, quote = FALSE, right = TRUE)
    cat("\nInferences based on the weighted mean of k-closest analogues:\n\n")
    print(x$tbl.w, quote = FALSE, right = TRUE)
    k <- attr(x, "k")
    x <- x$summ
    class(x) <- "data.frame"
    cat("\nResults for training set\n")
    cat(paste("\n  * (W.)Est and (W.)Resi are based on k=",
              k, "-closest analogues", sep = ""))
    cat("\n  * minDC is the minimum distance to another sample in the training set")
    cat("\n  * min(W.)Resi is the minimum residual for a k-closest model,")
    cat(paste("\n    where k = 1,...,", k,
              ". Column k(.W) displays which k has minResi\n\n", sep = ""))
    print(x, digits = digits, print.gap = 2)
    cat("\n")
    invisible(x)
  }
