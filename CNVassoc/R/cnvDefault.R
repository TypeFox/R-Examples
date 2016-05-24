cnvDefault  <-
function (x, num.copies, num.class, cnv.tol = 0.001, mix.method = "mixdist", check.probs = TRUE,
    threshold.0, threshold.k, mu.ini, sigma.ini, pi.ini, cutoffs = NULL, check.alpha = 0.05, check.cnv = TRUE, var.equal)
{
    if (any(is.na(x)))
        stop("missing data in intensities or probabilities is not allowed")
    if (missing(num.class))
        num.class <- 2:6
    if (!missing(num.copies))
        num.class <- length(num.copies)
    copynum.range <- num.class
    if (!missing(num.copies))
        num.class = length(num.copies)
    if (!is.vector(x) & !is.matrix(x) & !is.data.frame(x))
        stop("x must be either a vector or a matrix or a data.frame")
    if (NCOL(x) == 1) {
        x <- as.vector(x)
        if (is.null(cutoffs)){
          out <- mixture(x, num.class, mix.method, threshold.0, threshold.k, mu.ini, sigma.ini, pi.ini, var.equal)
          if (missing(num.copies))
            num.copies <- attr(out, "num.copies")
        } else{
          out <- as.integer(cut(x,c(-Inf,cutoffs,Inf)))
          k <- length(cutoffs) + 1
          if (missing(num.copies))
            num.copies <- 0:(k-1)
          attr(out, "probabilities") <- sapply(1:k, function(j) ifelse(out == j, 1,0 ))
        }
    } else {
        if (NCOL(x) > 1) {
            rs <- apply(x, 1, sum)
            if (check.probs)
              if (any(!sapply(rs, function(temp) isTRUE(all.equal(temp, 1, tolerance = cnv.tol)))))
                stop("Rows of x (probabilities) must sum 1")
            k <- NCOL(x)
            out <- apply(x, 1, which.max)
            attr(out, "k") <- k
            attr(out, "probabilities") <- x
            if (missing(num.copies))
                num.copies <- 1:k - 1
            if (k != length(num.copies))
                stop("number of columns of x must be equal to length of num.copies")
        }
        else {
            stop("x must have positive NCOL")
        }
    }
    attr(out, "k") <- ncol(attr(out,"probabilities"))
    attr(out, "num.copies") <- num.copies
    attr(out, "pi") <- colMeans(attr(out, "probabilities"))
    attr(out, "copynum.range") <- copynum.range
    attr.old <- attributes(out)
    out <- attr(out, "num.copies")[out]
    attributes(out) <- attr.old
    class(out) <- "cnv"
    if (check.cnv)
      if (!is.null(attr(out,"meanRatio")))
        if (!checkCNV(out,check.alpha)){
          plot(out)
          stop("CNV intensities are not fitted properly\ntry\n   1-. cnv defined by threshold, use locator (see 'plotSignal')\n   2-. change initial values (mu.ini, sigma.ini or pi.ini)\n   3-. change mix.method\n   4-. change number of classes (num.class)\n   ...")
        }
    out
}
