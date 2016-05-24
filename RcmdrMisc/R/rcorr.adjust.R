# the following function is adapted from a suggestion by Robert Muenchen

# uses rcorr in the Hmisc package

# last modified 2014-09-04 by J. Fox

rcorr.adjust <- function (x, type = c("pearson", "spearman"), 
            use = c("complete.obs", "pairwise.complete.obs")) {
    opt <- options(scipen = 5)
    on.exit(options(opt))
    type <- match.arg(type)
    use <- match.arg(use)
    x <- if (use == "complete.obs") 
      as.matrix(na.omit(x))
    else as.matrix(x)
    R <- rcorr(x, type = type)
    P <- P.unadj <- R$P
    p <- P[lower.tri(P)]
    adj.p <- p.adjust(p, method = "holm")
    P[lower.tri(P)] <- adj.p
    P[upper.tri(P)] <- 0
    P <- P + t(P)
    P <- ifelse(P < 1e-04, 0, P)
    P <- format(round(P, 4))
    diag(P) <- ""
    P[c(grep("0.0000", P), grep("^ 0$", P))] <- "<.0001"
    P[grep("0.000$", P)] <- "<.001"
    P.unadj <- ifelse(P.unadj < 1e-04, 0, P.unadj)
    P.unadj <- format(round(P.unadj, 4))
    diag(P.unadj) <- ""
    P.unadj[c(grep("0.0000$", P.unadj), grep("^ 0$", P.unadj))] <- "<.0001"
    P.unadj[grep("0.000$", P.unadj)] <- "<.001"
    result <- list(R = R, P = P, P.unadj = P.unadj, type = type)
    class(result) <- "rcorr.adjust"
    result
  }

print.rcorr.adjust <- function(x, ...){
    cat("\n", if (x$type == "pearson") "Pearson" else "Spearman", "correlations:\n")
    print(round(x$R$r, 4))
    cat("\n Number of observations: ")
    n <- x$R$n
    if (all(n[1] == n)) cat(n[1], "\n")
    else{
        cat("\n")
        print(n)
    }
    cat("\n Pairwise two-sided p-values:\n")
    print(x$P.unadj, quote=FALSE)
    cat("\n Adjusted p-values (Holm's method)\n")
    print(x$P, quote=FALSE)
}
