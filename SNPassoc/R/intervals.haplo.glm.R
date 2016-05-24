intervals.haplo.glm<-function (o, level = 0.95, sign = 1, FUN = exp, ...) 
{
    if (o$family$family != "binomial") 
        FUN = function(x) x
    z <- abs(qnorm((1 - level)/2))
    co <- summary(o)$coef
    strsplit2 <- function(x, split) {
        ans <- unlist(strsplit(x, split))
        return(ifelse(length(ans) == 2, ans[2], ans))
    }
    control0 <- sapply(dimnames(summary(o)$coef)[[1]], FUN = strsplit2, 
        split = "[.]")
    aux <- grep("(Intercept)", dimnames(summary(o)$coef)[[1]])
    control.geno <- c(1:length(control0))[-aux]
    control <- control0[c(1, control.geno)]
    n.control <- length(control)
    nombres <- rep(NA, n.control)
    freqs <- rep(NA, n.control)
    for (i in 1:n.control) {
        if (control[i] != "rare" & control[i] != "(Intercept)") {
            nombres[i] <- paste(o$haplo.unique[as.numeric(control[i]), 
                ], collapse = "")
            freqs[i] <- o$haplo.freq[as.numeric(control[i])]
        }
        else if (control[i] == "(Intercept)") {
            nombres[i] <- "(Intercept)"
            freqs[i] <- -1
        }
        else {
            nombres[i] <- "rare"
            freqs[i] <- sum(o$haplo.freq[o$haplo.rare])
        }
    }
    or <- FUN(co[, 1] * sign)
    li <- FUN(co[, 1] * sign - z * co[, 2])
    ls <- FUN(co[, 1] * sign + z * co[, 2])
    if (o$family$family != "binomial") 
        or <- c(or[1], or[1], or[-1])
    else or <- c(or[1], 1, or[-1])
    li <- c(li[1], NA, li[-1])
    ls <- c(ls[1], NA, ls[-1])
    pvals <- co[, 4]
    pvals <- c(pvals[1], NA, pvals[-1])
    nombre.ref <- paste(o$haplo.unique[o$haplo.base, ], collapse = "")
    nombre.cov <- dimnames(summary(o)$coef)[[1]][-c(1:n.control)]
    nombres <- c(nombres[1], nombre.ref, nombres[-1], nombre.cov)

    n.control.2<-length(o$haplo.freq)-length(o$haplo.rare)
    if (n.control.2!=length(o$coef))
     {
      nombres[(n.control.2+2):length(nombres)]<-names(o$coef)[-c(1:n.control.2)]
     }



    ncov <- length(nombre.cov)
    freqs <- c(freqs[1], o$haplo.freq[o$haplo.base], freqs[-1], 
        rep(NA, ncov))
    names(freqs) <- names(or)
    r <- cbind(freqs, or, li, ls, pvals)
    if (o$family$family != "binomial") 
        dimnames(r) <- list(nombres, c("freq", "diff", paste(level * 
            100, "%", sep = ""), "C.I.", "P-val"))
    else dimnames(r) <- list(nombres, c("freq", "or", paste(level * 
        100, "%", sep = ""), "C.I.", "P-val"))
    class(r) <- "intervals"
    r
}
