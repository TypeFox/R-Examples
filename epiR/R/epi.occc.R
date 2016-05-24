## x is a matrix line object, rows are cases, columns are raters
## na.rm: logical, if NAs should be excluded
## pairs: logical, if pairwise statistic values should be returned as 
##        part of the return value

epi.occc <- function(dat, na.rm = FALSE, pairs = FALSE){
  
  ## Create a list to hold all variables:
  elements <- list()
  
  ## Do all data manipulation within the list:
  elements <- within(elements, {
    
  if (!na.rm) {
        m <- apply(dat, 2, mean)
        s <- apply(dat, 2, sd)
        COV <- cov(dat)
    } else {
        m <- apply(dat, 2, mean, na.rm = TRUE)
        s <- apply(dat, 2, sd, na.rm = TRUE)
        COV <- cov(dat, use = "pairwise.complete.obs")
    }
   
    J <- ncol(dat)
    j <- col(matrix(0,J,J))[lower.tri(matrix(0,J,J))]
    k <- row(matrix(0,J,J))[lower.tri(matrix(0,J,J))]
    n <- (J * J - J) / 2
    v <- numeric(n)
    u <- numeric(n)
    ksi <- numeric(n)
    ccc <- numeric(n)
    for (i in seq_len(n)) {
        v[i] <- s[j[i]] / s[k[i]]
        u[i] <- (m[j[i]] - m[k[i]]) / sqrt(s[j[i]] * s[k[i]])
        ksi[i] <- s[j[i]]^2 + s[k[i]]^2 + (m[j[i]] - m[k[i]])^2
        ccc[i] <- (2 * COV[j[i], k[i]]) / ksi[i]
    }
    
    accu <- ((v + 1/v + u^2) / 2)^-1
    prec <- ccc / accu
    occc <- sum(ksi * ccc) / sum(ksi)
    oaccu <- sum(ksi * accu) / sum(ksi)
    oprec <- occc / oaccu
    prs <- if (pairs) {
        list(ccc = ccc, prec = prec, accu = accu, ksi = ksi, scale = v, location = u)
    } else NULL
  })  
  
    rval <- list(occc = elements$occc, oprec = elements$oprec, oaccu = elements$oaccu, pairs = elements$prs, data.name = deparse(substitute(dat)))
    class(rval) <- "epi.occc"
    return(rval)
}

# https://cran.r-project.org/web/packages/knitr/vignettes/knit_print.html
print.epi.occc <- function(x, ...) {
    # cat("Overall concordance correlation coefficients\n")
    cat(sprintf("\nOverall CCC           %.4f", x$occc))
    cat(sprintf("\nOverall precision     %.4f", x$oprec))
    cat(sprintf("\nOverall accuracy      %.4f", x$oaccu))
    cat("\n")
    # print(data.frame(Value = c("Overall CCC" = x$occc, "Overall precision" = x$oprec, "Overall accuracy" = x$oaccu)), ...)
}

## Summary method for epi.occc:
summary.epi.occc <- function(object, ...) {
   out <- data.frame(occc = object$occc, oprec = object$oprec, oaccu = object$oaccu) 
   return(out)
}