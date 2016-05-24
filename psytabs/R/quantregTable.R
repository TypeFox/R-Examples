quantregTable <-
  function (x, digits = 2, significance="none") {
    taus <- function(x) x$tau
    taus <- unlist(lapply(x, taus))
    taus <- format(round(taus, digits))
    coef <- lapply(x, stats::coefficients)
    p <- nrow(coef[[1]])
    k <- ncol(coef[[1]])
    m <- length(taus)
    rlab <- dimnames(coef[[1]])[[1]]
    clab <- taus
    a <- array(unlist(coef), c(p, k, m))
    table <- matrix("", p, m)
    for (i in 1:m) {
      for (j in 1:p) {
        if (k == 3) {
          table[j, i] <- paste(round(a[j, 1, i], digits), " (", round(a[j, 2, i], digits), ", ", round(a[j, 3, i], digits), ")", sep = "")
        }
        else if (k == 4) {
          table[j, i] <- paste(if(significance=="bold"){ifelse(a[j, 4, i] < .05, "\\b ", "")}, round(a[j, 1, i], digits), if(significance=="bold"){ifelse(a[j, 4, i] < .05, " \\b0", "")}, if(significance=="stars"){ifelse(a[j, 4, i] < .001, "***", ifelse(a[j, 4, i] < .01, "**", ifelse(a[j, 4, i] < .05, "* ", "")))}, " (", round(a[j, 2, i], digits), ")", "", sep = "")
        }
      }
    }
    rowlabel <- "Covariates"
    dimnames(table) <- list(rlab, clab)
    table
  }

