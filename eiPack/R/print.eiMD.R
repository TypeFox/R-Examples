print.eiMD <- function(x, digits = min(3, getOption("digits") - 5), short = TRUE, ...) {
  cat("\nFormula: ", deparse(x$call$formula), "\n")
  cat("Total sims: ", (x$call$burnin) + (x$call$sample * x$call$thin), 
"\n")
  cat("Burnin discarded: ", x$call$burnin, "\n")
  cat("Sims saved: ", x$call$sample, "\n\n")

  "%w/o%" <- function(x,y) x[!x %in% y]

  if (is.mcmc(x$draws$Cell.counts)) { 
    tnames <- strsplit(colnames(x$draws$Cell.counts), "ccount.")
    get2 <- function(x) x[2]
    idx <- strsplit(sapply(tnames, get2), ".", fixed = TRUE)
    idx <- as.list(as.data.frame(matrix(unlist(idx), byrow = TRUE,
                                        nrow = length(idx), ncol = length(idx[[1]]))))
    idx <- lapply(idx, as.character)
    idx <- lapply(idx, unique)
  }
  else {
    idx <- dimnames(x$draws$Cell.counts)[1:2]
  }
  names(idx) <- c("rows", "columns")
  
  for (ii in names(x$draws) %w/o% c("Beta")) {
    if (is.mcmc(x$draws[[ii]])) {
      nr <- length(idx[[1]])
      nc <- ncol(x$draws[[ii]]) / nr
      cc <- array(t(x$draws[[ii]]), dim = c(nr, nc, nrow(x$draws[[ii]])),
                  dimnames = list(idx[[1]], idx[[2]][1:nc], NULL))
    }
    else {
      cc <- x$draws[[ii]]
      nr <- dim(x$draws[[ii]])[1]
      nc <- dim(x$draws[[ii]])[2] 
      nz <- dim(x$draws[[ii]])[3]
      if (nr == length(idx[[1]]) & is.na(nz)) {
        nc <- 1
      }
    }
    cat(paste("Mean ", ii, ": (averaged over simulations)\n", sep = ""))
    if (nc > 1) 
      print.default(format(apply(cc, c(1,2), mean), digits = digits),
                    print.gap = 2, quote = FALSE)
    else {
      print.default(format(apply(cc, 1, mean), digits = digits),
                    print.gap = 2, quote = FALSE)
    }
    cat("\n")
  }
  invisible(x)
}
  
