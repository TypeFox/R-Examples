'summary.DEoptim' <- function(object, ...){
  digits <- max(5, getOption('digits') - 2)

  cat("\n***** summary of DEoptim object *****",
      "\nbest member   : ", round(object$optim$bestmem, digits),
      "\nbest value    : ", round(object$optim$bestval, digits),
      "\nafter         : ", round(object$optim$iter), "generations",
      "\nfn evaluated  : ", round(object$optim$nfeval), "times",
      "\n*************************************\n")

  invisible(object)
}

plot.DEoptim <- function (x, plot.type = c("bestmemit", "bestvalit",
"storepop"), ...) {
    z <- x$member
    niter <- length(z$bestvalit)
    npar <- length(z$lower)
    nam <- names(z$lower)
    nstorepop <- length(z$storepop)
    if (identical(plot.type[1], "bestmemit")) {
        plot.new()
        par(mfrow = c(min(3, npar), 1))
        for (i in 1:npar) {
            if (identical(i%%4, 0)) {
                cat("new plot\n")
                devAskNewPage(ask = TRUE)
            }
            plot(1:niter, z$bestmemit[, i], xlim = c(1, niter),
                las = 1, xlab = "iteration", ylab = "value",
                main = nam[i], ...)
            abline(h = c(z$lower[i], z$upper[i]), col = "red")
        }
    }
    else if (identical(plot.type[1], "bestvalit")) {
        plot(1:niter, z$bestvalit, xlim = c(1, niter), las = 1,
            xlab = "iteration", ylab = "function value", main =
"convergence plot",
            ...)
    }
    else if (identical(plot.type[1], "storepop") && nstorepop > 0) {
        plot.new()
        par(mfrow = c(min(3, npar), 1))
        for (i in 1:npar) {
            if (identical(i%%4, 0)) {
                cat("new plot\n")
                devAskNewPage(ask = TRUE)
            }
            tmp <- NULL
            for (j in 1:nstorepop) {
                tmp <- cbind(tmp, z$storepop[[j]][, i])
            }
            matplot(t(tmp), col = "black", pch = 20, las = 1, xlab =
"stored population", ylab = "value",
                main = nam[i], ...)
            abline(h = c(z$lower[i], z$upper[i]), col = "red")
            par(new = FALSE)
        }
    }
    else {
        warning("'plot.type' does not correspond to any plotting type",
immediate. = TRUE)
    }
}
