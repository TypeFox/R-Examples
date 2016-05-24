## print method
`print.residLen` <- function(x,
                             digits = min(4, getOption("digits") - 4),
                             probs = c(0.5, 0.75, 0.9, 0.95, 0.99), ...) {
    cat("\n")
    writeLines(strwrap("Squared residual lengths",
        prefix = "\t"))
    cat("\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\n")
    writeLines(strwrap(paste("Ordination Method:",
                             attr(x, "method"))))
    cat("\nQuantiles of residual lengths:\n\n")
    quant <- with(x,
                  data.frame(rbind(quantile(train,
                                            probs = probs),
                                   quantile(passive,
                                            probs = probs)
                                   )
                             )
                  )
    names(quant) <- as.character(paste(probs * 100, "%",
                                       sep = ""))
    rownames(quant) <- c("Training Set:", "Passive:")
    print(quant, digits = digits)
    invisible(x)
}
