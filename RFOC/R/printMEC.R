`printMEC` <-
function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("Plane 1: ")
    print.default(format(c(x$F$az, x$F$dip), digits = digits), print.gap = 2,
            quote = FALSE)
    cat("Plane 2: ")
    print.default(format(c(x$G$az, x$G$dip), digits = digits), print.gap = 2,
            quote = FALSE)

    cat("Vector 1: ")
    print.default(format(c(x$U$az, x$U$dip), digits = digits), print.gap = 2,
            quote = FALSE)
    cat("Vector 2: ")
    print.default(format(c(x$V$az, x$V$dip), digits = digits), print.gap = 2,
            quote = FALSE)

    cat("P-axis: ")
    print.default(format(c(x$P$az, x$P$dip), digits = digits), print.gap = 2,
            quote = FALSE)
    cat("T-axis: ")
    print.default(format(c(x$T$az, x$T$dip), digits = digits), print.gap = 2,
            quote = FALSE)

    

    cat("\n")


  }

