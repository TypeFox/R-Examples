"print.wiII" <- function(x, ...)
{
    if (!inherits(x,"wiII"))
        stop("x should be of class \"wiII\"")
    cat("\n\n************** Manly's Selection ratios for design II ********\n\n")
    cat("1. Test of identical use of habitat by all animals\n")
    cat("   (Classical Khi-2 performed on the used matrix):\n")
    print(x$Khi2L1)

    cat("2. Test of overall habitat selection:\n")
    print(x$Khi2L2)
    cat("3. Test of hypothesis that animals are on average using resources\n")
    cat("   in proportion to availability, irrespective of whether they are\n")
    cat("   the same or not (Khi2L2 - Khi2L1):\n")
    print(x$Khi2L2MinusL1)
    cat("\n\nTable of selection ratios:\n")
    print(data.frame(Available=x$avail.prop, Used=x$used.prop, Wi=x$wi,
                     SE=x$se.wi, IClower=x$ICwilower, ICupper=x$ICwiupper),
          ...)
    cat("\n\nBonferroni classement \nBased on", (1 - x$alpha) *
        100, "% confidence intervals on the differences of Wi :\n")
    print(x$profile, quote = FALSE)
    cat("\n")
}

