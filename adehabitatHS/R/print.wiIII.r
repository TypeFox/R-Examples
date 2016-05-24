"print.wiIII" <- function(x, ...)
{
    if (!inherits(x,"wiIII"))
        stop("x should be of class \"wiIII\"")
    cat("\n\n************** Manly's Selection ratios for design III ********\n\n")
    cat("1. Test of habitat selection for each animal:\n\n")
    print(x$Khi2Lj)

    cat("\n\n2. Test of overall habitat selection:\n")
    print(x$Khi2L)
    cat("\n\nTable of selection ratios:\n")
    print(data.frame(Wi=x$wi,
                     SE=x$se.wi, IClower=x$ICwilower,
                     ICupper=x$ICwiupper), ...)
    cat("\n\nBonferroni classement \nBased on", (1 - x$alpha) *
        100, "% confidence intervals on the differences of Wi :\n")
    print(x$profile, quote = FALSE)
    cat("\n")
}

