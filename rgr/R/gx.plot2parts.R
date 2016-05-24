gx.plot2parts <-
function (xx1, xx2, x1lab = deparse(substitute(xx1)), x2lab = deparse(substitute(xx2)), 
     cex = 0.8, ifwarn = TRUE, ...) 
{
     # Function to support the investigation of 2 parts of a composition following 
     # the procedures outlined in Filzmoser et al. (2010).  The stability measure
     # (replacing the correlation coefficient) measuring the consistency of the
     # ratio of part xx1 to part xx2, is estimated robustly with the MAD.
     #
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variables, xx1 and xx2, must be run through ltdl.fix.df to convert any <dl
     # -ve values to positive half that value, and set zero2na = TRUE if it is
     # required to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     frame()
     oldpar <- par()
     on.exit(par(oldpar))
     par(mfrow = c(2, 2), cex.main = 0.9)
     if (length(xx1) != length(xx2)) 
         stop("Input vectors must be of equal length\n")
     if(ifwarn) cat("  ** Are the data in the same measurement units? **\n" )
     temp.x <- remove.na(cbind(xx1, xx2))
     x1 <- temp.x$x[1:temp.x$n, 1]
     x2 <- temp.x$x[1:temp.x$n, 2]
     plot(x1, x2, xlab = x1lab, ylab = x2lab, log = "xy", ...)
     par(pty = "m")
     x1x2.ilr <- log(x1/x2)/1.4142
     ilr.mad <- mad(x1x2.ilr)
     ilr.stab <- exp(-ilr.mad * ilr.mad)
     ilr.name <- paste("ilr(", deparse(substitute(xx1)), ",", 
         deparse(substitute(xx2)), ")", sep = "")
     ratio.name <- paste("Log10(", deparse(substitute(xx1)), "/", 
         deparse(substitute(xx2)), ")", sep = "")
     bxplot(x1x2.ilr, xlab = ilr.name, main = "Tukey Boxplot", cex = cex, ...)
     limits <- par("usr")
     xpos <- limits[1] + (limits[2] - limits[1]) * 0.05
     ypos <- limits[4] - (limits[4] - limits[3]) * 0.15
     text(xpos, ypos, labels = paste("Robust ilr stability =", signif(ilr.stab, 3)),
         adj = 0, cex = cex)
     plot(seq(1:temp.x$n), log10(x1/x2), ylim = c(-6, 6), xlab = "Index", 
         ylab = ratio.name, main = "Index Plot", ...)
     gx.ecdf(x1x2.ilr, xlab = ilr.name, ylab = " ", ifqs = T,
         main = "Empirical Cumulative Distribution\nFunction (ECDF)", cex = cex, ...)
     invisible()
}
