
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2011-11-23 11:11:32 +0100 (Wed, 23 Nov 2011) $
# $Rev: 183 $
# ----------------------------------------------------------------

plot.copyright <- function(){
  ## Prints copyright message to plot device

  this.year <- strftime(Sys.time(), format = "%Y")
  copyr.message <- paste("Copyright \u00A9", this.year,
                         "by Wegener Center")
  mtext(copyr.message, side=1, outer = TRUE, adj = 1, cex = 0.5)

  invisible()
}
