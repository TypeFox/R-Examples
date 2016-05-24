mult.fig.p <- function(nr.plots, mfrow, mfcol,
         marP = rep(0,4), mgp = c(1.5,.6,0),
         mar = marP + .1 + c(4,4,2,1),
         main = NULL, sub = NULL, adj.sub = 0.5,
         tit.wid = if (is.null(main)) 0 else 1 + 1.5*cex.main,
         quiet = .Device == "postscript",
         cex.main = par("cex.main"),
         col.main = par("col.main"),
         font.main = par("font.main"),
         ...)
{
  ## Purpose: 'MULTiple FIGures' incl. TITLE and other good defaults
  ## -------------------------------------------------------------------------
  ## Arguments: 
  ##             -- Either ONE of the first 3 arguments --
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, 1990 (UW, Seattle) -- 1995
  ## -------------------------------------------------------------------------

  use.row <- missing(mfcol)
  if (use.row)
    if (missing(mfrow)) {
      if (missing(nr.plots))
        stop("must either specify 'nr.plots', 'mfrow' or 'mfcol' !")
      else  mfrow <- n2mfrow (nr.plots)
    }
  oma <- c(tit.wid, 0, tit.wid, 0)
  old.par <-
    if(use.row) par(mfrow = mfrow, oma= oma, mar = mar, mgp= mgp)
    else        par(mfcol = mfcol, oma= oma, mar = mar, mgp= mgp)
  if(!quiet) cat("Execute\n\t par(result$old.par) \n later to restore graphical par\n")
  ##---- now go ahead :
  if(!is.R())
      frame()
  if (!is.null(main)) {# Do title *before* first plot!
      if(is.R()) plot.new()
      mtext(sub, side = 1, outer = TRUE,
            line = 0,
            cex = cex.main/2,
            font = font.main, col = col.main, adj=adj.sub, ...)
      mtext(main, side = 3, outer = TRUE,
            line = cex.main, # was tit.wid - 4,
            cex = cex.main,
            font = font.main, col = col.main, ...)
      if(is.R()) par(new=TRUE)# reverse `plot.new()' above
  }
  invisible(list(new.par = par(c("mfrow","mfcol","oma","mar","mgp")),
                 old.par = old.par))
}
