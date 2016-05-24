# Define plot() function for class "PerFit":
plot.PerFit <- function (x, #x = an object from 'PerFit' class
                         cutoff.obj=NULL, #cutoff.obj = an object from 'PerFit.cutoff' class
                         ModelFit="NonParametric", Nreps=1000, 
                         IP=x$IP, IRT.PModel=x$IRT.PModel, Ability=x$Ability, Ability.PModel=x$Ability.PModel, mu=0, sigma=1, 
                         Blvl = 0.05, Breps = 1000, CIlvl = 0.95, 
                         UDlvl = NA, 
                         # 
                         Type="Density", Both.scale=TRUE, Cutoff=TRUE, Cutoff.int=TRUE, Flagged.ticks = TRUE, 
                         Xlabel=NA, Xcex=1.5, title=NA, Tcex=1.5,
                         col.area="lightpink", col.hist="lightblue", col.int="darkgreen", col.ticks="red", ...)
{  
  # Sanity check - Class PerFit:
  Sanity.cls(x)  
  # 
  upp.PFS  <- c("Cstar", "C.Sato", "U3", "ZU3", "G", "Gnormed", "Gpoly", "Gnormed.poly", "U3poly", "D.KB")
  low.PFS  <- c("r.pbis", "NCI", "Ht", "A.KB", "E.KB", "lz", "lzstar", "lzpoly")
  # 
  pf.scores       <- x$PFscores[[1]]
  PFS.NA          <- is.na(pf.scores)
  pf.scores.noNAs <- pf.scores[!PFS.NA]
  # 
  if (is.null(cutoff.obj))
  {
    cutoff.res <- cutoff(x, ModelFit, Nreps, IP, IRT.PModel, Ability, Ability.PModel, mu, sigma, Blvl, Breps, CIlvl, UDlvl)
  } else
  {
    Sanity.clsPO(cutoff.obj)
    cutoff.res <- cutoff.obj
  }
  # 
  x.line       <- cutoff.res$Cutoff
  perc.flagged <- round(100 * cutoff.res$Prop.flagged, 2)
  direction    <- paste(", ", cutoff.res$Tail, " ", sep = "")
  # 
  PFS.flagged  <- flagged.resp(x, cutoff.res, scores = FALSE)[[1]][,2]
  # Find correct scale for y-axis:
  ymax.hist    <- max(hist(pf.scores.noNAs, plot = FALSE)$density)
  ymax.dens    <- max(density(pf.scores.noNAs)$y)
  ymax         <- switch(Type,
                         Density   = ymax.dens,
                         Histogram = ymax.hist,
                         Both      = if (Both.scale == TRUE) {max(ymax.dens,ymax.hist)} else {min(ymax.dens,ymax.hist)})
  par(mar=c(4,3.5,2,1)+.1,las=1)
  hist(pf.scores.noNAs, freq = FALSE, border = "white", ann = FALSE, ylim = c(0, ymax))
  #
  if (Cutoff == TRUE)
  {
    if (any(x$PFStatistic == upp.PFS))
    {
      rect(x.line, 0, par("usr")[2], par("usr")[4],col = col.area, border = NA)
    }
    if (any(x$PFStatistic == low.PFS))
    {
      rect(par("usr")[1], 0, x.line, par("usr")[4], col = col.area, border = NA)
    }
  }
  #
  if (Type == "Histogram")
  {
    par(new = TRUE)
    hist(pf.scores.noNAs, freq = FALSE, col = col.hist, ann = FALSE, ylim = c(0, ymax))
  }
  #
  if (Type == "Density")
  {
    points(density(pf.scores.noNAs), type = "l", lwd = 2, ann = FALSE, ylim = c(0, ymax))
  }
  #
  if (Type == "Both")
  {
    par(new=TRUE)
    hist(pf.scores.noNAs, freq = FALSE, col = col.hist, ann = FALSE, ylim = c(0,ymax))
    points(density(pf.scores.noNAs), type = "l",lwd = 2,ann = FALSE,ylim = c(0,ymax))
  }
  box(col="black")
  #
  if (Cutoff == FALSE)
  {
    tmp <- if (is.na(Xlabel)) {x$PFStatistic} else {Xlabel}
    mtext(side = 1, text = tmp, line = 2.5, col = "black", cex = 1.5, font = 1)
  }
  # Add flagged respondents:
  if (Flagged.ticks == TRUE)
  {
    axis(3, at = PFS.flagged, labels = FALSE, tick = TRUE, lwd.ticks = 2, col.ticks = col.ticks)
  }
  # Add bootstrap CIlvl% CI to x-axis:
  if (Cutoff.int == TRUE)
  {
    segments(x0 = cutoff.res$Cutoff.CI[1], y0 = par("usr")[3], 
             x1 = cutoff.res$Cutoff.CI[2], y1 = par("usr")[3], lwd = 6, col = col.int, xpd = TRUE)
  }
  # 
  if (Cutoff == TRUE)
  {
    abline(v = x.line, lwd = 2)
    tmp <- if (is.na(Xlabel))
    {
      paste(x$PFStatistic, " (cutoff = ",round(x.line,3),direction,perc.flagged,"%)", sep = "")
    } else 
    {
      Xlabel
    }
    mtext(side = 1, text = tmp, line = 2.5, col = "black", cex = Xcex, font = 1) 
  }
  #
  tmp <- if (is.na(title))
  {
    "Distribution"
  } else 
  {
    title
  }
  mtext(side = 3, text = tmp, line = .5, col = "black", cex = Tcex, font = 2)
}
