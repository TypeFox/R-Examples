### sda.ranking.R  (2013-11-21)
###
###    Shrinkage discriminant analysis (feature ranking)
###
### Copyright 2008-13 Miika Ahdesmaki, Verena Zuber, Sebastian Gibb,
### and Korbinian Strimmer
###
###
### This file is part of the `sda' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


sda.ranking = function(Xtrain, L, lambda, lambda.var, lambda.freqs, 
 ranking.score=c("entropy", "avg", "max"),
 diagonal=FALSE, fdr=TRUE, plot.fdr=FALSE, verbose=TRUE)
{
  ranking.score = match.arg(ranking.score)

  cat = catscore(Xtrain, L,
     lambda=lambda, lambda.var=lambda.var, lambda.freqs=lambda.freqs,
      diagonal=diagonal, verbose=verbose)

  cl.count = dim(cat)[2]
  freqs=attr(cat, "freqs")

  if (ranking.score == "entropy")
    score = as.vector(cat^2 %*% as.matrix(1-freqs))  # weighted sum of squared CAT scores

  if( ranking.score == "avg")
    score = apply(cat^2, 1, sum)/cl.count # average of squared CAT-scores

  if( ranking.score == "max")
    score = apply(cat^2, 1, max)          # max of squared CAT-scores

  names(score) = rownames(cat)
  idx = order(score, decreasing = TRUE)

  if (fdr)
  {
    if (verbose) cat("\nComputing false discovery rates and higher cricitism scores for each feature\n")

    if (cl.count == 2)
    {
      fdr.out = fdrtool(cat[,1], plot=plot.fdr, verbose=FALSE)
    }
    else
    {
      z = score^(1/3) # Wilson-Hilferty transformation to normality

      # center before feeding into fdrtool
      #offset = median(z)
      d = density(z)
      offset = d$x[which.max(d$y)]
      z = z-offset
      fdr.out = fdrtool(z, plot=plot.fdr, verbose=FALSE)
    }
    lfdr = fdr.out$lfdr # local false discovery rates
    pval = fdr.out$pval # p-values

    # compute HC score for each p-value
    HC = hc.score(pval) # function from fdrtool

    ranking = cbind(idx, score[idx], cat[idx, , drop=FALSE], lfdr[idx], HC[idx])
    colnames(ranking) = c("idx", "score", colnames(cat), "lfdr", "HC")
    rm(fdr.out)
  }
  else
  {
    ranking = cbind(idx, score[idx], cat[idx, , drop=FALSE])
    colnames(ranking) = c("idx", "score", colnames(cat))
  }
  rm(cat)

  attr(ranking, "class") = "sda.ranking"
  attr(ranking, "diagonal") = diagonal
  attr(ranking, "cl.count") = cl.count

  return(ranking)
}


plot.sda.ranking = function(x, top=40, arrow.col="blue", zeroaxis.col="red", 
  ylab="Features", main, ...)
{
  if (class(x) != "sda.ranking")
  {
    stop("sda.ranking x needed as input!")
  }

  cl.count = attr(x, "cl.count")
  diagonal = attr(x, "diagonal")

  top = min(nrow(x), top) # just to be sure ...
  if(missing(main)) main = paste("The", top, "Top Ranking Features")

  idx = 2+(1:cl.count)

  cn = colnames(x)[idx]
  rn = rownames(x)[1:top]

  if (diagonal)
  {
    cn = substr(cn, 3, nchar(cn))
    xlab = "t-Scores (Centroid vs. Pooled Mean)"
  }
  else
  {
    cn = substr(cn, 5, nchar(cn))
    xlab = "Correlation-Adjusted t-Scores (Centroid vs. Pooled Mean)"
  }

  if (is.null(rn))
  {
    rn = x[1:top, 1]
  }
  else
  {
    if (any(duplicated(rn)))
    {
      warning("There are duplicated row names! These are converted to unique labels in the dotplot.")
      rn = make.unique(rn)
    }
  }

  ## calculate breaks on x-axis
  xBreaks = pretty(range(x[1:top, idx]))

  ## calculate limits (for each subplot)
  dXlim = sum(diff(xBreaks))
  ylim = c(0, top+1)

  ## save current device parameter
  oldPar = par(no.readonly=TRUE)
  ## and restore them after leaving these function
  on.exit(par(oldPar[c("mar", "mgp", "las", "xaxs", "yaxs")]))

  ## default character width/height in inch
  cin = par("cin") * par("mex")

  ## calculate max rowname width for left margin
  maxRowNameWidth = max(strwidth(as.character(rn), units="inches"))
  ## calculate margin lines (for title(ylab))
  lines = round(maxRowNameWidth/cin[2] + (nchar(ylab) > 0))
  ## left margin: lines+1 (because line counting starts at 0)
  ## + some additional margin
  leftMargin = (lines+1)*cin[2] + 0.05

  ## set device geometry
  par(mai=c(0.7, leftMargin, 0.9, 0.25),
      mgp=c(3, 0.5, 0), las=1, yaxs="i", xaxs="i")

  ## calculate x-limits by combining all subplots
  fullXlim = c(0, dXlim * cl.count)

  ## create empty plot area
  plot(NA, type="n", xaxt="n", yaxt="n", xlim=fullXlim, ylim=ylim,
       main="", xlab="", ylab="", ...)

  ## title
  title(main=main, line=3)

  ## xlab
  title(xlab=xlab, line=2)

  ## ylab
  title(ylab=ylab, line=lines)

  ## calculate x values for border of subplots and x==0 lines
  xlimLeft = dXlim*(0:(cl.count-1))
  xZero = abs(xBreaks[1]) + xlimLeft

  ## border
  abline(v=xlimLeft[-1])

  ## scale y-axis labels if needed
  yusr = diff(par("usr")[3:4])
  sh = sum(strheight(as.character(rn), units="user"))*1.05

  if (yusr < sh)
  {
    cex = yusr/sh
  }
  else
  {
    cex = par("cex.axis")
  }

  ## y-axis
  axis(side=2, at=top:1, labels=rn[1:top], tick=FALSE, cex.axis=cex)

  ## x-axis bottom
  xLabels = head(xBreaks[-1], -1)
  nLabels = length(xLabels)
  xAxisLabels = rep(xLabels, times=cl.count)
  at = xAxisLabels + rep(xZero, each=nLabels)
  above = rep(as.logical(1:cl.count %% 2), each=nLabels)

  ## x-axis above
  axis(side=1, at=at[!above], labels=rep("", sum(!above)))
  axis(side=1, at=at[above], labels=xAxisLabels[above])

  ## x-axis bottom
  axis(side=3, at=at[above], labels=rep("", sum(above)))
  axis(side=3, at=at[!above], labels=xAxisLabels[!above])

  ## title
  mtext(cn, side=3, at=xZero, line=1.5, col=1, font=2)

  ## horizontal gray lines
  abline(h=1:top, col="#808080", lwd=0.25)

  ## vertical gray lines
  abline(v=xZero, col=zeroaxis.col, lty=3, lwd=0.25)

  ## values
  values = x[1:top, idx]
  x0 = rep(xZero, each=top)
  x1 = values + x0
  y0 = y1 = rep(top:1, times=cl.count)
  arrows(x0=x0, y0=y0, x1=x1, y1=y1, length=0, col=arrow.col)
  points(x1, y1, col=arrow.col, pch=19)
}


