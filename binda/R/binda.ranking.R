### binda.ranking.R  (2015-02-28)
###
###    Multi-class discriminant analysis with binary predictors
###
### Copyright 2013-2015  Sebastian Gibb and Korbinian Strimmer
###
###
### This file is part of the `binda' library for R and related languages.
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




binda.ranking = function(Xtrain, L, lambda.freqs, verbose=TRUE)
{
  if (!is.binaryMatrix(Xtrain)) stop("Training data must be a matrix containing only 0 or 1!")
  if (missing(L)) stop("Class labels are missing!")
  L = factor(L) # make sure L is a factor

  if (verbose) reportDetails(Xtrain, L)

  # class frequencies
  freqs = getClassFreqs(L, lambda.freqs=lambda.freqs, verbose=verbose)
  attr(freqs, "lambda.freqs") = NULL

  # means (columns=classes)
  mu = getClassMeans(Xtrain, L)

  # (approximate) mutual information between each feature and the response
  score = rankingScore(mu, freqs, tscores=TRUE, n=length(L))

  idx = order(score[1,], decreasing = TRUE)

  ranking = cbind(idx, t(score[,idx]))
  colnames(ranking) = c("idx", rownames(score))
  rownames(ranking) = colnames(Xtrain)[idx]

  attr(ranking, "class") = "binda.ranking"
  attr(ranking, "cl.count") = dim(score)[1]-1

  return(ranking)
}

# modified from sda package
plot.binda.ranking = function(x, top=40, arrow.col="blue", zeroaxis.col="red",
  ylab="Variables", main, ...)
{
  if (class(x) != "binda.ranking")
  {
    stop("binda.ranking x needed as input!")
  }

  cl.count = attr(x, "cl.count")

  top = min(nrow(x), top) # just to be sure ...
  if(missing(main)) main = paste("The", top, "Top Ranking Variables")

  idx = 2+(1:cl.count)

  cn = colnames(x)[idx]
  rn = rownames(x)[1:top]

  cn = substr(cn, 3, nchar(cn))
  xlab = "t-Scores (Class Mean vs. Pooled Mean)"


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




## private functions
rankingScore = function(mu, freqs, tscores=FALSE, n=1)
{
  # relative entropy ( mutual information between a binary feature and the response)
  # approximation (LR chi squared statistic)
  mu = avoidBoundaries(mu)
  mu.pooled = as.vector(mu %*% freqs)

  k1 = (mu-mu.pooled)^2/( mu.pooled*(1-mu.pooled) )
  t2 = t(t(k1) * freqs/(1-freqs))
  s = colSums( t(t2) * 1/2*(1-freqs) )
  s = 2*n*s # adjust scale with sample size (so that s is the sum of squared t)
  names(s) = "score"

  if(tscores) # also return the group-specific t-scores
  {
    t = sqrt(t2)*sign(mu-mu.pooled)
    t = sqrt(n)*t # adjust scale with sample size to match conventional t-score
    s = t(cbind(s, t))
    rownames(s) = c("score", paste("t", colnames(mu), sep="."))
  }

  return(s)
}

