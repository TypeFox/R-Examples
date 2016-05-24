#' Functions for positioning tick labels on axes
#'
#' \tabular{ll}{
#' Package: \tab labeling\cr
#' Type: \tab Package\cr
#' Version: \tab 0.2\cr
#' Date: \tab 2011-04-01\cr
#' License: \tab Unlimited\cr
#' LazyLoad: \tab yes\cr
#' }
#' 
#' Implements a number of axis labeling schemes, including those 
#' compared in An Extension of Wilkinson's Algorithm for Positioning Tick Labels on Axes
#' by Talbot, Lin, and Hanrahan, InfoVis 2010.
#'
#' @name labeling-package
#' @aliases labeling
#' @docType package
#' @title Axis labeling
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @references
#' Heckbert, P. S. (1990) Nice numbers for graph labels, Graphics Gems I, Academic Press Professional, Inc.
#' Wilkinson, L. (2005) The Grammar of Graphics, Springer-Verlag New York, Inc.
#' Talbot, J., Lin, S., Hanrahan, P. (2010) An Extension of Wilkinson's Algorithm for Positioning Tick Labels on Axes, InfoVis 2010.
#' @keywords dplot
#' @seealso \code{\link{extended}}, \code{\link{wilkinson}}, \code{\link{heckbert}}, \code{\link{rpretty}}, \code{\link{gnuplot}}, \code{\link{matplotlib}}, \code{\link{nelder}}, \code{\link{sparks}}, \code{\link{thayer}}, \code{\link{pretty}}
#' @examples
#' heckbert(8.1, 14.1, 4)	# 5 10 15
#' wilkinson(8.1, 14.1, 4)	# 8 9 10 11 12 13 14 15
#' extended(8.1, 14.1, 4)	# 8 10 12 14

#' # When plotting, extend the plot range to include the labeling
#' # Should probably have a helper function to make this easier
#' data(iris)
#' x <- iris$Sepal.Width
#' y <- iris$Sepal.Length
#' xl <- extended(min(x), max(x), 6)
#' yl <- extended(min(y), max(y), 6)
#' plot(x, y, 
#'     xlim=c(min(x,xl),max(x,xl)), 
#'     ylim=c(min(y,yl),max(y,yl)), 
#'     axes=FALSE, main="Extended labeling")
#' axis(1, at=xl)
#' axis(2, at=yl)
c()



#' Heckbert's labeling algorithm
#'
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @return vector of axis label locations
#' @references
#' Heckbert, P. S. (1990) Nice numbers for graph labels, Graphics Gems I, Academic Press Professional, Inc.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
heckbert <- function(dmin, dmax, m)
{
    range <- .heckbert.nicenum((dmax-dmin), FALSE)
    lstep <- .heckbert.nicenum(range/(m-1), TRUE)
    lmin <- floor(dmin/lstep)*lstep
    lmax <- ceiling(dmax/lstep)*lstep
    seq(lmin, lmax, by=lstep)
}
    
.heckbert.nicenum <- function(x, round)
{
	e <- floor(log10(x))
	f <- x / (10^e)
	if(round)
	{
		if(f < 1.5) nf <- 1
		else if(f < 3) nf <- 2
		else if(f < 7) nf <- 5
		else nf <- 10
	}
	else
	{
		if(f <= 1) nf <- 1
		else if(f <= 2) nf <- 2
		else if(f <= 5) nf <- 5
		else nf <- 10
	}
	nf * (10^e)
}



#' Wilkinson's labeling algorithm
#'
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @param Q set of nice numbers
#' @param mincoverage minimum ratio between the the data range and the labeling range, controlling the whitespace around the labeling (default = 0.8)
#' @param mrange range of \code{m}, the number of tick marks, that should be considered in the optimization search
#' @return vector of axis label locations
#' @note Ported from Wilkinson's Java implementation with some changes.
#'	Changes: 	1) m (the target number of ticks) is hard coded in Wilkinson's implementation as 5. 
#'						Here we allow it to vary as a parameter. Since m is fixed, 
#'						Wilkinson only searches over a fixed range 4-13 of possible resulting ticks.
#'						We broadened the search range to max(floor(m/2),2) to ceiling(6*m), 
#'						which is a larger range than Wilkinson considers for 5 and allows us to vary m,
#'						including using non-integer values of m.
#'				2) Wilkinson's implementation assumes that the scores are non-negative. But, his revised
#'						granularity function can be extremely negative. We tweaked the code to allow negative scores.
#'						We found that this produced better labelings.
#'				3) We added 10 to Q. This seemed to be necessary to get steps of size 1.
#'	It is possible for this algorithm to find no solution.
#'					In Wilkinson's implementation, instead of failing, he returns the non-nice labels spaced evenly from min to max.
#'					We want to detect this case, so we return NULL. If this happens, the search range, mrange, needs to be increased.
#' @references
#' Wilkinson, L. (2005) The Grammar of Graphics, Springer-Verlag New York, Inc.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
wilkinson <-function(dmin, dmax, m, Q = c(1,5,2,2.5,3,4,1.5,7,6,8,9), mincoverage = 0.8, mrange=max(floor(m/2),2):ceiling(6*m))
{
	best <- NULL
	for(k in mrange)
	{
		result <- .wilkinson.nice.scale(dmin, dmax, k, Q, mincoverage, mrange, m)
		if(!is.null(result) && (is.null(best) || result$score > best$score))
		{
			best <- result
		}
	}
	seq(best$lmin, best$lmax, by=best$lstep)
}

.wilkinson.nice.scale <- function(min, max, k, Q = c(1,5,2,2.5,3,4,1.5,7,6,8,9), mincoverage = 0.8, mrange=c(), m=k)
{
	Q <- c(10, Q)

	range <- max-min
	intervals <- k-1
	granularity <- 1 - abs(k-m)/m

	delta <- range / intervals
	base <- floor(log10(delta))
	dbase <- 10^base

	best <- NULL
	for(i in 1:length(Q))
	{
		tdelta <- Q[i] * dbase
		tmin <- floor(min/tdelta) * tdelta
		tmax <- tmin + intervals * tdelta

		if(tmin <= min && tmax >= max)
		{
			roundness <- 1 - ((i-1) - ifelse(tmin <= 0 && tmax >= 0, 1, 0)) / length(Q)
			coverage <- (max-min)/(tmax-tmin)
			if(coverage > mincoverage)
			{
				tnice <- granularity + roundness + coverage

				## Wilkinson's implementation contains code to favor certain ranges of labels
				## e.g. those balanced around or anchored at 0, etc.
				## We did not evaluate this type of optimization in the paper, so did not include it.
				## Obviously this optimization component could also be added to our function.
				#if(tmin == -tmax || tmin == 0 || tmax == 1 || tmax == 100)
				#	tnice <- tnice + 1
				#if(tmin == 0 && tmax == 1 || tmin == 0 && tmax == 100)
				#	tnice <- tnice + 1

				if(is.null(best) || tnice > best$score)
				{
					best <- list(lmin=tmin,
							 lmax=tmax,
							 lstep=tdelta,
					      	 score=tnice
						)
				}
			}
		}
	}
	best
}




## The Extended-Wilkinson algorithm described in the paper.

## Our scoring functions, including the approximations for limiting the search
.simplicity <- function(q, Q, j, lmin, lmax, lstep)
{
	eps <- .Machine$double.eps * 100

	n <- length(Q)
	i <- match(q, Q)[1]
	v <- ifelse( (lmin %% lstep < eps || lstep - (lmin %% lstep) < eps) && lmin <= 0 && lmax >=0, 1, 0)

	1 - (i-1)/(n-1) - j + v
}

.simplicity.max <- function(q, Q, j)
{
	n <- length(Q)
	i <- match(q, Q)[1]
	v <- 1

	1 - (i-1)/(n-1) - j + v
}

.coverage <- function(dmin, dmax, lmin, lmax)
{
	range <- dmax-dmin
	1 - 0.5 * ((dmax-lmax)^2+(dmin-lmin)^2) / ((0.1*range)^2)
}

.coverage.max <- function(dmin, dmax, span)
{
	range <- dmax-dmin
	if(span > range)
	{
		half <- (span-range)/2
		1 - 0.5 * (half^2 + half^2) / ((0.1 * range)^2)
	}
	else
	{
		1
	}
}

.density <- function(k, m, dmin, dmax, lmin, lmax)
{
	r <- (k-1) / (lmax-lmin)
	rt <- (m-1) / (max(lmax,dmax)-min(dmin,lmin))
	2 - max( r/rt, rt/r )
}

.density.max <- function(k, m)
{
	if(k >= m)
		2 - (k-1)/(m-1)
	else
		1
}

.legibility <- function(lmin, lmax, lstep)
{
	1			## did all the legibility tests in C#, not in R.
}


#' An Extension of Wilkinson's Algorithm for Position Tick Labels on Axes
#'
#' \code{extended} is an enhanced version of Wilkinson's optimization-based axis labeling approach. It is described in detail in our paper. See the references.
#' 
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @param Q set of nice numbers
#' @param only.loose if true, the extreme labels will be outside the data range
#' @param w weights applied to the four optimization components (simplicity, coverage, density, and legibility)
#' @return vector of axis label locations
#' @references
#' Talbot, J., Lin, S., Hanrahan, P. (2010) An Extension of Wilkinson's Algorithm for Positioning Tick Labels on Axes, InfoVis 2010.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
extended <- function(dmin, dmax, m, Q=c(1,5,2,2.5,4,3), only.loose=FALSE, w=c(0.25,0.2,0.5,0.05))
{
	eps <- .Machine$double.eps * 100
	
	if(dmin > dmax) {
		temp <- dmin
		dmin <- dmax
		dmax <- temp
	}

	if(dmax - dmin < eps) {
		#if the range is near the floating point limit,
		#let seq generate some equally spaced steps.
		return(seq(from=dmin, to=dmax, length.out=m))
	}

	n <- length(Q)

	best <- list()
	best$score <- -2
	
	j <- 1
	while(j < Inf)
	{
		for(q in Q)
		{
			sm <- .simplicity.max(q, Q, j)

			if((w[1]*sm+w[2]+w[3]+w[4]) < best$score)
			{
				j <- Inf
				break
			}
		
			k <- 2
			while(k < Inf)													# loop over tick counts
			{		
				dm <- .density.max(k, m)

				if((w[1]*sm+w[2]+w[3]*dm+w[4]) < best$score)
					break
			
				delta <- (dmax-dmin)/(k+1)/j/q
				z <- ceiling(log(delta, base=10))

				while(z < Inf)
				{			
					step <- j*q*10^z

					cm <- .coverage.max(dmin, dmax, step*(k-1))

					if((w[1]*sm+w[2]*cm+w[3]*dm+w[4]) < best$score)
						break
					
					min_start <- floor(dmax/(step))*j - (k - 1)*j
					max_start <- ceiling(dmin/(step))*j

					if(min_start > max_start)
					{
						z <- z+1
						next
					}

					for(start in min_start:max_start)
					{
						lmin <- start * (step/j)
						lmax <- lmin + step*(k-1)
						lstep <- step

						s <- .simplicity(q, Q, j, lmin, lmax, lstep)
						c <- .coverage(dmin, dmax, lmin, lmax)						
						g <- .density(k, m, dmin, dmax, lmin, lmax)
						l <- .legibility(lmin, lmax, lstep)						

						score <- w[1]*s + w[2]*c + w[3]*g + w[4]*l

						if(score > best$score && (!only.loose || (lmin <= dmin && lmax >= dmax)))
						{
							best <- list(lmin=lmin,
								 lmax=lmax,
								 lstep=lstep,
						         score=score)
						}
					}
					z <- z+1
				}				
				k <- k+1
			}
		}
		j <- j + 1		
	}

	seq(from=best$lmin, to=best$lmax, by=best$lstep)
}



## Quantitative evaluation plots (Figures 2 and 3 in the paper)


#' Generate figures from An Extension of Wilkinson's Algorithm for Position Tick Labels on Axes
#'
#' Generates Figures 2 and 3 from our paper.
#' 
#' @param samples number of samples to use (in the paper we used 10000, but that takes awhile to run).
#' @return produces plots as a side effect
#' @references
#' Talbot, J., Lin, S., Hanrahan, P. (2010) An Extension of Wilkinson's Algorithm for Positioning Tick Labels on Axes, InfoVis 2010.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
extended.figures <- function(samples = 100)
{
	oldpar <- par()
	par(ask=TRUE)
	
	a <- runif(samples, -100, 400)
	b <- runif(samples, -100, 400)
	low <- pmin(a,b)
	high <- pmax(a,b)
	ticks <- runif(samples, 2, 10)

	generate.labelings <- function(labeler, dmin, dmax, ticks, ...)
	{
		mapply(labeler, dmin, dmax, ticks, SIMPLIFY=FALSE, MoreArgs=list(...))
	}
	
	h1 <- generate.labelings(heckbert, low, high, ticks)
	w1 <- generate.labelings(wilkinson, low, high, ticks, mincoverage=0.8)
	f1 <- generate.labelings(extended, low, high, ticks, only.loose=TRUE)
	e1 <- generate.labelings(extended, low, high, ticks)
	
	figure2 <- function(r, names)
	{
		for(i in 1:length(r))
		{
			d <- r[[i]]
			
			#plot coverage
			cover <- sapply(d, function(x) {max(x)-min(x)})/(high-low)
			hist(cover, breaks=seq(from=-0.01,to=1000,by=0.02), xlab="", ylab=names[i], main=ifelse(i==1, "Density", ""), col="darkgray", lab=c(3,3,3), xlim=c(0.5,3.5), ylim=c(0,0.12*samples), axes=FALSE, border=FALSE)
			#hist(cover)
			axis(side=1, at=c(0,1,2,3,4), xlab="hello", line=-0.1, lwd=0.5)
			
			# plot density
			dens <- sapply(d, length) / ticks
			hist(dens, breaks=seq(from=-0.01,to=10,by=0.02), xlab="", ylab=names[i], main=ifelse(i==1, "Density", ""), col="darkgray", lab=c(3,3,3), xlim=c(0.5,3.5), ylim=c(0,0.06*samples), axes=FALSE, border=FALSE)
			axis(side=1, at=c(0,1,2,3,4), xlab="hello", line=-0.1, lwd=0.5)
		}
	}

	par(mfrow=c(4, 2), mar=c(0.5,1.85,1,0), oma=c(1,0,1,0), mgp=c(0,0.5,-0.3), font.main=1, font.lab=1, cex.lab=1, cex.main=1, tcl=-0.2)
	figure2(list(h1,w1, f1, e1), names=c("Heckbert", "Wilkinson", "Extended\n(loose)", "Extended\n(flexible)"))

	figure3 <- function(r, names)
	{
		for(i in 1:length(r))
		{
			d <- r[[i]]
			steps <- sapply(d, function(x) round(median(diff(x)), 2))
			steps <- steps / (10^floor(log10(steps)))
			tab <- table(steps)
			barplot(rev(tab), xlim=c(0,0.4*samples), horiz=TRUE, xlab=ifelse(i==1,"Frequency",""), xaxt='n', yaxt='s', las=1, main=names[i], border=NA, col="gray")
		}
	}
	
	par(mfrow=c(1,4), mar=c(0.5, 0.75, 2, 0.5), oma=c(0,2,1,1), mgp=c(0,0.75,-0.3), cex.lab=1, cex.main=1)
	figure3(list(h1,w1, f1, e1), names=c("Heckbert", "Wilkinson", "Extended\n(loose)", "Extended\n(flexible)"))
	par(oldpar)
}



#' Nelder's labeling algorithm
#'
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @param Q set of nice numbers
#' @return vector of axis label locations
#' @references
#' Nelder, J. A. (1976) AS 96. A Simple Algorithm for Scaling Graphs, Journal of the Royal Statistical Society. Series C., pp. 94-96.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
nelder <- function(dmin, dmax, m, Q = c(1,1.2,1.6,2,2.5,3,4,5,6,8,10))
{
	ntick <- floor(m)
	tol <- 5e-6
	bias <- 1e-4

	intervals <- m-1
	x <- abs(dmax)
	if(x == 0) x <- 1
	if(!((dmax-dmin)/x > tol))
	{
		## special case handling for very small ranges. Not implemented yet.
	}

	step <- (dmax-dmin)/intervals
	s <- step

	while(s <= 1)
		s <- s*10
	while(s > 10)
		s <- s/10

	x <- s-bias
	unit <- 1
	for(i in 1:length(Q))
	{
		if(x < Q[i])
		{
			unit <- i
			break
		}
	}
	step <- step * Q[unit] / s
	range <- step*intervals

	x <- 0.5 * (1+ (dmin+dmax-range) / step)
	j <- floor(x-bias)
	valmin <- step * j

	if(dmin > 0 && range >= dmax)
		valmin <- 0
	valmax <- valmin + range

	if(!(dmax > 0 || range < -dmin))
	{
		valmax <- 0
		valmin <- -range
	}

	seq(from=valmin, to=valmax, by=step)
}


#' R's pretty algorithm implemented in R
#'
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @param n number of axis intervals (specify one of \code{m} or \code{n})
#' @param min.n nonnegative integer giving the \emph{minimal} number of intervals. If \code{min.n == 0}, \code{pretty(.)} may return a single value.
#' @param shrink.sml positive numeric by a which a default scale is shrunk in the case when \code{range(x)} is very small (usually 0).
#' @param high.u.bias non-negative numeric, typically \code{> 1}. The interval unit is determined as \code{\{1,2,5,10\}} times \code{b}, a power of 10. Larger \code{high.u.bias} values favor larger units.
#' @param u5.bias non-negative numeric multiplier favoring factor 5 over 2. Default and 'optimal': \code{u5.bias = .5 + 1.5*high.u.bias}.
#' @return vector of axis label locations
#' @references
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
rpretty <- function(dmin, dmax, m=6, n=floor(m)-1, min.n=n%/%3, shrink.sml = 0.75, high.u.bias=1.5, u5.bias=0.5 + 1.5*high.u.bias)
{
	ndiv <- n
	h <- high.u.bias
	h5 <- u5.bias

	dx <- dmax-dmin
	if(dx==0 && dmax==0)
	{
		cell <- 1
		i_small <- TRUE
		U <- 1
	}
	else
	{
		cell <- max(abs(dmin), abs(dmax))
		U <- 1 + ifelse(h5 >= 1.5*h+0.5, 1/(1+h), 1.5/(1+h5))
		i_small = dx < (cell * U * max(1, ndiv) * 1e-07 * 3)
	}

	if(i_small)
	{
		if(cell > 10)
		{
			cell <- 9+cell/10
		}
    	cell <- cell * shrink.sml
		if(min.n > 1) cell <- cell/min.n
	}	
	else
	{
		cell <- dx
		if(ndiv > 1) cell <- cell/ndiv
	}

	if(cell < 20 * 1e-07)
		cell <- 20 * 1e-07
	
	base <- 10^floor(log10(cell))

	unit <- base

	if((2*base)-cell < h*(cell-unit))
	{
		unit <- 2*base
		if((5*base)-cell < h5*(cell-unit))
		{
			unit <- 5*base
			if((10*base)-cell < h*(cell-unit))
				unit <- 10*base
		}
	}

	# track down lattice labelings...

	## Maybe used to correct for the epsilon here??
	ns <- floor(dmin/unit + 1e-07)
	nu <- ceiling(dmax/unit - 1e-07)

	## Extend the range out beyond the data. Does this ever happen??
	while(ns*unit > dmin+(1e-07*unit)) ns <- ns-1
	while(nu*unit < dmax-(1e-07*unit)) nu <- nu+1


	## If we don't have quite enough labels, extend the range out to make more (these labels are beyond the data :( )
	k <- floor(0.5 + nu-ns)
	if(k < min.n)
	{
		k <- min.n - k
		if(ns >=0)
		{
			nu <- nu + k/2
			ns <- ns - k/2 + k%%2
		}
		else
		{
			ns <- ns - k/2
			nu <- nu + k/2 + k%%2
		}
		ndiv <- min.n
	}
	else
	{
		ndiv <- k
	}

	graphmin <- ns*unit
	graphmax <- nu*unit

	seq(from=graphmin, to=graphmax, by=unit)
}

#' Matplotlib's labeling algorithm
#'
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @return vector of axis label locations
#' @references
#' \url{http://matplotlib.sourceforge.net/}
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
matplotlib <- function(dmin, dmax, m)
{
	steps <- c(1,2,5,10)
	nbins <- m
	trim <- TRUE

	vmin <- dmin
	vmax <- dmax
	params <- .matplotlib.scale.range(vmin, vmax, nbins)
	scale <- params[1]
	offset <- params[2]

	vmin <- vmin-offset
	vmax <- vmax-offset

	rawStep <- (vmax-vmin)/nbins
	scaledRawStep <- rawStep/scale

	bestMax <- vmax
	bestMin <- vmin

	scaledStep <- 1
	chosenFactor <- 1

	for (step in steps)
	{
		if (step >= scaledRawStep)
		{
			scaledStep <- step*scale
			chosenFactor <- step
			bestMin <- scaledStep * floor(vmin/scaledStep)
			bestMax <- bestMin + scaledStep*nbins
			if (bestMax >= vmax)
				break
		}
	}
	if (trim)
	{
		extraBins <- floor((bestMax-vmax)/scaledStep)
		nbins <- nbins-extraBins
	}
	graphMin <- bestMin+offset
	graphMax <- graphMin+nbins*scaledStep

	seq(from=graphMin, to=graphMax, by=scaledStep)
}

.matplotlib.scale.range <- function(min, max, bins)
{
	threshold <- 100
	dv <- abs(max-min)
	maxabsv<-max(abs(min), abs(max))
	if (maxabsv == 0 || dv/maxabsv<10^-12)
		return(c(1, 0))

	meanv <- 0.5*(min+max)

	if ((abs(meanv)/dv) < threshold)
		offset<- 0
	else if (meanv>0)
	{
		exp<-floor(log10(meanv))
		offset = 10.0^exp
	} else
	{
		exp <- floor(log10(-1*meanv))
		offset <- -10.0^exp
	}
	exp <- floor(log10(dv/bins))
	scale = 10.0^exp
	c(scale, offset)
}



#' gnuplot's labeling algorithm
#'
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @return vector of axis label locations
#' @references
#' \url{http://www.gnuplot.info/}
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
gnuplot <- function(dmin, dmax, m)
{
	ntick <- floor(m)
	power <- 10^floor(log10(dmax-dmin))
	norm_range <- (dmax-dmin)/power
	p <- (ntick-1) / norm_range

	if(p > 40)
		t <- 0.05
	else if(p > 20)
		t <- 0.1
	else if(p > 10)
		t <- 0.2
	else if(p > 4)
		t <- 0.5
	else if(p > 2)
		t <- 1
	else if(p > 0.5)
		t <- 2
	else
		t <- ceiling(norm_range)

	d <- t*power
	graphmin <- floor(dmin/d) * d
	graphmax <- ceiling(dmax/d) * d

	seq(from=graphmin, to=graphmax, by=d)
}



#' Sparks' labeling algorithm
#'
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @return vector of axis label locations
#' @references
#' Sparks, D. N. (1971) AS 44. Scatter Diagram Plotting, Journal of the Royal Statistical Society. Series C., pp. 327-331.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
sparks <- function(dmin, dmax, m)
{
	fm <- m-1
	ratio <- 0
	key <- 1
	kount <- 0
	r <- dmax-dmin
	b <- dmin
	
	while(ratio <= 0.8)
	{
		while(key <= 2)
		{
			while(r <= 1)
			{
				kount <- kount + 1
				r <- r*10
			}
			while(r > 10)
			{
				kount <- kount - 1
				r <- r/10
			}

			b <- b*(10^kount)
			if( b < 0 && b != trunc(b)) b <- b-1
			b <- trunc(b)/(10^kount)
			r <- (dmax-b)/fm
			kount <- 0
			key <- key+2
		}
	
		fstep <- trunc(r)
		if(fstep != r) fstep <- fstep+1
		if(r < 1.5) fstep <- fstep-0.5
		fstep <- fstep/(10^kount)
		ratio <- (dmax - dmin)*(fm*fstep)
		kount <- 1
		key <- 2
	}
	fmin <- b
	c <- fstep*trunc(b/fstep)
	if(c < 0 && c != b) c <- c-fstep
	if((c+fm*fstep) > dmax) fmin <- c
	
	seq(from=fmin, to=fstep*(m-1), by=fstep)
}


#' Thayer and Storer's labeling algorithm
#'
#' @param dmin minimum of the data range
#' @param dmax maximum of the data range
#' @param m number of axis labels
#' @return vector of axis label locations
#' @references
#' Thayer, R. P. and Storer, R. F. (1969) AS 21. Scale Selection for Computer Plots, Journal of the Royal Statistical Society. Series C., pp. 206-208.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
#' @export
thayer <- function(dmin, dmax, m)
{
	r <- dmax-dmin
	b <- dmin
	kount <- 0
	kod <- 0

	while(kod < 2)
	{
		while(r <= 1)
		{
			kount <- kount+1
			r <- r*10
		}
		while(r > 10)
		{
			kount <- kount-1
			r <- r/10
		}
		b <- b*(10^kount)
		if(b < 0)
			b <- b-1
		ib <- trunc(b)
		b <- ib
		b <- b/(10^kount)
		r <- dmax-b
		a <- r/(m-1)
		kount <- 0
		while(a <= 1)
		{
			kount <- kount+1
			a <- a*10
		}
		while(a > 10)
		{
			kount <- kount-1
			a <- a/10
		}
		ia <- trunc(a)
		if(ia == 6) ia <- 7
		if(ia == 8) ia <- 9

		aa <- 0
		if(a < 1.5) aa <- -0.5
		a <- aa + 1 + ia
		a <- a/(10^kount)
				
		test <- (m-1) * a
		test1 <- (dmax-dmin)/test
		if(test1 > 0.8)
			kod <- 2

		if(kod < 2)
		{
			kount <- 1
			r <- dmax-dmin
			b <- dmin
			kod <- kod + 1
		}
	}

	iab <- trunc(b/a)
	if(iab < 0) iab <- iab-1
	c <- a * iab
	d <- c + (m-1)*a
	if(d >= dmax)
		b <- c

	valmin <- b
	valmax <- b + a*(m-1)	

	seq(from=valmin, to=valmax, by=a)
}

