### $Id: baseline.peakDetection.R 170 2011-01-03 20:38:25Z bhm $
baseline.peakDetection <- function(spectra, left, right, lwin, rwin, snminimum, mono=0, multiplier=5, left.right, lwin.rwin){
# peakDetection(spectra, L, R, LB, RB, sn, mono, multiplier)
#
# Given a spectrum, we want to remove the baseline and identify the
# location of peaks. The main return from this function is a list,
# PEAKS, of the peak locations in ticks. The second return value is
# the modified spectrum with the baseline removed.
#
# The second argument, snminimum, is the minimum acceptable value for
# the signal to noise ratio. The next pair of values are the left and
# right window sizes for peak resolution; see .mixedup for details. The
# next pair of values are the left and right window sizes for peak
# removal in preparation for baseline correction; see .localMinPet for
# details. The next value specifies whether the baseline should be made
# monotonically decreasing. The final value is a multiplier specifying
# the number of peak widths used to calculate the noise; this defaults to 5.
#
# Described as algorithm SPDBC in Clinical Chemistry 2003; 49:1615-1623.
#
# Original version by Kevin R. Coombes
# $Revision: 2.2 $
# Last modified by $Author: krc $ on $Date: 2003/03/31 19:01:58 $.

np <- dim(spectra)
Midspec <- corrected <- matrix(0,np[1],np[2])
Sn <- Peaks <- Y1 <- Y2 <- Y3 <- list()
if(!missing(left.right)){
	left <- left.right
	right <- left.right
}
if(!missing(lwin.rwin)){
	lwin <- lwin.rwin
	rwin <- lwin.rwin
}

for(a in 1:np[1]){
	spectrum <- spectra[a,]
	# step 1: Use our basic peak finder on the raw spectrum.
	mixr <- .mixedup(spectrum, left, right)
	y <- mixr$topper; yleft <- mixr$lotter; yright <- mixr$rotter

	# step 2: Remove the peaks and fill in with straight lines
	# across the base.
	slip <- .peakRemoval(spectrum, y, yleft, yright, mono)

	# step 3: Compute a local minimum in a window of width 256. (The
	# argument to the function is 128, which extends on each side from
	# the center.)
	mini1 <- .localMinPet(slip, lwin, rwin, mono)

	# step 4: Remove the local minimum to estimate the base line.
	midspec <- spectrum - mini1

	# step 5: Repeat steps one through four. The first pass typically
	# leaves a slight slant to the baseline in the exponentially decaying
	# region, and a second pass cleans up most of this.
	mixr <- .mixedup(midspec, left, right)
	y2 <- mixr$topper; yleft2 <- mixr$lotter; yright2 <- mixr$rotter
	slip2 <- .peakRemoval(midspec, y2, yleft2, yright2, mono)
	mini <- .localMinPet(slip2, lwin, rwin, mono)
	newspec <- midspec - mini
	corrected[a,] <- newspec

	if(!missing(snminimum)){
		# step 6: Make a final peak finding pass on the baseline-corrected
		# spectrum.
		mixr <- .mixedup(newspec, left, right)
		y3 <- mixr$topper; yleft3 <- mixr$lotter; yright3 <- mixr$rotter

		# step 7: In general, the peak finder is too inclusive, since it does
		# not take noise into consideration. As described in that file, the
		# general process found about 450 peaks in a spectrum of 18000 points.
		# Visual inspection revealed that many of them were still in the noise
		# region. So, we compute a signal to noise ratio using the MAD in a
		# local window, and only keep peaks above a threshold. Using S/N > 3,
		# we get down to about 200 peaks, most of which seem believeable.
		width <- multiplier*pmax(yright3-y3, y3-yleft3)
		noise <- 0*y3
		middle <- 0*y3
		for(i in 1:length(y3)){
		    lend <- max(1, (y3[i]-width[i]))
		    rend <- min(length(newspec), (y3[i]+width[i]))
		    slice <- c(newspec[lend:yleft3[i]], newspec[yright3[i]:rend])
		    middle[i] <- median(slice)
		    noise[i] <- median(abs(slice-middle[i]))
		    noise[i] <- median(abs(diff(newspec[lend:rend])))
		}
		signal <- newspec[y3]-middle

		sn <- pmax(0, signal/noise)

		peaks <- y3[sn > snminimum]

		Peaks[[a]] <- peaks
		Sn[[a]] <- sn
		Y1[[a]] <- y
		Y2[[a]] <- y2
		Y3[[a]] <- y3
		Midspec[a,] <- midspec
	}
}
# return
if(missing(snminimum)){
	list(baseline=spectra-corrected, corrected=corrected)
} else {
	list(baseline=spectra-corrected, corrected=corrected, peaks=Peaks, sn=Sn, y3=Y3, midspec=Midspec, y=Y1, y2=Y2)
}
}


#########################################################################################################
.peakRemoval <- function(w, y, yleft, yright, mono){
# .peakRemoval(spectrum, y, left, right, mono)
#
# Inputs are a spectrum, and a list of peaks (y) with their
# left (yleft) and right (yright) boundaries. There is another
# optional argument (mono). if(mono is greater than zero, then
# we believe the ultimate baseline should be monotonically
# decreasing. The default value of mono = 0, which does not
# enfore monotonicity. The output is a new spectrum of the
# same length with the peaks removed and replaced by straight
# lines across the base.
#
# Part of the SPF, SPDBC package described in Clinical Chemistry 2003; 49:1615-1623.
#
# Original version by Kevin Coombes
# $Revision: 2.2 $
# Last modified by $Author: krc $ on $Date: 2003/03/31 19:01:58 $.

if(missing(mono))
   mono <- 0

# The problem is what happens if(you have several real peaks
# close together. if(one peak doesn't come all the way down to
# baseline, you get big jumps in the "peak-removed" spectrum.
# We can start chopping those out, but have to adjust the
# left and right boundaries.

magic <- 3	# magic number of multiples of the median slope
slopes <- (w[yright]-w[yleft])/(yright-yleft)
m <- magic*median(abs(slopes))
counters <- 1:length(slopes)
if(mono > 0){
   startup <- min(counters[abs(slopes)<m]) # initial slide
} else {
   temp <- 1:length(w)
   screwy <- max(temp[w==max(w)])
   startup <- min(counters[y > screwy])	# end of saturation
}
m <- magic*median(abs(slopes[counters>startup]))	# ignore beginning
selector <- !(counters>startup & abs(slopes) > m)

# only apply this to the peaks with slopes at most 3 times the median
# however, we can have large slopes at the beginning of the spectrum
# in the matrix noise region.

# now we move the left and right boundaries accordingly
z  <- y
zr <- yright
zl <- yleft
for(i in 2:(length(y)-1)){
   if(selector[i] == 0){ # forget this peak
      if(w[yright[i-1]] > w[yright[i]])
         zr[i-1] <- zr[i]
      if(w[yleft[i+1]] > w[yleft[i]])
         zl[i+1] <- zl[i]
   }
}

# use the selector
z  <- z[selector]
zl <- zl[selector]
zr <- zr[selector]

# finally, remove the peaks
newspec <- w
for(i in 1:length(z)){
	yr <- zr[i]
	yl <- zl[i]
	slope <- (w[yr] - w[yl])/(yr-yl)
	newspec[yl:yr] <- slope*((yl:yr) - yl) + w[yl]
}

# still not done... have to make sure we don't get negative values.
# return
pmin(newspec, w)
}


#########################################################################################################
.mixedup <- function(w, left, right){
# .mixedup(spectrum, smallWindow, largeWindow)
#
# Input consists of a spectrum and two window sizes. The smallWindow is
# used to collapse potential peaks together at the left end of the
# spectrum; ther default size is 5 ticks. The largeWindow is used to
# collapse peaks together at the right end of the spectrum; the default
# size is 30. The window sizes are smoothly interpolated along a
# quadratic curve, which makes them a fixed percentage of the mass on
# that scale.
#
# The return values are:
#	PEAKS: a list of peak locations, in ticks
#	LEFTMIN:	a list of the closest minimum to the left of each peak
#	RIGHTMIN: a list of the closest minimum to the right of each peak
#
# Described as algorithm SPF in Clinical Chemistry 2003; 49:1615-1623.
#
# Original version by Kevin R. Coombes

# set up default values for optional arguments

n <- length(w)-1

if (missing(right)) right <- 30
if (missing(left)) left <- 5

# perform a coherency check in case someone gave us nonsensical values
if(left <= 0)
	left <- 5
if(right < left){
	temp <- right
	right <- left
	left <- temp
}


## Step 1: Find all local maxima and associated local minima.
x <- diff(w)		# first difference
xl <- x[1:(n-1)]	# left side
xr <- x[2:n]		# right side
neartop <- (xl > 0 & xr <= 0) | (xl == 0 & xr < 0)
nearbot <- (xl < 0 & xr >= 0) | (xl == 0 & xr > 0)
xx <- 1:length(x)
nb <- c(1,(xx[nearbot]+1))	# tentative minimum
nt <- xx[neartop]+1			# tentative maximum

## Step 2: Coalesce rising plateaus to find the right-most one.
curt <- 1
curb <- 1
tops <- 0
bots <- 1
numtop <- 1
numbot <- 1
while(curt < length(nt) && curb < length(nb)){
	while(nb[curb]<nt[curt] && curb < length(nb)){
		curb <- curb+1
	}
	bots[numbot] <- nb[curb]
	numbot <- numbot+1
	while(nt[curt] < nb[curb] && curt < length(nt)){
	   curt <- curt + 1
	}
	tops[numtop] <- nt[curt];
    numtop <- numtop + 1
}

# At this point, "tops" and "bots" contain local maxima and minima
# including those that go up for one tick and then immediately go
# back down. The arrays should be the same size, with each minimum
# occuring strictly to the left of each maximum. In a test example,
# there were about 3500 of these local maxima in a spectrum that
# contained 18000 points.

## Step 3: Cull some obviously inadequate candidates. Anything where the
# distance from minimum to maximum is less than the median single step
# size gets removed. In our example, this step got us down to about
# 1500 potential peaks.
fuzz <- median(abs(x))		# typical first difference; a noise estimate
flux <- w[tops] - w[bots]	# size of the step up from bottom to top
tops <- tops[flux > 2*fuzz]

## Step 4: Further culling. We combine maxima that are obviously too
# close together to be resolved. The resolution is determined by
# the arguments specifying window sizes. In our test example, this
# step got us down from about 1500 peaks to about 500.

dt <- diff(tops)	# distance in ticks between successive tops
a <- (right-left)/length(x)^2	# define "close" using y = a*x^2 + b
resolve <- left + (a*tops)*tops	# using parameters left and right.

temtops <- tops
topper <- 1	# list of resolvable maxima.
i <- 1
while(i < length(temtops)-1){
   if(dt[i] < resolve[i]){
      if(w[temtops[i]] < w[temtops[i+1]])
         temtops[i] <- temtops[i+1]
      else
         temtops[i+1] <- temtops[i]
      dt[i+1] <- temtops[i+2]-temtops[i]
	  }
   else
      topper[length(topper)+1] <- temtops[i]
   i <- i + 1
}

# After combining maxima, we then figure out where the closest minima
# are to both left and right. These might not agree for adjacent peaks
# because of the presence of plateaus.
lotter <- rep(1,length(topper)-1)	# list of resolvable minima to left of peak.
rotter <- rep(1,length(topper)-1)	# similar list to right of peak.
for(i in 2:length(topper)){
   temp <- topper[i-1]:topper[i]
   lotter[i] <- max(temp[w[temp]==min(w[temp])])
   rotter[i] <- min(temp[w[temp]==min(w[temp])])
}
rotter <- c(rotter[2:length(rotter)], length(x))

# Step 5: The primary goal of this step is to strip out candidate
# peaks that lie on the slope up to a true peak. We insist that every
# peak include a minimum number of "up ticks" and "down ticks" from
# the maximum to the neighboring minima.

uptick <- topper-lotter
upper <- median(uptick)-median(abs(uptick-median(uptick)))
dntick <- rotter-topper
downer <- median(dntick)-median(abs(dntick-median(dntick)))
magic <- min(upper, downer) + 1
keeper <- !(uptick < magic | dntick < magic)

topper <- topper[keeper]
# rotter <- rotter[keeper]
lotter <- lotter[keeper]

## Step 6: Refine the peak location. We may have culled the wrong
# peak when two peaks were close together. Thus, we again locate the
# local maximum between two successive minima.
bottle <- c(lotter, length(x));
for(i in 1:length(topper)){
   temp <- bottle[i]:bottle[i+1];
   topper[i] <- max(temp[w[temp]==max(w[temp])])
}

# Of course, if the peaks move, we may need to adjust the bases as well.
lotter <- rep(1,length(topper)-1)	# list of resolvable minima to left of peak.
rotter <- rep(1,length(topper)-1)	# similar list to right of peak.
for(i in 2:length(topper)){
   temp <- topper[i-1]:topper[i]
   lotter[i] <- max(temp[w[temp]==min(w[temp])])
   rotter[i] <- min(temp[w[temp]==min(w[temp])])
}
maxright <- min(topper[length(topper)]+100, length(w))
rotter <- c(rotter[2:length(rotter)], maxright)
# return
list(topper=topper, lotter=lotter, rotter=rotter)
}


#########################################################################################################
.localMinPet <- function(spectrum, window, rightWindow, mono){
# .localMinPet(spectrum, leftWindow, rightWindow, mono)
#
# Arguments are a spectrum, optional window sizes that scale
# quadratically from the left window to the right window, and
# an optional flag to tell us if the baseline should be
# monotonically decreasing. Default window sizes both equal 100.
# Default value of mono = 0, meaning
#
# Part of the SPF, SPDBC package described in Clinical Chemistry 2003; 49:1615-1623.
#
# Original version by Kevin Coombes

if(missing(window))
   window <- 100
if(missing(rightWindow))
   rightWindow <- window
if(missing(mono))
   mono <- 0

len <- length(spectrum)

tops <- 1:len
a <- (rightWindow-window)/len^2		# define "close" using y = a*x^2 + b
resolve <- floor(window + (a*tops)*tops)	# using window parameters.

# compute the local minima in the windows. Edge effects to avoid running
# off the ends of the array.
mini <- numeric(len)
lend <- pmax(1,(tops-resolve))
rend <- pmin(len, (tops+resolve))

for(i in 1:len){
   s <- spectrum[lend[i]:rend[i]]
   mini[i] <- min(s)
}

# monotonize
if(mono > 0){
   mono <- numeric(len)
   mono[1] <- mini[1]
   for(i in 2:len){
      mono[i] <- mono[i-1]
      if(mini[i] < mono[i])
         mono[i] <- mini[i]
   }
   mini <- mono
}

# return
mini
}
