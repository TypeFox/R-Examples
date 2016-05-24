"hftrialdatagen" <-
function (nreps=4, nps=128,plot.it=FALSE,uvp=0.8) 
{
#
# Number of replicates and `probe sets' (although this is one
# channel cDNA, not dual) is set by arguments
#
#
# Simulate nps probe sets with nreps replications each. Basically, this is
# simulating nps \mu values (where each \mu is Gamma distributed with shape
# and scale parameter specified).
#
m <- genesimulator(nreps=nreps, nps=nps,shape=0.8, scale=2000000)
#
# Setup DHHR parameters
#
alpha <- 24800
seta <- 0.227
seps <- 4800
Seta <- sqrt(exp(seta^2)*( exp(seta^2)-1))
cc <- seps^2 / Seta^2
#
# Use the Durbin and Rocke mixed log-normal/Gaussian distribution with the
# parameters specified in their J. Comp Biol. 2001 and Bioinformatics 2002
# papers to simulate y values using the \mu values already obtained.
#
y <- simdurbin2(mu=m[,1],alpha=alpha,seta=seta, seps=seps)
#
# Attach these simulated y values to the matrix of \mu, replicate index,
# probe index.
#
m <- cbind(m,y)
#
# Now for each probe set work out the mean of the y for each replicate. I.e.
# fix a probe set and collect together all replicates for that probe set and
# compute their mean (we assume that all experiments underlying a replicate
# set originate from the same distribution).
#
# The matrix mm is computed first. Each row indexs a probe set and the columns
# for a row contain the replicates.
#
mm <- matrix(0, nrow=nps, ncol=nreps)
for(i in 1:nps)
	mm[i,] <- m[ m[,3]==i, 4]
#
# Work out the mean and standard deviation for the replicates for each
# probe set.
#
mps <- apply(mm,1,mean)
msd <- sqrt(apply(mm,1,var))
#
# If plot.it==TRUE then this plot reproduces Figure 1 in Rocke and Durbin
# in J. Comp. Biol. 2001 paper (with simulated data, not real data). 
#
if (plot.it==TRUE)	{
	plot(mps, msd)
	scan()
	}
#
# Now reorder the rows of the mm matrix according to increasing mean value
# stored in a sorted mps 
#
neworder <- sort.list(mps)
mmrowsort <- mm[neworder,]
#
# Now sort the elements within each row of mmrowsort
# I would have liked to use apply(mmrowsort, 1, sort) but this doesn't work
#
mmallsort <- t(apply(t(mmrowsort), 2, sort))
#
# Now create a matrix with equivalent format to mmallsort and mmrowsort but
# containing the underlying true \mu values for comparison
#
muasmm <- matrix(mm[,1], nrow=nps, ncol=nreps)[neworder,]
#
# Now plot the y data that will be smoothed
#
n <- nreps*nps	# Total number of observations
#
# Create the data to smooth with the data driven Haar Fisz transform DDHFT
# (basically create vectors out of the matrices)
#
yhf <- as.vector(t(mmrowsort))	# Noisy data
muhf <- as.vector(t(muasmm))	# The truth
ymasmm <- rep(sort(mps), rep(nreps,nps))	# Y mean over replicate ints
if (plot.it==TRUE)	{
	plot(1:n, yhf, xlab="Mu", ylab="Y", sub="Red=True Mu, Green=Ymean")
	lines(1:n, muhf, col=2)
	lines(1:n, ymasmm, col=3)
	scan()
	} 
#
# Do the DDHFT
#
yddhft <- ddhft.np.2(yhf)
#
# Plot the DDHFT estimated mean/variance function
#
if (plot.it==TRUE)	{
	px <- yddhft$mu
	ply <- yddhft$sigma^2
	plot(px, yddhft$sigma2, log="xy", xlab="Log(mu)", ylab="Log(sigma)", sub="Line is best increasing sigma" )
	lines(px, ply, col=2)
	lines(px[c(1, length(px))], ply[c(1,length(ply))], col=3)
	scan()
	}
#
# Plot the DDHFT transformed `signal plus noise'
#
yhftm <- matrix(yddhft$hft, nrow=nps, ncol=nreps, byrow=TRUE)
yhftmm <- apply(yhftm, 1, mean)
yhftms <- sqrt(apply(yhftm, 1, var))
if (plot.it==TRUE)	{
	plot(yhftmm, yhftms, xlab="(Transformed) mean of reps",
		ylab="(Transformed) scale of reps",
		sub="cf Figure 4, Durbin et al., Bioinformatics, (2002)")
	scan()
	}
if (plot.it==TRUE)	{
	plot(1:n, yddhft$hft, ylab="DDHF transformed signal", sub="Red=Haar Wavelet, Green=Smooth Spline")
	}

#
# Now denoise the signal
#
# With non-decimated Haar wavelets with universal thresholding by 0.8
#
# Reflection is to make series periodic for denoising
#

perhft <- c(yddhft$hft, rev(yddhft$hft))

yhftwst <- wst(perhft, filter.number=1, family="DaubExPhase")
yhftwstT <- uvp*threshold(yhftwst, by.level=TRUE, policy="universal", return=TRUE)
yhftwstT <- threshold(yhftwst, by.level=TRUE, policy="manual", value=yhftwstT)

#return(yhftwstT)
yhftwr <- AvBasis(yhftwstT)

yhftwr <- yhftwr[1:n]

#
# Alternative smoothing using smooth.spline
#
yss <- smooth.spline(x=1:n, y=yddhft$hft)
	
if (plot.it==TRUE)	{
	lines(1:n, yhftwr, col=2)
	lines(1:n, yss$y, col=3)
	scan()
	}
#
# Do inverse Haar-Fisz transform
#
yddhft$hft <- yhftwr
yihf <- ddhft.np.inv(yddhft)

if (plot.it==TRUE)	{
	plot(1:n, yihf, xlab="Mu", ylab="Y", sub="Red=True Mu, Green=Ymean")
	lines(1:n, muhf, col=2)
	lines(1:n, ymasmm, col=3)
	} 
#
# Convert answer back into array like mm
#
yihfm <- matrix(yihf, nrow=nps, ncol=nreps, byrow=TRUE)
yihfmm <- apply(yihfm, 1, mean)
#
#cat("Length yihfmm" , length(yihfmm), "\n")
#cat("Length muasmm[,1]" , length(muasmm[,1]), "\n")



ansm <- cbind(yihfmm, muasmm[,1])
dimnames(ansm) <- list(NULL, c("HF", "TrueMu"))
hftssq <- (yihfmm-muasmm[,1])^2
l <- list(ansm=ansm, hftssq=hftssq, yhf=yhf)
l
}

