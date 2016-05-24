hwtos2 <-
function(x, alpha=0.05, filter.number=1, family="DaubExPhase", lowlev=3,
    WTscale=NULL, maxSD=NULL, verbose=FALSE, silent=FALSE,
    UseCForVarip2=TRUE, OPLENGTH=100000){

if (any(is.na(x)))
	stop("Some values are NA. This is not permitted")

lx <- length(x)

if (!IsPowerOfTwo(lx))
	stop("Data length is not power of two. Please consider using the hwtos function instead")

xS <- ewspecHaarNonPer(x, filter.number=filter.number, family=family, WPsmooth=FALSE)

xWP <- xS$WavPer
xS <- xS$S

J <- xS$nlevels

if (is.null(WTscale))
	WTscale <- ceiling(J/2)-1 

xSD <- apply(matrix(xS$D, nrow=J, byrow=TRUE),1, mean)

xSD[xSD < 0] <- 0

#print("TEMPORARY MEASURE: USING PREDEFINED xSD")
#xSD <- 2^(-(1:length(xSD)))
#print(xSD)

if (verbose==TRUE)	{
	cat("xSD: ")
	print(xSD)
	}

if (is.null(maxSD))
	maxSD <- length(xSD)

AllTS <- AllPVal <- vector("list", J)

P <- PsiJ(-maxSD, filter.number = filter.number, family = family, OPLENGTH=OPLENGTH)
Pmat <- PsiJmat(-maxSD, filter.number = filter.number, family = family, OPLENGTH=OPLENGTH)


for(i in (J-1):lowlev)	{

	if (silent==FALSE)
		cat(" ", i, " ")

	if (verbose==TRUE)
		cat("Looking at periodogram level: ", i, "\n")

	#
	# Get the raw periodogram coefficients at a level
	#
	the.wp <- accessD(xWP, lev=i)
	#
	# Do the Haar WT  (note tmp.wp.pval is a placeholder)
	#
	the.wp.pval <- the.wp.hwt <- wd(the.wp, filter.number=1, family="DaubExPhase")

	the.wp.pval$D <- rep(1, length(the.wp.pval$D))

	for(k in (J-WTscale-1):0)	{

		if (verbose==TRUE)
			cat("Looking at HWT level: ", k, "\n")

		v <- accessD(the.wp.hwt, level=k)

		#cat("varip2: i=", J-k, " ll=", J-i, "\n")
		if (UseCForVarip2==TRUE)
			the.varip <- sqrt(Cvarip2(i=J-k, p=5, ll=J-i, S=xSD[1:maxSD], Pmat=Pmat, PsiJL=sapply(P, "length")))
		else
			the.varip <- sqrt(varip2(i=J-k, p=5, ll=J-i, S=xSD[1:maxSD], P=P)$ans)
		the.varip.new <- sqrt(littlevar(xWP, ll=i))

		if (verbose==TRUE)
			cat("k: ", k, " New Var Est: ", the.varip.new, " (Old is: ", the.varip, ")\n")

		if (length(v) > 4)
			data.sd <- sd(v)
		else
			data.sd <- the.varip

		#cat("SD: Th1:Th2:Dat is:", the.varip, the.varip.new, data.sd, "\n")

		the.varip <- max(the.varip.new, the.varip, data.sd)

		if (the.varip==0)
			stop("Variance is zero. Is your time series a constant function? This function only works on stochastic series")

		StudT <- v/the.varip

		if (any(abs(StudT) > 1000))	{
			cat("Pgram level i: ", i,", HWT lev k: ", k, " max: ", max(abs(StudT)), "\n")
			print(v)
			cat("the.varip: ", the.varip, "\n")
			cat("--------------\n")
			}


		the.wp.hwt <- putD(the.wp.hwt, level=k, v=StudT)
		the.wp.pval <- putD(the.wp.pval, level=k, v=2*pnorm(-abs(StudT)))
		}
	AllTS[[i+1]] <- the.wp.hwt
	AllPVal[[i+1]] <- the.wp.pval
	}

if (silent==FALSE)
	cat("\n")

#
# Now collect all relevant pvals together
#

allTS <- allpvals <- NULL
for(i in (J-1):lowlev)	{
	for(j in (J-WTscale-1):0)	{
		allTS <- c(allTS, accessD(AllTS[[i+1]], lev=j))
		}
	}

allpvals <- 2*pnorm(-abs(allTS))

m <- length(allpvals)

if (verbose==TRUE)
	cat("There are ", m, " p values \n")

opvals <- sort.list(allpvals)
spvals <- allpvals[opvals]
sTS <- allTS[opvals]
#
# Find largest k such that spvals[k] <= k*alpha/m
#

kstore <- 0
k <- m
while(kstore == 0 && k > 0)	{
	if (spvals[k] <= k*alpha/m)
		kstore <- k
	else
		k <- k-1
	}
nreject <- kstore
rejpval <- spvals[kstore]

l <- list(nreject=nreject, rejpval=rejpval, spvals=spvals, sTS=sTS,
	AllTS=AllTS, AllPVal=AllPVal, alpha=alpha, x=x, xSD=xSD)

class(l) <- "tos"
l

}
