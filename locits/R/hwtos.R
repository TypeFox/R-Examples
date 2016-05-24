hwtos <-
function(x, alpha=0.05, lowlev=1, WTscale=NULL, maxSD=NULL, verbose=FALSE,
	silent=FALSE, UseCForVarip2 = TRUE, OPLENGTH=1e+05, mc.method=p.adjust.methods){

lx <- length(x)

if (lx < 20) 
        stop("Time series has to be length 20 or greater. Power might not be good for short series")

if (!is.na(IsPowerOfTwo(lx)))
	if (silent==FALSE)	{
		warning("The length of your data is a power of two. Please use hwtos2")
		}


#
# First take non-decimated Haar transform for arbitary n of x 
#
x.Hndwt <- hwt(x, type="station")
#
#
# Note: coarsest scale is 1. 
#
#
# Now square all of the values of the wavelet coefficients, to obtain
# "any n" raw wavelet periodogram
#
# And also compute time averaged periodogram (assumes stationarity)
#
J <- nlevelsWT(x.Hndwt)
#
# Storage for T statistics and p-values for each level of the WP
#
AllLV <- AllRecord <- AllTS <- AllPVal <- vector("list", J)

if (lowlev < 0)
	stop("lowlev has to be >= 0")
else if (lowlev+1 > J)
	stop(paste("lowlev has to be <= ", J))

if (is.null(WTscale))
	WTscale <- ceiling(J/2) - 1

x.WP <- vector("list", J)
xWP.stat <- rep(0, J)

for(i in 1:J)	{
	vv <- x.Hndwt$d[[i]]^2
	vv <- vv[!is.na(vv)]
	x.WP[[i]] <- vv
	xWP.stat[i] <- mean(x.WP[[i]], na.rm=TRUE)
	}

#
# Now compute time-averaged ewspec
#
A.HaarJ <- ipndacw(-J, filter.number=1, family="DaubExPhase")
xS.stat <- solve(A.HaarJ) %*% rev(xWP.stat)
xS.stat[xS.stat < 0] <- 0
if (verbose==TRUE)	{
	cat("xS.stat: ")
	print(xS.stat)
	}
if (is.null(maxSD))
	maxSD <- length(xS.stat)

#cat("xS.stat: ")
#print(xS.stat)
#print(xWP.stat)

#
# Generate autocorrelation wavelets
# 
P <- PsiJ(-maxSD, filter.number=1, family="DaubExPhase", OPLENGTH=OPLENGTH)
Pmat <- PsiJmat(-maxSD, filter.number=1, family="DaubExPhase", OPLENGTH=OPLENGTH)
#
# Now examine each scale of the wavelet periodogram, x.WP
#
for(i in J:(lowlev+1))	{

	if (silent==FALSE)
		cat(" ", i, " ")
	if (verbose==TRUE)
		cat("Looking at periodogram level: ", i, "\n")
	the.wp <- x.WP[[i]]

	#cat("Length of WP at scale ", i, ", is: ", length(the.wp), "\n")

	#cat("Sum of NA: ", sum(is.na(the.wp)), "\n")

	# Takes any-n Haar WT of wavelet periodogram (& make storage for pvals)
	the.wp.lv <- the.wp.record <- the.wp.pval <- the.wp.hwt <- hwt(the.wp,
		reindex=TRUE)

	

	if (verbose==TRUE)
		cat("Number of levels in the.wp.hwt: ", nlevelsWT(the.wp.hwt), "\n")


	for(k in (J-WTscale-2):0)	{
		if (verbose==TRUE)
			cat("Looking at HWT level: ", k, "\n")

		v <- the.wp.hwt$d[[k+1]]
		lv <- length(v)

		#cat("P: ", i, ", H: ", k, " Length: ", lv, "\n")
		#cat("V values: ")
		#print(v)

		if (lv>0)	{

		  if (verbose==TRUE)
			cat("Length of coef vector: ", lv, "\n")

		  if (UseCForVarip2 == TRUE)	
			the.varip <- sqrt(Cvarip2(i=J-k,
				p=5, ll=J-(i-1), S=xS.stat[1:maxSD],
				Pmat=Pmat, PsiJL = sapply(P, "length")))
		  else
			the.varip <- sqrt(varip2(i=J-k, p=5, ll=J-(i-1),
				S=xS.stat[1:maxSD], P=P)$ans)

		  the.varip.new <- sqrt(2*mean(the.wp^2, na.rm=TRUE))

		  if (length(v[!is.na(v)]) > 4)
			data.sd <- sd(v, na.rm=TRUE)
		  else
			data.sd <- the.varip

		  if (verbose==TRUE)
			cat("k: ", k, " New Sd Est: ", the.varip.new,
				" (Old is: ", the.varip, ", data.sd: ",
				data.sd, ")\n")

		  #the.varip <- max(the.varip.new, the.varip, data.sd)
		  the.varip <- max(the.varip, data.sd)

		

		  #cat("The varip is: ", the.varip, "\n")

		if (the.varip==0)
                        stop("Variance is zero. Is your time series a constant function? This function only works on stochastic series")


		  StudT <- v/the.varip

		  #cat("StudT values: ")
		  #  print(StudT)
		  #  cat("-----\n")

		  if (any(abs(StudT) > 1000)) {
			cat("Pgram level i: ", i, ", HWT lev k: ", k, 
			" max: ", max(abs(StudT)), "\n")
			print(v)
			cat("the.varip: ", the.varip, "\n")
			cat("--------------\n")
			}

		  the.wp.hwt$d[[k+1]] <- StudT

		  if (verbose==TRUE)	{
			cat("Length StudT: ", length(StudT), "And values:\n")
			if (length(StudT) < 4)
				print(StudT)
			}

		  the.wp.pval$d[[k+1]] <- 2*pnorm(-abs(StudT))
		  the.wp.lv$c[[k+1]] <- rep(lv, length(StudT))

		  the.wp.record$c[[k+1]] <- rep(k, length(StudT))
		  the.wp.record$d[[k+1]] <- 1:length(StudT)

		  if (verbose==TRUE)	{
			cat("lv: ", lv, "\n")
			cat("Length StudT: ", length(StudT), "And values:\n")
			if (length(StudT) < 4)	{
				print(StudT)
				print(the.wp.record$d[[k+1]])
				}
			}
		}
		else	{	# lv <= 0
			the.wp.hwt$d[k+1] <- list(NULL)
			the.wp.pval$d[k+1] <- list(NULL)
			the.wp.lv$c[[k+1]] <- 0
			the.wp.record$c[k+1] <- list(NULL)
			the.wp.record$d[k+1] <- list(NULL)
			}
		
		}
	AllTS[[i]] <- the.wp.hwt
	AllPVal[[i]] <- the.wp.pval
	AllLV[[i]] <- the.wp.lv
	AllRecord[[i]] <- the.wp.record
	}
	if (silent==FALSE)
		cat("\n")

	alllv <- allindex <- alllitscale <- allbigscale <- allTS <- NULL

	for(i in J:(lowlev+1))	{
		for(k in (J-WTscale-2):0)	{

			if (length(AllTS[[i]]$d[[k+1]]) != length(AllRecord[[i]]$d[[k+1]]))	{
				s <- paste("i: ", i, "k: ", k, "l-AllTS:", length(AllTS[[i]]$d[[k+1]]), " l-AllRecord: ", length(AllRecord[[i]]$d[[k+1]]))
				stop(s)
			}


			if (!is.null(AllTS[[i]]$d[[k+1]]))	{
			  allTS <- c(allTS, AllTS[[i]]$d[[k+1]])
			  allindex <- c(allindex, AllRecord[[i]]$d[[k+1]])
			  lscale <- AllRecord[[i]]$c[[k+1]]
			  allbigscale <- c(allbigscale, rep(i, length(lscale)))
			  alllitscale <- c(alllitscale, lscale)
			  alllv <- c(alllv, AllLV[[i]]$c[[k+1]])	

			  if (verbose==TRUE)	{
				if (length(lscale) != length(AllLV[[i]]$c[[k+1]]))	{
				  if (!is.null(lscale))	{
					print(lscale)
					s <- paste("i: ", i, "k: ", k, "length(lscale):", length(lscale), " l-AllLV: ", length(AllLV[[i]]$c[[k+1]]))
					stop(s)
					}
				  else
					cat("Warning: lscale is NULL\n")
				 }
				}
				
			  }
		}
	}
allpvals.unadjust <- allpvals <- 2*pnorm(-abs(allTS))

allpvals <- p.adjust(allpvals, method=mc.method)

nreject <- sum(allpvals < alpha)

l <- list(nreject=nreject, mc.method=mc.method, AllTS=AllTS,
	AllPVal=AllPVal, alpha=alpha, x=x, xSD=xS.stat,
	allTS = allTS,
	allpvals=allpvals, allbigscale=allbigscale, alllitscale=alllitscale,
	allindex=allindex, alllv=alllv, allpvals.unadjust=allpvals.unadjust)
class(l) <- "tosANYN"
return(l)
}
