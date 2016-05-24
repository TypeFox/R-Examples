EWSsmoothRM <-
function (S, s) 
{
#
# Apply running mean smoothing to spectrum S with bandwidth s
#
#
# Prepare smoothing filter
#
the.filter <- rep(1/(2 * s + 1), length = (2 * s + 1))
#
# Work out dimensions of spectrum
#
J <- S$nlevels
TT <- 2^J
#
# Now smooth each level of the spectrum one at a time
#
for (j in 0:(J - 1)) {
	#
	# Get the spectrum at level j
	#
        sl <- accessD(S, level = j)
	#
	# Apply smoothing filter
	#
        sl2 <- filter(x = sl, filter = the.filter)
	#
	# Get rid of NAs and augment with the value of the first answer	
	# to remake the series of the correct length
	#
        sl2 <- sl2[!is.na(sl2)]
        thetimesarg <- TT - length(sl2)
        sl2 <- c(rep(sl2[1], TT - length(sl2)), sl2)
	#
	# Store the answer
	#
        S <- putD(S, level = j, v = sl2)
    }
#
# Return the smoothed spectrum
#
return(S)
}
