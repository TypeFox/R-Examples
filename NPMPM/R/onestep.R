onestep <-
function(cfu_in,pstep,pdistpsteps,apsomvmethod,intexmethod){

# DRAW RANDOM VALUES FOR THE PARAMETERS OF THE PSTEP
pstepparameters <- pdistpsteps(pstep)
psteptime<- pstepparameters[1]

###################################
# GET THE FACTOR FOR DILUTION OR CONCENTRATION
factor <- pstep$factor
###################################

# GET ALL APPROPRIATE SERIES OF MEASURED VALUES
mvalues <- apsomvmethod(pstep)
###################################
countme <- length(mvalues)
if (countme<30) {
countwarnings <- countwarnings+1
} # END IF


# VECTOR FOR THE CALCULATED CFUs
cfucalc <- numeric()

# THE LOOP FOR ALL APPROPRIATE SERIES OF MEASURED VALUES
seriesnumber <- 1

repeat {if (seriesnumber>countme) break else {

series <- mvalues[seriesnumber]

lasttime <- max(measured_values[(measured_values$raw_data_key==as.vector(series)),]$time)

# CALCULATION OF THE LOG CFU 
if (lasttime<psteptime) {
cfucalc[seriesnumber] <- NA
warning("No Extrapolation")
} else {
if (cfu_in==0) {
cfucalc[seriesnumber] <- -10^(35)
} else { 
cfucalc[seriesnumber] <- intexmethod(psteptime, series, cfu_in)
}} # END ELSE,ELSE
seriesnumber <- seriesnumber+1
}} # END REPEAT

# ESTIMATORS FOR THE PARAMETERS OF THE LOGNORMAL DISTRIBUTION
# !!! FROM THE NUMBERS, NOT FROM THE LOG-VALUES !!!
meandist <- mean(cfucalc,na.rm=TRUE)
sddist <- sqrt(var(cfucalc,na.rm=TRUE))

# STOPP CRITERION - IF INTERPOLATION LEADS TO NA FOR ALL SOMV: STOPP
if (is.na(meandist)) {
stop("Simulation stopped! There are not enough data for interpolation with",
 intexmethod ," in process step ",psteps$stepnumber," ! ")
} # END IF

# DRAW ONE RANDOM VALUE, 
logcfu <- rnorm(1,meandist,sddist)

#  
cfu_out <- 10^{logcfu}

# "TRUNCATE" THE CFU AND STORE IT IN THE VECTOR cfu
if (cfu_out<1) {cfu_out <- 0} 
if (cfu_out>10^9) {cfu_out <- 10^9} 
cfu_out <- cfu_out*factor

return(cfu_out)

} # END FUNCTION

