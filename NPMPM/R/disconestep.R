disconestep <-
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

countme <- length(mvalues)

##########
# CALCULATION OF ONE CFU FROM A RANDOM DRAWN SOMV
# LOOP: REPEAT UNTIL ONE SOMV IS CALCULATED (BEWARE NA)
##########
cfucalc <- NA

repeat{if(is.na(cfucalc)==FALSE)break else{

# draw random value in (1,countme)
seriesnumber <- ceiling(countme*runif(1,0,1))

series <- mvalues[seriesnumber]

lasttime <- max(measured_values[(measured_values$raw_data_key==as.vector(series)),]$time)

# CALCULATION OF THE LOG CFU 
if (lasttime<psteptime) {
cfucalc <- NA
warning("No Extrapolation")
} else {
if (cfu_in==0) {
cfucalc <- -10^(35)
} else { 
cfucalc <- intexmethod(psteptime, series, cfu_in)
}} # END ELSE,ELSE

# IF NO CFU WAS CALCULATED DISCARD THE SOMV AND TRY ANOTHER
if(is.na(cfucalc)) {
mvalues <- mvalues[-seriesnumber]
countme <- countme-1
} # END IF

# STOPP CRITERION - IF INTERPOLATION LEADS TO NA FOR ALL SOMV: STOPP
if (countme==0 & is.na(cfucalc)) {
stop("Simulation stopped! There are not enough data for interpolation with",
 intexmethod ," in process step ",psteps$stepnumber," ! ")
} # END IF

}} # END REPEAT, ELSE

# COUNT IF THERE WERE REALLY MORE THAN 30 APPROPRIATE SOMV 
if (countme<30) {
countwarnings <- countwarnings+1
} # END IF

#  
cfu_out <- 10^{cfucalc}

# "TRUNCATE" THE CFU AND STORE IT IN THE VECTOR cfu
if (cfu_out<1) {cfu_out <- 0} 
cfu_out <- cfu_out*factor

return(cfu_out)

} # END FUNCTION

