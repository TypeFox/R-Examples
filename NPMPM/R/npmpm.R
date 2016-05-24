npmpm <-
function(inoculum,errorb=0.01,intexmethod=cfuinterpolation,
apsomvmethod=apsomv,pdistpsteps=psteppar,lastiteration=numberiterations,calconestep=onestep){
errorb<<-errorb
intexmethod <<- intexmethod
apsomvmethod <<- apsomv
pdistpsteps <<- psteppar
lastiteration <<- numberiterations
calconestep <<- onestep

inputsettings <<- match.call()

# CHECK IF THE INOCULUM IS NOT NEGATIVE  
if (min(inoculum)<0) {
stop("inoculum must be nonnegative")}

if (length(inoculum[inoculum>0])==0) {
warning("no contamination")
stop(return(0))}

# ELSE GO ON

# CHECK IF THE INPUT MATCHES THE DATA.FRAMES
# THE NUMBER OF INOCULA SHOULD EQUAL THE NUMBER OF PROCESS STEPS
pnumber <- length(psteps$id)
if (length(inoculum) != pnumber)
stop("The number of inocula does not match the number of process steps.")

# ELSE GO ON

# COUNTER FOR NUMBER OF ITERATIONS
iterations <<- 0

# NUMBER OF THE FIRST STEP IN WHICH A CONTAMINATION OCCURS
framestep <- data.frame(a=c(1:length(inoculum)),b=inoculum)
minstepnumber <- min(framestep[framestep$b>0]$a)
maxstepnumber <- max(framestep[framestep$b>0]$a)

# VECTOR WITH THE END-CFUs
cfu <<- c(0)

# VECTOR WITH THE AVERAGES CALCULATED FROM THE END-CFUs
cfu_average <- c(0)

# NUMBER OF PSTEPS WHITH LESS THAN 30 SERIES OF MEASURED VALUES FOR CALCULATION
countwarnings <<- 0

#######################################################
# START THE ITERATIONS
#######################################################

# OUTER LOOP
while (lastiteration(iterations,errorb,cfu_average))
{
# NUMBER OF THE CURRENT PROCESS STEP, BEGIN WITH THE FIRST STEP IN WHICH 
# A CONTAMINATION OCCURS
stepnr <- minstepnumber-1

# COUNT THE ITERATIONS
iterations <- iterations+1
print(iterations)

# THE CONTAMINATION IN THE PREVIOUS STEP (NO CONTAMINATION!)
cfu_in <- 0

##################
# LOOP FOR THE PROCESS STEPS (IN ONE ITERATION)
# INNER LOOP
##################
repeat {if (stepnr==pnumber) break else {

# COUNT THE PROCESS STEPS CALCULATED IN THE CURRENT ITERATION
stepnr <- stepnr+1

# IF THE CONTAMINATION FROM THE PREVIOUS STEP IS ZERO AND THERE WILL BE
# NO MORE CONTAMINATION IN A FUTURE STEP THE CFU REMAINS ZERO
if (cfu_in==0 & (stepnr>maxstepnumber)) {cfu_in <- 0} else {
#
# ELSE TAKE THE CONTAMINATION FROM THE PREVIOUS STEP AND ADD THE NEW INOCULUM
cfu_in <- cfu_in + inoculum[stepnr]

#
# CURRENT PROCESS STEP
pstep <- psteps[(psteps$stepnumber==stepnr),]


# CALCULATE THE CFU AT THE END OF THE PROCESS STEP
cfu_in <- calconestep(cfu_in,pstep,pdistpsteps,apsomvmethod,intexmethod)
} # END ELSE

}} # END ELSE + REPEAT (PROCESS STEPS IN ONE ITERATION)

# SAVE THE CALCULATED CFU FROM THIS ITERATION
cfu[iterations] <<- cfu_in

# CALCULATE THE MEAN OF THE HITHERTO CALCULATED CFUs
cfu_average[iterations] <- mean(cfu)


} # END WHILE (OUTER LOOP)

return(cfu)

} # END FUNCTION

