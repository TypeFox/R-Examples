siarsolomultigrouprun <-
function(siardata) {
# This function runs the multi group MCMC for siar

if(siardata$SHOULDRUN==FALSE || siardata$GRAPHSONLY ==TRUE) {
    cat("You must load in some data first (via option 1) in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    return(siardata)
}

cat("Run SIAR for multiple groups with a single observation in each group. \n")
cat("For this you will need to have input successfully some data in option 1. \n")
cat("In this instance, the target isotope file must be a three column file with \n")
cat("a group label in the first column. \n \n")
cat("See the demo for more details on this function. \n")
cat("\n")
cat("Press <Enter> to continue...")
readline()
invisible()

if(siardata$numgroups > 1) {

# Run size
runchoices <- c("Standard","Long","Very long")
runtitle <- "Choose the size of the model run:"
BADRUN <- TRUE
while(BADRUN ==TRUE) {
    runchoose <- menu(runchoices,title = runtitle)
    if(any(runchoose==seq(1,3))) BADRUN <- FALSE
}

# Now run the code
if(runchoose == 1) {
    siardata$iterations <- 200000   
    siardata$burnin <- 50000
    siardata$howmany <- 10000
    siardata$thinby <- 15
}
if(runchoose == 2) {
    siardata$iterations <- 400000   
    siardata$burnin <- 200000
    siardata$howmany <- 10000
    siardata$thinby <- 100
}
if(runchoose == 3) {
    siardata$iterations <- 1000000   
    siardata$burnin <- 400000
    siardata$howmany <- 20000
    siardata$thinby <- 300
}
siardata <- siarsolomcmcv4(data=siardata$targets,sources=siardata$sources,corrections=siardata$corrections,concdep=siardata$concdep,siardata=siardata)

return(siardata)

} else {

cat("This data has only 1 group - choose the single group option instead. \n \n")
return(list(targets=siardata$targets,sources=siardata$sources,corrections=siardata$corrections,concdep=siardata$concdep,PATH=siardata$PATH,TITLE=siardata$TITLE,numgroups=siardata$numgroups,numdata=siardata$numdata,numsources=siardata$numsources,numiso=siardata$numiso,SHOULDRUN=siardata$SHOULDRUN,GRAPHSONLY=siardata$GRAPHSONLY,EXIT=siardata$EXIT,output=siardata$output))

}

}
