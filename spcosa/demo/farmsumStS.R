# example spcosa package: stratified simple random sampling

# check if required packages are available
if (suppressWarnings(!require(rgdal))) {
    stop("This demo requires package 'rgdal'.\nThis package is currently not available. Please install 'rgdal' first.", call. = FALSE)
}    

# initialize pseudo random number generator
set.seed(700124)

# read vector representation of the Farmsum paddock
shpFarmsum <- readOGR(dsn = system.file("maps", package = "spcosa"), layer = "farmsum")

# stratify Farmsum into 50 strata
myStratification <- stratify(shpFarmsum, nStrata = 50)

# plot stratification
plot(myStratification)

# sample two sampling units per stratum
mySamplingPattern <- spsample(myStratification, n = 2)

# plot sampling pattern
plot(myStratification, mySamplingPattern)

# extract sampling points
myData <- as(mySamplingPattern, "data.frame")

# simulate data (in real world cases these data have to be obtained by field work)
myData$observation <- rnorm(n = nrow(myData), mean = 10, sd = 1)

# design-based inference
estimate("spatial mean", myStratification, mySamplingPattern, myData["observation"])
estimate("standard error", myStratification, mySamplingPattern, myData["observation"])
