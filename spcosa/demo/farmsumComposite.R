# example spcosa package: stratified simple random sampling

# check if required packages are available
if (suppressWarnings(!require(rgdal))) {
    stop("This demo requires package 'rgdal'.\nThis package is currently not available. Please install 'rgdal' first.", call. = FALSE)
}    

# initialize pseudo random number generator
set.seed(700124)

# read vector representation of the Farmsum paddock
shpFarmsum <- readOGR(dsn = system.file("maps", package = "spcosa"), layer = "farmsum")

# stratify Farmsum into 40 strata of equal size
# try 10 random starting configurations (*very* slow!)
myStratification <- stratify(shpFarmsum, nStrata = 40, equalArea = TRUE,
    nGridCells = 5000, nTry = 5, verbose = TRUE)

# plot stratification
plot(myStratification)

# sample two sampling units per stratum
mySamplingPattern <- spsample(myStratification, n = 2, type = "composite")

# plot sampling pattern
plot(myStratification, mySamplingPattern)

# simulate data (in real world cases these data have to be obtained by field work)
myData <- data.frame(observation = rnorm(n = 2, mean = 10, sd = 1))

# design-based inference
estimate("spatial mean",   myStratification, mySamplingPattern, myData)
estimate("standard error", myStratification, mySamplingPattern, myData)

