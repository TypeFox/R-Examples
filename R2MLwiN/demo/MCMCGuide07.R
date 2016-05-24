############################################################################
#     MLwiN MCMC Manual
#
# 7   Using the WinBUGS Interface in MLwiN . . . . . . . . . . . . . . . .83
#
#     Browne, W.J. (2009) MCMC Estimation in MLwiN, v2.13. Centre for
#     Multilevel Modelling, University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

# 7.1 Variance components models in WinBUGS . . . . . . . . . . . . . . . 84

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

## Read tutorial data
data(tutorial, package = "R2MLwiN")

## openbugs executable
if (!exists("openbugs")) openbugs <- "C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
while (!file.access(openbugs, mode = 0) == 0 || !file.access(openbugs, mode = 1) == 0 || !file.access(openbugs, mode = 4) == 
  0) {
  cat("Please specify the path for the OpenBUGS executable:\n")
  openbugs <- scan(what = character(0), sep = "\n")
  openbugs <- gsub("\\", "/", openbugs, fixed = TRUE)
}

# User's input if necessary

## winbugs executable winbugs <- 'C:/Program Files (x86)/WinBUGS14/WinBUGS14.exe'

## The highest level comes first, then the second highest and so on
## Uses the results from IGLS to create initial values for bugs
## Fit the model by calling openbugs using the rbugs package
mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, show.file = TRUE), 
  BUGO = c(version = 4, n.chains = 1, debug = FALSE, seed = 1, bugs = openbugs, OpenBugs = TRUE), data = tutorial)

summary(mymodel1)
summary(mymodel1[, "beta[2]"])
sixway(mymodel1[, "beta[2]", drop = FALSE])

# 7.2 So why have a WinBUGS interface ? . . . . . . . . . . . . . . . . . 92
# 7.3 t distributed school residuals . . . . . . . . . . . . . . . . . . .92

## Download the model, initial, data files
modelfile <- paste0(tempdir(), "/tutorial1_model.txt")
download.file("http://www.bristol.ac.uk/cmm/media/r2mlwin/tutorial1_model.txt", modelfile, method = "auto")
file.show(modelfile)

initfile <- paste0(tempdir(), "/tutorial1_inits.txt")
download.file("http://www.bristol.ac.uk/cmm/media/r2mlwin/tutorial1_inits.txt", initfile, method = "auto")
file.show(initfile)

datafile <- paste0(tempdir(), "/tutorial1_data.txt")
download.file("http://www.bristol.ac.uk/cmm/media/r2mlwin/tutorial1_data.txt", datafile, method = "auto")

bugEst <- paste0(tempdir(), "/tutorial1_log.txt")


chains.bugs1 <- mlwin2bugs(D = "t", levID = c("school", "student"), datafile, initfile, modelfile, bugEst, fact = NULL, 
  addmore = NULL, n.chains = 1, n.iter = 5500, n.burnin = 500, n.thin = 1, debug = TRUE, bugs = openbugs, bugsWorkingDir = tempdir(), 
  OpenBugs = TRUE)
## Close winbugs manually
summary(chains.bugs1)
sixway(chains.bugs1[, "df", drop = FALSE])

chains.bugs2 <- mlwin2bugs(D = "t", levID = c("school", "student"), datafile, initfile, modelfile, bugEst, fact = NULL, 
  addmore = NULL, n.chains = 1, n.iter = 12000, n.burnin = 2000, n.thin = 1, debug = TRUE, bugs = openbugs, bugsWorkingDir = tempdir(), 
  OpenBugs = TRUE)
## Close winbugs manually
summary(chains.bugs2)
sixway(chains.bugs2[, "df", drop = FALSE])

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 96





############################################################################
