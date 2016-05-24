.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste(pkgname, RFver))
  #packageStartupMessage("Ensemble Learning and Integration")
  packageStartupMessage("Heart and Lung Institute, Imperial College London &\nScientific Computing Group, Sentrana Inc.")
  #packageStartupMessage("Loading packages: gbm,nnet,e1071,randomForest,kknn,glmnet,doParallel")
}

.onLoad <- function(libname, pkgname) {
  #suppressMessages(library(gbm))
  #suppressMessages(library(nnet))
  #suppressMessages(library(e1071))
  #suppressMessages(library(randomForest))
  #suppressMessages(library(kknn))
  #suppressMessages(library(glmnet))
  #suppressMessages(library(doParallel))
}

