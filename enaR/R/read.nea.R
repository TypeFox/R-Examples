#' read.nea.RData
#' INPUT = Model Data (flows, inputs, outputs, storage) formatted as for NEA.m, saved as CSV file
#'        S=  |[F][z][X]|
#'            |[y][0][0]|
#' OUPUT = R Network data object for use with enaR
#'
#' Borrett | July 15, 2013
#' --------------------------------------------------

read.nea <- function(file="file name",sep=',',warn=TRUE){
  dat <- read.table(file,header=FALSE,sep=sep)  # assumes 
  n <- max(dim(dat)) - 2
  Flow <- t(dat[1:n,1:n])   # NEA.m stores flows col to row, so here we transpose
  z <- dat[1:n,(n+1)]  # inputs
  y <- dat[(n+1),1:n]  # outputs
  X <- dat[1:n,(n+2)]  # storage
  if (warn){
    model <- pack(flow=Flow,input=z,respiration=y,storage=X)  # create network data object
  }else{
    suppressWarnings(model <- pack(flow=Flow,input=z,respiration=y,storage=X))   # create network data object
  }
  return(model)
}
