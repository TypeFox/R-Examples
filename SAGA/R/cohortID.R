cohortID <- function(sexed=T){
  Cmatrix <- read.csv(file = system.file("cmatrix.xy.csv", package = "SAGA"), row.names=1)[, -1]
  ids <- cbind(1:72,row.names(Cmatrix))
  if(sexed==F){
    ids <- ids[seq(from=3, to=72, by=3),]
  }
  return(ids)
}