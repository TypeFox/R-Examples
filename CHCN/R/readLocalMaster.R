readLocalMaster <- function(){
  f <- system.file("external/EnvCanadaMaster.csv", package = "CHCN")
  return(read.csv(f, stringsAFactors = FALSE))
}