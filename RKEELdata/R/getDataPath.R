#Returns datasets path
getDataPath <- function(){

  dataPath = system.file("datasets", "", package="RKEELdata")

  return(dataPath)
}
