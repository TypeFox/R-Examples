work.dir <-
function(dir.name){
  WorkinDir=paste(getwd(),"/", dir.name, "/",sep="")
  dir.create(WorkinDir)
  file = list.files()
  file.copy(file, WorkinDir)
  setwd(WorkinDir)
  unlink("muma.R")
}
