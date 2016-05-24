writeData <- function(Data,filename = "EnvCanadaData.dat",directory=DATA.DIRECTORY){
  
  if (!file.exists(directory)) dir.create(directory)
  fname <- file.path(directory,filename,fsep =.Platform$file.sep)
   
  write.csv(Data,fname)
}