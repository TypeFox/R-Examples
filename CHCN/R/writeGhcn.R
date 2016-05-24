writeGhcn <- function(data,directory = DATA.DIRECTORY,filename = "TaveCHCN.dat"){
  fname <- file.path(directory,filename,fsep = .Platform$file.sep)
  write.table(data,fname)
  
}