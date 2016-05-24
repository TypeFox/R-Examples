tdir <- system.file("testfiles",package="inlinedocs")
testfiles <- Sys.glob(file.path(tdir,"*.R"))
library(inlinedocs)
for(f in testfiles){
  print(f)
  test.file(f,verbose=FALSE)
}

