# Build against mingw-w64 build of redland
if(!file.exists("../windows/redland-1.0.17/include/redland.h")){
  if(getRversion() < "3.3.0") setInternet2()
  download.file("https://github.com/rwinlib/redland/archive/v1.0.17.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
