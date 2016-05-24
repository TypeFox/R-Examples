# Build against openssl libraries that were compiled with the Rtools gcc toolchain.
if(!file.exists("../windows/webp-0.5.0/include/webp/encode.h")){
  if(getRversion() < "3.3.0") setInternet2()
  download.file("https://github.com/rwinlib/webp/archive/v0.5.0.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
