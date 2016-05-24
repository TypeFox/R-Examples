
.onAttach <- function(...) { 

  library(help=bgeva)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("\nThis is bgeva ",version,".\nFor overview type 'help(\"bgeva-package\")'.\n",sep="")
  packageStartupMessage(hello)
  
}






