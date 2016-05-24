
.onAttach <- function(...) { 

  library(help=SemiParBIVProbit)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("\nThis is SemiParBIVProbit ",version,".\nFor overview type 'help(\"SemiParBIVProbit-package\")'.\n",sep="")
  packageStartupMessage(hello)
  
}






