.onAttach <- function(...) { 

  library(help=SemiParSampleSel)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("\nThis is SemiParSampleSel ",version,".\nFor overview type 'help(\"SemiParSampleSel-package\")'.\n",sep="")
  packageStartupMessage(hello)
  
}






