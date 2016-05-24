
# this function is internal, not meant for the user


.onAttach <- function(...)  {
  library(help=lm.br)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste( " lm.br  version ", version, 
    ",  '?lm.br' starts help", sep="" )
  packageStartupMessage( hello )
}

