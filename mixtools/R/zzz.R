.onAttach <- function(lib, pkg){
  info <- packageDescription("mixtools")
  packageStartupMessage(
    paste('mixtools package, version ', info$Version, ', Released ', info$Date, '\n',
          'This package is based upon work supported by the National Science ',
          'Foundation under Grant No. SES-0518772.\n', sep="")
 )
}


