#####
## Collection of internal functions used
##
#####
.onAttach <- function(...)
   {
   packageStartupMessage(paste(" **********   MAMS Version",
packageDescription("MAMS")$Version), "********** \n")
   packageStartupMessage("Type MAMSNews() to see new features/changes/bug fixes.\n")
   #if(as.numeric(R.Version()$major) >= 2 & as.numeric(R.Version()$minor) < 10) {packageStartupMessage(paste("The help functions for this package might not work properly. Please upgrade R to Version 2.10 or above to fix this problem.\n"))}
   }

