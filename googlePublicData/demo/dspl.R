################################################################################
# rdspl demo
# In this demo, we'll use the Public Data Explorer tutorial data available at:
#   https://developers.google.com/public-data/docs/tutorial
#
################################################################################
pause <- function() {  
  invisible(readline("\nPress <return> to continue: ")) 
}

# Complete directory path where the data is saved (in this case, where the
# package is installed)
data.path <-try(paste(.libPaths()[1],'/googlePublicData/data',sep=''), silent=T)
if (class(data.path) == "try-error") data.path <-paste(.libPaths()[2],'/googlePublicData/data',sep='')
pause()

# First Simplest example, 
mydspl <- dspl(path=data.path, sep=";")

# Printing the data
pause()
mydspl

# Summary of the dspl class object
pause()
summary(mydspl)

# End of the demo