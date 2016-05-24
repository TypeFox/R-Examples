# TODO: Add comment
# 
# Author: Giorgio
###############################################################################

#adding a startup message
#for future version, should test : if(verbose <- getOption("verbose")) 

.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage('Package:  ', desc$Package, '\n',
                        'Version:  ', desc$Version, '\n', 
                        'Date:     ', desc$Date, '\n',
                        'BugReport: ', desc$BugReports, '\n\n')
}




# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("mbbefd", libpath)
}