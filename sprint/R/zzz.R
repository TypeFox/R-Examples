##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################


# These functions are called automatically by the R extensions
# system. They perform some R-side initialisation, then call
# down to the library initialisation.

## Called after the library has finished loading, and all NAMESPACE
## exports have been processed
.onAttach <- function(lib, pkg) {
  ## We start the worker after attaching the library so that
  ## lazy-loaded R code in SPRINT is available to the slave processes
  invisible(return_val <- .Call("worker"))

  ## If the return value is not 0, which indicates a 'worker' MPI
  ## process, then immediately quit. This will prevent any messages
  ## printed out, or any other unwanted behavior.
  if(return_val != 0) {
    q(save = "no")
  }

  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  ver <- as.character(ver)
  packageStartupMessage("SPRINT ", ver, " loaded\n")
  packageStartupMessage(paste("Welcome to SPRINT\n",
                              "Please help us fund SPRINT by filling in \n",
                              "the form at http://www.r-sprint.org/ \n",
                              "or emailing us at sprint@ed.ac.uk and letting \n",
                              "us know whether you use SPRINT for commercial \n",
                              "or academic use."))
}

.onUnload <- function(libpath) {
	.Call("sprint_shutdown")
}

## Called when the extension is unloaded. This is expected to happen
## when the script terminates and R shuts down.
.Last.lib <- function(libpath) {
  .Call("sprint_shutdown")
}

