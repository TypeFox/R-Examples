##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright � 2008,2009 The University of Edinburgh                     #
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

# The R stub for the pcor function. This does some rudimentary
# argument type checking and then hands off to the C stub.

pcor <- function(
  data_x                       # input numerical matrix
, data_y = NULL                # matrix with compatible dimensions to data_x.
, distance   = FALSE           # Return the distance matrix instead of the correlation coefficients
, caching_   = "mmeachflush"   # getOption("ffcaching")
, filename_  = NULL            # tempfile(pattern = pattern, tmpdir = getOption("fftempdir"))

)
  {
    # we only work on doubles
    vmode_ = "double"

    # determine filename and finalizer
    # if user choses to work on a temporary file it will be deleted whenn all
    # references to the ff object are closed

    if (is.null(filename_)){
      # delete if temporary ff object
      finalizer_<- "delete"
    } else {
      finalizer_<- "close"
    }
     
    if (is.null(filename_)){
      # temporary ff object
      filename_<- tempfile(pattern =  "ff" , tmpdir = getwd())
    }

     # determine length of the result correlation matrix
    
    if (is.matrix(data_x) && is.numeric(data_x)) {
      height = dim(data_x)[2]
      length_ <- height * height
    }
    else
      stop(..sprintMsg$error["non.numeric"])

    if(!is.null(data_y)) {
      if (!is.matrix(data_y) && !is.numeric(data_y))
        stop(..sprintMsg$error["non.numeric"])
      if (dim(data_x)[0] != dim(data_y)[0] && dim(data_x)[1] != dim(data_y)[1])
        stop(..sprintMsg$error["no.dims"])
    }

    if (is.null(caching_))
      caching_<- getOption("ffcaching")
    else
      caching_<- match.arg(caching_, c("mmnoflush", "mmeachflush"))

    # Check the value of the "distance" option
    if ( (!is.logical(distance)) || (length(distance)>1) ) {
        warning(paste("Value of option \"distance\" must be a scalar logical (TRUE of FALSE). You supplied : ", distance))
        return(FALSE)
    }

    # Call C interface
    return_val <- .Call("pcor", data_x, data_y, filename_, distance)
	
	colnames1 <- colnames(data_x)
	if(is.null(data_y))
	  colnames2 <- colnames(data_x)
	else
	  colnames2 <- colnames(data_y)
	  
    # Return values from the interface have meaning.
    #  0    -->     success
    # -1    -->     MPI is not initialized
    # -2    -->     Only the master process exists, no workers
    if ( return_val == 0 ) {
      # Open result binary file and return as ff object
      result = ff(
        dim=c(height,height),
		dimnames=list(colnames1,colnames2),
        , filename=filename_
        , vmode=vmode_
        , caching=caching_
        , finalizer=finalizer_
        , length=length_
        )
    } else {

        if ( return_val == -1 )
            warning(paste("MPI is not initialized. Function is aborted.\n"))
        if ( return_val == -2 )
            warning(paste("No worker processes exist. Function pcor() is aborted.\n"))
        result <- FALSE
    }

    return(result)
}

