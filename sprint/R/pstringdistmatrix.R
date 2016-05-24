##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
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

# Usage:

# stringdistmatrix(a, b,
# method = "h",
# weight = c(d = 1, i = 1, s = 1, t = 1), maxDist = 0,
# ncores = 1)

# Arguments:

# a: R object (target); will be converted by ‘as.character’.

# b: R object (source); will be converted by ‘as.character’.

# method: Method for distance calculation (see details)

# weight: The penalty for deletion, insertion, substitution and
# transposition, in that order.  Weights must be positive and
# not exceed 1. ‘weight[4]’ is ignored when ‘method='lv'’ and
# ‘weight’ is ignored completely when ‘method='h'’.

# maxDist: Maximum string distance before calculation is stopped,
# ‘maxDist=0’ means calculation goes on untill the distance is
# computed.

# ncores: number of cores to use. Parallelisation is over ‘b’, so the
# speed gain by parallelisation is highest when ‘b’ is shorter
# than ‘a’.

# filename: 

#(x, method="hamming", filename="output_file")
pstringdistmatrix <- function (a, b, method="h", filename=NULL, weight=NULL, maxDist=0, ncores=NULL) {
	
# TODO use a and b.
	if(!(identical(as.character(a),as.character(b)))){
		stop(..sprintMsg$error["not.supported.diff.strings"])
	}
	
	data <- as.character(a)
	
	if(!method=="h"){
		stop(..sprintMsg$error["not.supported.non.hamming"])	}
	
# ncores should not be set.
	if(!(is.null(ncores)||(1 == ncores))){
		stop(..sprintMsg$error["not.supported.ncores"])
	}

	if(!(maxDist == 0)){
		stop(..sprintMsg$error["not.supported.maxDist"])
	}
	   
# determine filename and finalizer
# if user choses to work on a temporary file it will be deleted whenn all
# references to the ff object are closed
	
    if (is.null(filename)){
# delete if temporary ff object
		finalizer_<- "delete"
    } else {
		finalizer_<- "close"
    }
	
    if (is.null(filename)){
# temporary ff object
		filename <- tempfile(pattern =  "ff" , tmpdir = getwd())
    }
	
    
  if(!length(data)) stop(..sprintMsg$error["empty"])
  
	flatData <- paste(data, collapse = '')

	  dataNames <- names(data)


	sample_width <- nchar(data[1]) #NB. Expect each string to be the same length. Should add a check for this.
  number_of_samples <- length(data)
	
	if(!exists("dataNames")||is.null(dataNames)){
		dimnames_ <- NULL
	}
	
  if(sample_width<1 || number_of_samples<2) stop(..sprintMsg$error["empty"])

  return_val <- .C("pstringDist",
                   as.character(flatData),
                   as.character(filename),
                   as.integer(sample_width),
                   n=as.integer(number_of_samples)                   
                   )

	# The number_of_samples is overloaded to also indicate whether MPI is initialized.
	# -1    -->     MPI is not initialized
	
	vmode_ <- "integer"
	caching_ <- "mmeachflush"
	filename_ <- as.character(filename)
	if(!exists("dimnames_")){dimnames_ <- list(dataNames,dataNames)}
	
		if ( return_val$n == -1 )  {
			warning(paste("MPI is not initialized. Function is aborted.\n"))
			result <- FALSE
		} else {
			# Open result binary file and return as ff object
		result = 	ff(
		   dim=c(number_of_samples,number_of_samples)
		   , dimnames=dimnames_
		   , filename=filename_
		   , vmode=vmode_
		   , caching=caching_
		   , finalizer=finalizer_
		   , length=(number_of_samples*number_of_samples)
		   )
    } 
  return(result)
}
