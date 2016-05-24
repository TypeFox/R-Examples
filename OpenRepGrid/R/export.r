#' Save grid in a text file (txt).
#'
#' \code{saveAsTxt} will save the grid as a \code{.txt} file
#' in format used by \pkg{OpenRepGrid}. This file format can also 
#' easily be edited by hand (see \code{\link{importTxt}} for a 
#' description).
#' The funtion will open an interactive dialog box to let the user 
#' enter a filename if no \code{file} argument is supplied
#' in the function call. 
#'
#' @param x     \code{repgrid} object.
#' @param file  Filename to save the grid to. The name should have 
#'              the suffix .txt. If the function is called without
#'              specifying this argumnet a dialog box is opened.
#' @return      Invisibly returns the name of the file.
#'
#' @note
#' Structure of a txt file that can be read by \code{\link{importTxt}}.
#'
#' \code{---------------- .txt file -----------------}
#'
#' \code{anything not contained within the tags will be discarded}
#'
#' \tabular{l}{
#' \code{ELEMENTS}         \cr
#' \code{element 1}        \cr
#' \code{element 2}        \cr
#' \code{element 3}        \cr
#' \code{END ELEMENTS}     \cr
#' \cr
#' \code{CONSTRUCTS}                 \cr
#' \code{left pole 1 : right pole 1} \cr
#' \code{left pole 2 : right pole 2} \cr
#' \code{left pole 3 : right pole 3} \cr
#' \code{left pole 4 : right pole 4} \cr
#' \code{END CONSTRUCTS}             \cr
#' \cr
#' \code{RATINGS}        \cr
#' \code{1 3 2}          \cr
#' \code{4 1 1}          \cr
#' \code{1 4 4}          \cr
#' \code{3 1 1}          \cr
#' \code{END RATINGS}    \cr
#' \cr
#' \code{RANGE}          \cr
#' \code{1 4}            \cr
#' \code{END RANGE}      \cr
#' }
#' \code{---------------- end of file ----------------}
#'
#' @author    Mark Heckmann
#' @export
#'
#' @seealso     \code{\link{importTxt}}
#'
#' @examples \dontrun{
#' 
#'  x <- randomGrid()
#'  saveToTxt(x, "random.txt")
#'
#' }
#'
saveAsTxt <- function(x, file=NA){  
	if (is.na(file)){
	  #require(tcltk)
  	fileName <- tclvalue(tkgetSaveFile())     # open filename dialog  
	} else 
	  fileName <- file
	
	enames <- getElementNames(x)
	cnames <- getConstructNames(x)
	scores <- getRatingLayer(x)
	# write txt file
	con <- file(fileName, "w")                # open an output file connection
	
	cat("=========================\n", file = con)
	cat("Data File for OpenRepGrid\n", file = con)
	cat("=========================\n", file = con)
	
	# write element names
	cat("\nELEMENTS\n", file = con)
	for (ename in enames)
	  cat(ename, "\n", sep="", file=con)
	cat("END ELEMENTS\n", file = con)
	
	# write construct names
	cat("\nCONSTRUCTS\n", file = con)
	cnames.string <- paste(cnames[,1], ":", cnames[, 2])
	for (cname in cnames.string)
	  cat(cname, "\n", sep="", file=con)
	cat("END CONSTRUCTS\n", file = con)
	
	# write ratings to file
	cat("\nRATINGS\n", file = con)
	write.table(scores, file = con, sep = " ", 
	            row.names = FALSE, col.names=FALSE)      
	cat("END RATINGS\n", file = con)
  
  # write scale range to file
  cat("\nRANGE\n", file = con)
  cat(x@scale$min, x@scale$max, "\n", file=con)
	cat("END RANGE\n", file = con)
	
	close(con)
	cat("grid succesfully written to file: ", unlist(fileName))
	
	invisible(fileName)
}


### .TXT FILE FORMAT ###

# everything not contained within the tags will be discarded
# 
# ELEMENTS
# element 1
# element 2
# element 3
# element 4
# END ELEMENTS
# 
# CONSTRUCTS
# left pole 1 : right pole 1
# left pole 2 : right pole 2
# left pole 3 : right pole 3
# left pole 4 : right pole 4
# END CONSTRUCTS
# 
# RATINGS
# 1 3 NA 2
# 4 1 3 NA
# 1 4 1 3
# 3 1 1 6
# END RATINGS
# 
# RANGE
# 1 6
# END RANGE


# exportTxt <- function(x, file=NULL)
# {
#   file <- "test.txt"
#   # redirect output to connection
#   sink(file)
#   
#   # elements
#   cat("\n")
#   cat("ELEMENTS\n")
#   for (name in getElementNames(g))
#     cat(name, "\n")
#   cat("END ELEMENTS\n\n")
#   
#   # constructs
#   cat("CONSTRUCTS\n")
#   for (name in getConstructNames2(g, sep=" : ", trim=NA))
#     cat(name, "\n")
#   cat("END CONSTRUCTS\n\n")
#   
#   # ratings
#   cat("RATINGS\n")
#   r <- getRatingLayer(g)
#   for (i in 1:nrow(r))
#     cat(r[i, ], "\n")
#   cat("END RATINGS\n\n")
#   
#   # range
#   cat("RANGE\n")
#   cat(getScale(g), "\n")  
#   cat("END RANGE\n\n")
#   
#   # reset output stream
#   sink()
# }


