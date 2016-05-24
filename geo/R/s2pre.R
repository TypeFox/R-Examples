#' Writes out data.frame or matrix to a prelude-file.
#' 
#' Data.frame or matrix object is written to a prelude-file, that inherits
#' names/dimnames attributes from the object.
#' 
#' 
#' @param data Data.frame or matrix object.
#' @param file Name of the output file ("Splus.pre" by default).
#' @param na.replace A character to replace NA with in the output file ("" by
#' default).
#' @return A prelude-file representation of the data-object is written to a
#' file.
#' @section Side Effects: No warning is given if the filename "file" already
#' exists -- it is simply over-written.
#' @seealso \code{\link{cat}}, \code{\link{write}}.
#' @examples
#' 
#' \dontrun{Within Splus:
#'        > tmp.test.frame
#'           tolur1     tolur2 stafir1
#'         1     11 0.04625551       a
#'         2     12 0.04845815       a
#'         3     13 0.05066079      NA
#'         4     14 0.05286344       a
#'         5     15 0.05506608       a
#'         6     16 0.05726872       b
#'         7     17 0.05947137       b
#'         8     18         NA       b
#'         9     19 0.06387665       b
#'        10     20 0.06607930       b
#'        > s2pre(tmp.test.frame,file="prufa.pre",na.replace="-1")
#'        >
#' 
#'        From UNIX:
#' 
#'        hafbitur/home/reikn/gardar/Papers/Methods95 [435] cat prufa.pre
#'        linu_nofn       tolur1  tolur2  stafir1
#'        ---------       ------  ------  -------
#'        1       11      0.04625551      a
#'        2       12      0.04845815      a
#'        3       13      0.05066079      -1
#'        4       14      0.05286344      a
#'        5       15      0.05506608      a
#'        6       16      0.05726872      b
#'        7       17      0.05947137      b
#'        8       18      -1      b
#'        9       19      0.06387665      b
#'        10      20      0.06607930      b
#'        hafbitur/home/reikn/gardar/Papers/Methods95 [436]
#' }
#' @export s2pre
s2pre <-
function(data, file = "splus.pre", na.replace = "")
{
	# data       :matrix or data.frame.
	# na.replace :a character to replace NA with.
	#
	# VALUE      :a prelude file, named "Splus.pre" by default.
	if(is.data.frame(data)) data <- as.matrix.data.frame(data)
	data[is.na(data) | data == "NA"] <- na.replace
	col.names <- dimnames(data)[[2]]
	if(is.null(col.names) || length(col.names) == 0)
		col.names <- paste("dalkur", 1:ncol(data), sep = "")
	row.names <- dimnames(data)[[1]]
	if(!is.null(row.names) && length(row.names) > 0) {
		col.names <- c("linu_nofn", col.names)
		data <- cbind(row.names, data)
	}
	n.of.col <- length(col.names)
	# Write out rownames:
	cat(col.names, sep = c(rep("\t", n.of.col - 1), "\n"), file = file)
	strika.lina <- rep("", n.of.col)
	for(i in 1:n.of.col)
		strika.lina[i] <- paste(rep("-", nchar(col.names[i])), collapse
			 = "")
	# Write out the ------ line:
	cat(strika.lina, sep = c(rep("\t", n.of.col - 1), "\n"), file = file,
		append = T)
	# Write out the data:
	cat(t(data), sep = c(rep("\t", n.of.col - 1), "\n"), file = file, 
		append = T)
	return(invisible(NULL))
}

