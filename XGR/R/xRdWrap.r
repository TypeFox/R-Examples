#' Function to wrap texts from Rd files
#'
#' \code{xRdWrap} is supposed to wrap texts from Rd files under a given directory.
#'
#' @param path a directory containing Rd files
#' @param remove.dontrun logical to indicate whether to remove the restriction of not running examples. By default, it sets to FALSE without any modefications
#' @return 
#' none
#' @note This auxiliary function helps create a new package. The orignal Rd files will be replaced with new ones.
#' @export
#' @seealso \code{\link{xRdWrap}}
#' @include xRdWrap.r
#' @examples
#' # xRdWrap(path="./XGR/man", remove.dontrun=FALSE)

xRdWrap <- function(path="./XGR/man", remove.dontrun=FALSE)
{
	for(nm in list.files(path, pattern="\\.Rd$")) {
		input_file <- file.path(path, nm)
		message(sprintf(input_file), appendLF=T)
		x <- readLines(input_file)
		# make sure those lines starting with '#' are not wrapped
		y <- vector()
		flag_examples <- F
		flag_examples_dontrun <- F
		num_examples_brace <- 0
		for(i in 1:length(x)){
			if(length(grep("^#", x[i]))==0 & length(grep("^\\\\title",x[i]))==0){
				# line is wrapped
				y <- c(y, strwrap(x[i]))
			}else{
				# line starting with '#' is not wrapped
				y <- c(y, x[i])
			}
		
			######################
			if(remove.dontrun==T){
				# mark first-encountered examples
				if(x[i]=="\\examples{"){
					flag_examples <- T
				}
				if(flag_examples==T){
					# remove the dontrun mark (if any)
					if(x[i]=="\\dontrun{"){
						y <- y[-1*length(y)]
						flag_examples_dontrun <- T
					}
					# execute the following when detecting the dontrun mark
					if(flag_examples_dontrun==T){
						# count the number of brace within examples
						if(x[i]=="}"){
							num_examples_brace <- num_examples_brace+1
						}
						# remove the last brace within examples
						if(num_examples_brace==1){
							y <- y[-1*length(y)]
							flag_examples <- F
						}
					}
				}
			}
			######################        
		}
	
		fileConn <- file(input_file)
		writeLines(y, fileConn)
		close(fileConn)
	}
}
