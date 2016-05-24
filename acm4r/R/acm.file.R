#' acm.file.copy (the modified file.copy function in acm)
#'
#' acm.file.copy behaves similar to file.copy except that 
#' the program terminates and shows the already existing file 
#' which cannot be overwritten by "to" if "overwrite = FALSE"
#'
#' @param from      see file.copy
#' @param to        see file.copy
#' @param overwrite see file.copy
#' @param recursive see file.copy
#' @param copy.mode see file.copy
#' @return none
#'
#' @author XiaoFei Zhao \email{xiaofei.zhao@@mail.mcgill.ca}
acm.file.copy <- function(
		from, to, 
		overwrite = FALSE, 
		recursive = FALSE,
		copy.mode = TRUE)
{
	if (overwrite | !file.exists(to)) {
		file.copy(from = from, to = to , overwrite = overwrite, 
				recursive = recursive, copy.mode = copy.mode)
	}
	else {
		print(paste("In directory ", getwd(), " ",               sep = "'"))
		print(paste("File ",         to,      " already exists", sep = "'"))
		print("Set delete = TRUE to allow overwriting of already existed file")
		stop()
	}
}



#' acm.file.exists (the modified file.exists function in acm)
#'
#' acm.file.exists behaves similar to file.exists except that 
#' the program terminates and shows the already existing file 
#' if the file already exists
#'
#' @param filename the name of the file which is tested 
#' for its existence
#' @return none
#'
#' @author XiaoFei Zhao \email{xiaofei.zhao@@mail.mcgill.ca}
acm.file.exists <- function(filename)
{
  if (FALSE == file.exists(filename))
  {
    print(paste("In directory ", getwd(),  " ",               sep = "'"))
    print(paste("File ",         filename, " does not exist", sep = "'"))
    stop()
  }
}





#' acm.file.checkdelete
#'
#' acm.file.checkdelete abort the program if the file with the given filename already exists
#'
#' @param filename the name of the file which is tested 
#' for its existence
#' @return none
#'
#' @author XiaoFei Zhao \email{xiaofei.zhao@@mail.mcgill.ca}
acm.file.checkdelete <- function (filename) 
{
	if (file.exists(filename)) 
	{
		print(paste("In directory ", getwd(),  " ",               sep = "'"))
		print(paste("File ",         filename, " already exists", sep = "'"))
		print("Set delete = TRUE to allow overwriting of already existed file")
		stop()
	}
}