#' Source all file with extension .r or .R
#' @name sourceAll
#' @aliases sourceAll
#' @title Source all the R files of a directory
#' @param path path of the directory
#' @param ... other arguments to source
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \link{repCol}
#' @export
#' @examples
#' \dontrun{   
#' sourceAll()
#' }
sourceAll<-function(path=".",...){
	file_list <- list.files(path=path)
	file_list <-file_list[grep("*.r$|*.R$",file_list)]
	file_list<-path%+%"/"%+%file_list
	cat("Loading...\n")
	junk<-lapply(file_list,function(x) {cat(" ",x,"\n");source(x,...) })	
	cat("Done\n")
}


