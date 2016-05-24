#' Open the current working directory on mac
#' 
#' Opens the current working directory on mac.
#' 
#' @author Stephen Turner
#' @keywords NA
#' 
#' @examples
#' \dontrun{
#' o()
#' }
#' 
#' @export
o <- function() {
    if(Sys.info()[1]=="Darwin") {
        message(getwd())
        system("open .")
    } else {
        warning("This only works on Mac.\nKnow how to use Windows? Submit a PR at:\nhttps://github.com/stephenturner/Tmisc")
    }
}