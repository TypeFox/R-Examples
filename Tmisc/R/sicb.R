#' Write sessionInfo to the clipboard
#' 
#' Writes output of \code{sessionInfo()} to the clipboard. Only works on Mac.
#' 
#' @author Stephen Turner
#' @keywords sessioninfo
#' @import utils
#' 
#' @examples
#' \dontrun{
#' # Write sessionInfo() to the clipboard on mac.
#' sicb()
#' }
#' 
#' @export
sicb <- function() {
    # Check to make sure you're running on mac.
    if (Sys.info()[1]=="Darwin") {
        capture.output(sessionInfo(), file=pipe("pbcopy"))
    } else {
        warning("This only works on Mac.\nKnow how to use Windows? Submit a PR at:\nhttps://github.com/stephenturner/Tmisc")
    }
}
