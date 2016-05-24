#' clear console by function - does not work!
#' 
#' clear console (CTRL + L) using a function call. Does not work, as \code{rcom} is not available!
#' 
#' 
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2014
#' @references
#' \url{http://r.789695.n4.nabble.com/how-to-clear-screen-in-R-console-td793936.html}
#' @keywords programming utilities
#' @export
#' @examples
#' 
#' #cls()
#' 
cls <- function()

{
    "Does currently not Work"
    #if (!require(rcom)) stop("Package 'rcom' is required for 'cls()'")
    #wsh <- comCreateObject("Wscript.Shell")
    #comInvoke(wsh, "SendKeys", "\014")
    #invisible(wsh)
}

# cls() # test 
