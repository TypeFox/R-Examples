
#' @title userInput let the user type in a character string
#' 
#' @description this is useful for asking the user for a username or password, so that it goes directly to a variable and doesn't get inadvertently saved into an R script. There are no arguments. This will only return a character string. This is low key, don't bother using it for data entry. Just type characters, no need to put it in quotes, pressing enter will cause the function to return. Output will not be printed to the console, but it can be assigned directly. This is useful to have as an auxiliary function in case multiple calls to functions such as \code{readHMDweb()} are desired.
#' 
#' @param silent logical should a little prompt be given, telling the user to enter text in the console?
#' 
#' @export
#' 
#' @return a character string, as given by the user.
#' 
#' @examples
#' ### mypassword <- userInput()
#' ### myusername <- userInput()
#' ### DAT <- readHMDweb("USA","mltper_1x1",mypassword,myusername)
#' 
userInput <- function(silent = FALSE){
    if (interactive()){
      if(!silent){
        cat("\ntype in a single character string.")
      }
      out <- scan(file = "", n = 1, what = "character")
    } else {
        stop("User input only works for an interactive R session")
    }
    invisible(out)
}

   