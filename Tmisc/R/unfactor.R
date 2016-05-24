#' Unfactor a data.frame
#' 
#' Did you forget to pass \code{stringsAsFactors=FALSE}? This converts factor
#' variables to characters in a dataframe.
#' 
#' @author \url{https://github.com/Dasonk}
#' @keywords NA
#'   
#' @param df The dataframe you wish to change the factors into characters.
#'   
#' @return A data.frame with factors converted to characters.
#' 
#' @examples
#' df <- data.frame(a = letters[1:5], x = 1:5, y = LETTERS[1:5], stringsAsFactors = TRUE)
#' str(df)
#' df <- unfactor(df)
#' str(df)
#'   
#' @export

unfactor <- function(df) {
    # Find the factors
    id <- sapply(df, is.factor)
    # Convert to characters
    df[id] <- lapply(df[id], as.character)
    df
}