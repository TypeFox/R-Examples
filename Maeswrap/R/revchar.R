#' Reverse a character string.
#' 
#' Reverses the characters in a character string, unless a vector is supplied,
#' in which case reverses the element of the vector.
#' 
#' When a character vector of length > 1 is provided, reverses the elements of
#' the vector, not the characters itself.
#' 
#' @param x A character vector (typically of length 1).
#' @return A vector.
#' @author Remko Duursma
#' @seealso
#' \code{\link{rev.default}},\code{\link{substr}},\code{\link{strsplit}}
#' @references None.
#' @keywords methods
#' @examples
#' 
#' 
#' \dontrun{
#' # Take a substring, counting from the end:
#' substrfromend <- function(x,start,stop)revchar(substr(revchar(x),start,stop))
#' substrfromend('filename.txt', 1,3)
#' 
#' # Check if a word is a palindrome:
#' s <- 'saippuakivikauppias'
#' s == revchar(s)
#' 
#' # A semordnilap:
#' cat('I am so stressed, I need to eat', revchar('stressed'),'\n')
#' }
#' 
revchar <- function(x){
if(length(x) > 1){
    rev.default(x)
 } else {
paste(rev.default(strsplit(x,"")[[1]]),collapse="")
}
}
