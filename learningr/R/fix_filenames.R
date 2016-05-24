#' Convert filenames to match those in the book
#'
#' Some filenames have been altered in order to comply with portability
#' requirements on CRAN.  This function converts the filenames between 
#' the CRAN forms and the book forms.
#'
#' @param x Either ``tocran'' or ``tobook''.
#' @param dir Directory containing the files.
#' @return A logical vector of length 4, \code{TRUE} for each file whose
#' name was changed.
#'
#' @examples
#' \dontrun{
#' #To convert the files to the book form, use:
#' fix_filenames("tobook")
#' #The files were converted to CRAN form using:
#' fix_filenames("tocran", "learningr/inst/extdata")
#' }
#' @export
fix_filenames <- function(x = c("tobook", "tocran"), dir = system.file("extdata", packages = "learningr"))
{
  x <- match.arg(x)
  book_filenames <- c(
    "Alpe d'Huez.xls",
    "Jamaican Cities.json",
    "Shakespeare's The Tempest, from Project Gutenberg pg2235.txt",
    "multi-drug-resistant gonorrhoea infection.xls"
  )
  cran_filenames <- make.names(book_filenames)
  result <- if(x == "tobook")
  {
    file.rename(
      file.path(dir, cran_filenames), 
      file.path(dir, book_filenames)
    )
  } else
  {
    file.rename(
      file.path(dir, book_filenames), 
      file.path(dir, cran_filenames)
    )
  }
  names(result) <- book_filenames
  result
}
