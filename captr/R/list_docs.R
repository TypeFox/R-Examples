#' List Documents
#'
#' A list of document resources owned by the calling account.
#' 
#' @export
#' @references \url{https://shreddr.captricity.com/developer/api-reference/#v1-documents}
#' @examples \dontrun{
#' list_docs()
#' }

list_docs <- function() {
    
    res <- captr_GET("document/", NULL)

    cat("No. of document resources:", nrow(res), "\n")

    return(invisible(res))
}


