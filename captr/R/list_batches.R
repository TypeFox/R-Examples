#' List Batches
#'
#' A list of batches owned by the calling account.
#' 
#' @export
#' @references \url{https://shreddr.captricity.com/developer/api-reference/#v1-batches}
#' @examples \dontrun{
#' list_batches()
#' }

list_batches <- function() {
      
    res <- captr_GET("batch/", NULL)

    cat("No. of batches:", nrow(res), "\n")

    return(invisible(res))

}

