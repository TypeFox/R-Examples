#' List Jobs
#'
#' A list of jobs owned by the calling account. 
#' 
#' @export
#' @references \url{https://shreddr.captricity.com/developer/api-reference/#v1-jobs}
#' @examples \dontrun{
#' list_jobs()
#' }

list_jobs <- function() {
    
    res <- captr_GET("job/", NULL)

    cat("No. of jobs:", nrow(res), "\n")

    return(invisible(res))

}

 
