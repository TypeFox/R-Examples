#' Delete a Job
#' 
#' @param job_id ID for the job (which you get from related_job_id field of submit_batch)
#' @export
#' @references \url{https://shreddr.captricity.com/developer/api-reference/#v1-batches}
#' @examples \dontrun{
#' delete_job("job_id")
#' }

delete_job <- function(job_id="") {
   
    captr_CHECKAUTH()

    if ( is.null(job_id) | identical(job_id, "")) stop("Provide a Valid Job ID.")

    h <- new_handle()
    handle_setopt(h,  customrequest = "DELETE")
    handle_setheaders(h, "Captricity-API-Token" = Sys.getenv('CaptricityToken'))

    tag_con    <- curl_fetch_memory(paste0("https://shreddr.captricity.com/api/v1/job/", job_id), handle=h)
    tag        <- fromJSON(rawToChar(tag_con$content))
    
    cat(tag[[1]], "\n")
    
    return(invisible(tag))
    
}

