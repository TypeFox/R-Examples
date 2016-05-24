#' Delete a Batch
#' 
#' @param batch_id ID for the batch
#' @export
#' @references \url{https://shreddr.captricity.com/developer/api-reference/#v1-batches}
#' @examples \dontrun{
#' delete_batch("batch_id")
#' }

delete_batch <- function(batch_id="") {
   
    captr_CHECKAUTH()

    if ( is.null(batch_id) | identical(batch_id, "")) stop("Provide a Valid Batch ID.")

    h <- new_handle()
    handle_setopt(h,  customrequest = "DELETE")
    handle_setheaders(h, "Captricity-API-Token" = Sys.getenv('CaptricityToken'))

    tag_con    <- curl_fetch_memory(paste0("https://shreddr.captricity.com/api/v1/batch/", batch_id), handle=h)
    tag        <- fromJSON(rawToChar(tag_con$content))
    
    if (tag$status=="success") {
    	cat("Batch", batch_id, "successfully deleted \n")
    } else {
    	cat("Error:", tag$status, "\n")
    }

    return(invisible(tag))
    
}

