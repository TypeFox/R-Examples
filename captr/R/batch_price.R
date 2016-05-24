#' Price of a Batch
#' 
#' The price for processing the batch.
#' 
#' @param batch_id ID for the batch
#' @export
#' @references \url{https://shreddr.captricity.com/developer/}
#' @examples \dontrun{
#' batch_price("batch_id")
#' }

batch_price <- function(batch_id="") {

    captr_CHECKAUTH()
   
    if (is.null(batch_id) | identical(batch_id, "")) stop("Provide a Valid Batch ID.")

    h <- new_handle()
    handle_setopt(h,  customrequest = "GET")
    handle_setheaders(h, "Captricity-API-Token" = Sys.getenv('CaptricityToken'))

    tag_con    <- curl_fetch_memory(paste0("https://shreddr.captricity.com/api/v1/batch/", batch_id, "/price"), handle=h)
    tag        <- fromJSON(rawToChar(tag_con$content))
    tag
    return(invisible(tag))

}

