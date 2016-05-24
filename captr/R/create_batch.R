#' Create Batch
#'
#' Create a new batch. 
#' 
#' @param batch_name name of the batch; Required; character
#' 
#' @return List of length 26. Includes information like created_by, user_id, etc.
#' 
#' @export
#' @references \url{https://shreddr.captricity.com/developer/api-reference/#v1-batch}
#' @examples \dontrun{
#' create_batch(batch_name="name_of_batch")
#' }

create_batch <- function(batch_name="") {
    
    captr_CHECKAUTH()
  
    h <- new_handle()
    handle_setopt(h,  customrequest = "POST")
    handle_setheaders(h, "Captricity-API-Token" = Sys.getenv('CaptricityToken'))
    handle_setform(h, name= batch_name)

    tag_con    <- curl_fetch_memory("https://shreddr.captricity.com/api/v1/batch/", handle=h)
    tag        <- fromJSON(rawToChar(tag_con$content))
    tag
    return(invisible(tag))
}


