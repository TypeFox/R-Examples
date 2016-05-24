#' Upload Image
#' 
#' Upload an image to captricity.
#' 
#' @param batch_id ID for the batch
#' @param path_to_image Path to the image you want OCRd
#' @export
#' @references \url{https://shreddr.captricity.com/developer/}
#' @examples \dontrun{
#' upload_image("batch_id", "path_to_image")
#' }

upload_image <- function(batch_id="", path_to_image="") {
    
    captr_CHECKAUTH()
 
    if ( is.null(batch_id) | identical(batch_id, "")) stop("Provide a Valid Batch ID.")
        
    if (!file.exists(path_to_image)) stop("File Doesn't Exist. Please check the path.")

    h <- new_handle()
    handle_setopt(h,  customrequest = "POST")
    handle_setheaders(h, "Captricity-API-Token" = Sys.getenv('CaptricityToken'))
    handle_setform(h, uploaded_file = form_file(path_to_image))

    tag_con    <- curl_fetch_memory(paste0("https://shreddr.captricity.com/api/v1/batch/", batch_id, "/batch-file/"), handle=h)
    tag        <- fromJSON(rawToChar(tag_con$content))
    tag
    return(invisible(tag))

}

