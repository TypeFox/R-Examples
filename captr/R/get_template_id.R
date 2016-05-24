#' Set Template ID
#' 
#' To digitize documents, create a template using the Captricity Web UI at \url{https://shreddr.captricity.com/job/} 
#' The template tells Captricity which data to get from where in the document.
#' Get the template id of the template
#' 
#' @export
#' @references \url{https://shreddr.captricity.com/developer/}
#' @examples \dontrun{
#' get_template_id()
#' }

get_template_id <- function() {
    
    captr_CHECKAUTH()
 
    h <- new_handle()
    handle_setopt(h,  customrequest = "POST")
    handle_setheaders(h, "Captricity-API-Token" = Sys.getenv('CaptricityToken'))

    tag_con    <- curl_fetch_memory("https://shreddr.captricity.com/api/v1/document/", handle=h)
    tag        <- fromJSON(rawToChar(tag_con$content))
    tag
    return(invisible(tag))

}

