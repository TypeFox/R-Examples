#' User Profile
#'
#' Provides account information for the provided API Token.
#' 
#' @export
#' @references \url{https://shreddr.captricity.com/developer/api-reference/#v1-user_profile}
#' @examples \dontrun{
#' list_batches()
#' }

list_batches <- function()
{
    
    captr_CHECKAUTH()
  
    h <- new_handle()
    handle_setopt(h,  customrequest = "GET")
    handle_setheaders(h, "Captricity-API-Token" = Sys.getenv('CaptricityToken'))

    tag_con    <- curl_fetch_memory("https://shreddr.captricity.com/api/v1/account/user-profile", handle=h)
    tag        <- fromJSON(rawToChar(tag_con$content))
    tag

    return(invisible(tag))

}


