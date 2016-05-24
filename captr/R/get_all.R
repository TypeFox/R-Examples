#' Get All the Data For a Particular Job in a csv
#' 
#' Get all the data
#' 
#' @param job_id ID for the job
#' @param output_dir output directory
#' @export
#' @references \url{https://shreddr.captricity.com/developer/}
#' @examples \dontrun{
#' get_all(job_id ="job_id")
#' }

get_all <- function(job_id ="", output_dir="./") {
   
    captr_CHECKAUTH()
 
    if ( is.null(job_id) | job_id=="") stop("Provide a Valid Job ID.")

    h <- new_handle()
    handle_setopt(h,  customrequest = "GET")
    handle_setheaders(h, "Captricity-API-Token" = Sys.getenv('CaptricityToken'))

    tag_con    <- curl_fetch_memory(paste0("https://shreddr.captricity.com/api/v1/job/", job_id, "/csv/"), handle=h)
    curl_download(tag_con$url, destfile=paste0(output_dir, job_id, ".csv"))    
}

