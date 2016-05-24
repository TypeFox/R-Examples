#' Provide Feedback on Tags of an Image
#'
#' Suggest tags, suggest removal of tags, suggest (dis)similar documents, or let clarifai know which tag was clicked on 
#' 
#' @param feedback_type must be one of the following types:
#' add_tags vector of suggested tags; default
#' remove_tags vector of tags you suggest be removed
#' search_click which tag(s) where clicked on
#' similar_docids similar documents 
#' dissimilar_docids dissimilar documents
#' @param file_path image on which the feedback is being provided
#' @param feedback_value what tags to add, remove etc.
#' 
#' @export
#' @references \url{https://developer.clarifai.com/}
#' @examples \dontrun{
#' feedback(file_path="path_to_image", feedback_type="add_tags", feedback_value="suggested_tag")
#' }

feedback <- function(file_path = "", feedback_type='add_tags', feedback_value="") {
   
    clarifai_check_token()

    if(!file.exists(file_path)) stop("File Doesn't Exist. Please check the path.")
    if(identical(feedback_value, "")) stop("Provide something other than an empty string as feedback.")

    h <- new_handle()
	handle_setopt(h,  customrequest = "POST")
	handle_setheaders(h, "Authorization" = paste0("Bearer ", Sys.getenv("ClarifaiToken")))

	switch(feedback_type, 
						add_tags     = {handle_setform(h, encoded_image = form_file(file_path), add_tags          = feedback_value)},
						remove_tags  = {handle_setform(h, encoded_image = form_file(file_path), remove_tags       = feedback_value)},
						search_click = {handle_setform(h, encoded_image = form_file(file_path), search_click      = feedback_value)},
						similar_docids = {handle_setform(h, encoded_image = form_file(file_path), similar_docids    = feedback_value)},
						dissimilar_docids = {handle_setform(h, encoded_image = form_file(file_path), dissimilar_docids = feedback_value)})
	

	tag_con    <- curl_fetch_memory("https://api.clarifai.com/v1/feedback/", handle=h)
	tag        <- fromJSON(rawToChar(tag_con$content))
	tag

	return(invisible(tag))

}

