#' Tag Images on the Computer
#'
#' @param file_paths a vactor of paths to image file(s) for which you want tags
#' @inheritParams tag_image_urls
#' 
#' @return If meta is TRUE, a named list of length 4 with following elements is returned: 
#' \code{status_code} OK or not
#' \code{status_msg}  Successful or not
#' \code{meta} Named list of 1 containing another list named \code{tag}
#' Sublist \code{tag} has three elements: timestamp, model and config
#' \code{results} is a data.frame of length 6 and 1 row. Column names are:
#' \code{docid}, \code{status_code}, \code{status_msg}, \code{local_id} and 
#' a data.frame named tag which has a data.frame result which contains two columns: 
#' labels and probabilities
#' 
#' If meta is FALSE and simplify is TRUE,
#' a data.frame with three columns: file_paths, labels and probs returned
#'
#' If meta is FALSE and simplify is FALSE,
#' a data.frame with two columns carrying a vector of labels, vector of probs is returned
#' for each image
#' 
#' @export
#' @references \url{https://developer.clarifai.com/}
#' @seealso \code{\link{tag_image_urls}}
#' @examples \dontrun{
#' tag_images(file_paths="path_to_image")
#' }

tag_images <- function(file_paths=NULL, meta=FALSE, simplify=TRUE) {
	
    clarifai_check_token()
        
    if (! all(file.exists(file_paths))) stop("File Doesn't Exist. Please check the path.")

    paths <- lapply(file_paths, form_file)
    names(paths) <- rep("encoded_image", length(paths))

    h <- new_handle()
	handle_setopt(h,  customrequest = "POST")
	handle_setheaders(h, "Authorization" = paste0("Bearer ", Sys.getenv("ClarifaiToken")))
	handle_setform(h, .list=paths)

	tag_con    <- curl_fetch_memory("https://api.clarifai.com/v1/tag/", handle=h)
	tag        <- fromJSON(rawToChar(tag_con$content))
	
	
	if (!meta) {
		
		if (simplify) {
        
          # Assumes 20 out
		   tags <- lapply(tag$results$result$tag[,1], unlist)
		   probs <- lapply(tag$results$result$tag[,2], unlist)
		   tags_probs <- do.call(rbind, Map(cbind, tags, probs))
		   names(tags_probs) <- c("tags", "probs")
		   len <- sapply(probs, length)
		   tags_probs_imgs <- data.frame(file_paths=rep(file_paths, len), tags_probs)
		   return(invisible(tags_probs_imgs))
		}

		return(invisible(tag$results$result$tag))

	}

	return(invisible(tag))

}

