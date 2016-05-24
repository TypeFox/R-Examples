#' Tag Remote Images
#'
#' @param img_urls A string or character vector of URLs of images for which you want tags
#' @param meta Boolean that toggles whether or not you want the entire object. 
#' @param simplify Boolean that toggles whether or not you want a simplified data frame with
#' each label and associated probability in a separate row. Default is TRUE.
#' 
#' The entire object returned by the API contains a lot of meta data. 
#' By default a simplified data.frame with img_url, associated labels, and probabilities is returned. 
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
#' a data.frame with three columns: img_urls, labels and probs returned
#'
#' If meta is FALSE and simplify is FALSE,
#' a data.frame with two columns carrying a vector of labels, vector of probs is returned
#' for each image
#' 
#' @export
#' @references \url{https://developer.clarifai.com/}
#' @seealso \code{\link{tag_images}}
#' @examples \dontrun{
#' tag_image_urls(img_url="url_of_image")
#' }

tag_image_urls <- function(img_urls=NULL, meta=FALSE, simplify=TRUE) {
    
    clarifai_check_token()
    
    urls <- as.list(img_urls)
    names(urls) <- rep("url", length(urls))

    h <- new_handle()
	handle_setopt(h,  customrequest = "POST")
	handle_setheaders(h, "Authorization" = paste0("Bearer ", Sys.getenv("ClarifaiToken")))
	handle_setform(h, .list = urls)
	tag_con    <- curl_fetch_memory("https://api.clarifai.com/v1/tag/", handle=h)
	tag        <- fromJSON(rawToChar(tag_con$content))

	status_msg <- tag$results$status_msg=="OK"

	if (!meta) {
		
		if (simplify) {
        
          # Assumes 20 out
		   tags <- lapply(tag$results$result$tag[,1], unlist)
		   probs <- lapply(tag$results$result$tag[,2], unlist)
		   tags_probs <- do.call(rbind, Map(cbind, tags, probs))
		   len <- sapply(probs, length)
		   tags_probs_imgs <- data.frame(img_urls=rep(img_urls, len), tags_probs)
		   names(tags_probs_imgs) <- c("img_url", "tags", "probs")
		   return(invisible(tags_probs_imgs))
		}

		return(invisible(tag$results$result$tag))

	}

	return(invisible(tag))
}

