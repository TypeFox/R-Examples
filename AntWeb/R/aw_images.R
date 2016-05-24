#' aw_images
#'
#' Download ant images based on time elapsed and/or type.
#' @param since number of days in the past to query
#' @param  img_type h for head, d for dorsal, p for profile, and l for label. If a img_type is not specified, all images are retrieved.
#' @export
#' @return data.frame
#' @examples \dontrun{
#' z <- aw_images(since = 5)
#' z1 <- aw_images(since = 5, img_type = "d")
#'}
aw_images <- function(since = NULL, img_type = NULL) {
	img <- "true"
	base_url <- "http://www.antweb.org/api/v2"
	args <- z_compact(as.list(c(since = since, img = img, img_type = img_type)))
	results <- GET(base_url, query = args)
	warn_for_status(results)
	data <- fromJSON(content(results, "text"))

	data_list <- lapply(data, function(z) {
			collection_date <- z[[1]][[1]]
			upload_date <- z[[2]][1]
			photo_list <- z[[2]][2]
			photo_data <- lapply(photo_list, function(il) {
				imgg <- list()
				for(i in 1:length(il)) {
					imgg[[i]] <-  data.frame(t(unlist(il[[i]]$img)))
					imgg[[i]]$img_type <- unlist(names(il[i]))
				}

				img_data <- do.call(rbind, imgg)
				img_data$collection_date <- collection_date
				img_data$upload_date <- upload_date
				names(img_data) <- c("high", "med", "low", "thumbnail", "img_type", "catalog_id", "upload_date")
				img_data
			})
			photo_data_df <- do.call(rbind, photo_data)

	})
  final_df <- do.call(rbind, data_list)
rownames(final_df) <- NULL
final_df

}


