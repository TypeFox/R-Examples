#' Quickly Creating BigML Sources
#' @export
#' @family source methods
#' @references \url{https://bigml.com/developers/sources}
#' @family quick methods
#' @param data A matrix or data frame containing data to upload to bigml.
#' @param name A string giving the name of the source.
#' @param header A logical value indicating whether to use the first row
#' 	of data as a header row.
#' @param locale A string indicating the desired locale.
#' @param missing_tokens A vector listing strings that should be treated as
#'	missing.
#' @param quote A string giving the quote character to use.
#' @param trim A logical value indicating whether to trim white space.
#' @param flatten A logical value indicating whether to flatten the response
#'	into a data frame.
#' @template dots
#' @template source_return
#' @details quickSource will take its "data" dataframe argument and attempt
#' 	to create an equivalent BigML source.  It does this by converting the
#'	dataframe to a csv file, compressing it, and uploading it directly to
#'	BigML.  Generally, it's better to use \code{\link{quickDataset}}, since
#'	this method attempts to preserve any type information in the data frame.
#' @note It is not currently possible to retrieve the original file from
#'	BigML, but it is possible to delete it.
#' @template author
quickSource <-
function (data, name = deparse(substitute(data)),
	header = !is.null(names(data)), locale = "en-US",
	missing_tokens = c("NA"), quote = "\"", trim = TRUE, flatten = TRUE, ...){

	file_name = paste(tempdir(), "/", name, ".csv.gz", sep = "")
	file_handle = gzfile(file_name, "w")
    write.table(data, file = file_handle, row.names = F,
        col.names = header, sep = ",")

	close(file_handle)
	return (createSource(file_name, name=name, header=header, locale=locale,
		missing_tokens=missing_tokens, quote=quote,
		trim=trim, flatten=flatten, ...))
}