#' Quickly Creating BigML Datasets
#' @export
#' @family dataset methods
#' @references \url{https://bigml.com/developers/datasets}
#' @family quick methods
#' @param data A matrix or data frame containing data to upload to bigml.
#' @param name A string giving the name for the dataset.
#' @param fields A vector of names in \code{data} that should be used for
#'	creating the dataset.
#' @param size A numeric value giving the amount (in bytes) of the source
#'	to use.
#' @template dots
#' @template dataset_return
#' @template author
#' @details quickDataset will take its "data" dataframe argument and attempt
#' 	to create an equivalent BigML dataset using \code{\link{quickSource}}.
#'	R "numeric" class fields will become "numeric" fields in the BigML
#'	dataset.  R "character" class fields become "text" fields.  R "factor"
#'	fields become "categorical" fields. However, if there are too many
#'	factors, BigML may convert the field to text.  It is possible to specify
#'	the fields to include using the \code{fields} argument.  This can be a
#'	a simple list of names that were present in the \code{data} argument.
#'	See references for more details.
#' @examples
#' \dontrun{
#' # simple example
#' iris.d = quickDataset(iris)
#' # configure a number of different parameters
#' iris.d2 = quickDataset(iris, fields = c('Species', 'Sepal.length'),
#'	name='test', size=10000)
#' }
quickDataset <-
function (data, fields = names(data), name = paste(deparse(substitute(data)),
    "'s dataset", sep = ""), size = NULL, ...)
{
    option = list()
    sresponse = quickSource(data,  name = deparse(substitute(data)),
		flatten = F, ...)
    if (is.null(name))
        name = paste(sresponse$file_name, "'s dataset", sep = "")
    classes = lapply(data, class)
    type_classes = lapply(classes, function(x) {
        if (x == "numeric")
            return("numeric")
        else if (x == "factor")
            return("categorical")
        else return("text")
    })
    option$name = name
    option$source = sresponse$resource
    if (is.null(size))
        option$size = size
    message("Dataset creation in progress...")
    return(.basic_api(.DATASET_URL)$postJson(option, ...))
}
