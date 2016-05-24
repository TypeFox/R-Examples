
#' @importFrom XML xmlTreeParse xmlRoot xmlToList
preprocess_response <- function(x) {
    y <- XML::xmlRoot(XML::xmlTreeParse(x))

    y <- list(
        request = XML::xmlToList(y[['request']]),
        message = XML::xmlToList(y[['message']]),
        response = y[['response']]
    )

    return(y)
}
