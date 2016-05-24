#'About the Berkeley Ecoinformatics Engine
#'
#' Function returns the current status of fast-evolving API. Returns endpoints and category. Default return is a \code{list} but one can also request a nicely formatted \code{data.frame} by setting the \code{as.df} argument to \code{TRUE}.
#' @param as.df  \code{FALSE} Returns a list unless this set to \code{TRUE}
#' @param type  The type of end point. Options include \code{data}, \code{meta-data}, and \code{actions}
#' @return \code{list}
#' @export
#' @importFrom httr GET content warn_for_status
#' @importFrom plyr ldply
#' @importFrom jsonlite fromJSON
#' @examples  
#' ee_about()
#' # set as.df = FALSE to return a list rather than a data.frame
#' ee_about(as.df = FALSE)
#' # You can also filter by methods by data, meta-data, and actions.
#' # ee_about(type = "data")
#' # ee_about(type = "meta-data")
#' # ee_about(type = "actions")
ee_about <- function(as.df = TRUE, type = NA) {
about_url <- paste0(ee_base_url(), "?format=json")
about_call <- GET(about_url)
warn_for_status(about_call)
about <- content(about_call)
if(!as.df) {
    return(about)
    } else {
        about_df <- lapply(about, function(f) {
             res <- data.frame(cbind(f))
            } )
        about_df <- ldply(about_df)
        names(about_df) <- c("type", "endpoint")
        about_df$endpoint <- unlist(about_df$endpoint)
        if(!is.na(type)) {
                about_df <- switch(type, 
                        data = subset(about_df, type == "data"),
                        actions = subset(about_df, type == "actions"),
                        "meta-data" = subset(about_df, type == "meta-data"),
                        )
        }
        return(about_df)

}
}