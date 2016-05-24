#' @title Lists all the dimensions and metrics for a particular report type
#'
#' @description
#' This dataset represents all of the dimensions and metrics for the reporting API with their attributes. Attributes returned include UI name, description, segments support, etc.
#'
#' @param reportType character. Report type. Allowed Values: "ga". Where "ga" corresponds to the Core Reporting API.
#'
#' @return A data.frame contains dimensions and metrics for a particular report type.
#' \item{id}{Parameter name.}
#' \item{type}{The type of column: \code{DIMENSION}, \code{METRIC}.}
#' \item{dataType}{The type of data this column represents: \code{STRING}, \code{INTEGER}, \code{PERCENT}, \code{TIME}, \code{CURRENCY}, \code{FLOAT}.}
#' \item{group}{The dimensions/metrics group the column belongs to.}
#' \item{status}{The status of the column: \code{PUBLIC}, \code{DEPRECATED}.}
#' \item{uiName}{The name/label of the column used in user interfaces (UI).}
#' \item{description}{The full description of the column.}
#' \item{allowedInSegments}{Indicates whether the column can be used in the segment query parameter.}
#' \item{addedInApiVersion}{API version with this param was added.}
#' \item{replacedBy}{The replacement column to use for a column with a \code{DEPRECATED} status.}
#' \item{calculation}{Only available for calculated metrics. This shows how the metric is calculated.}
#' \item{minTemplateIndex}{Only available for templatized columns. This is the minimum index for the column.}
#' \item{maxTemplateIndex}{Only available for templatized columns. This is the maximum index for the column.}
#' \item{premiumMinTemplateIndex}{Only available for templatized columns. This is the minimum index for the column for premium properties.}
#' \item{premiumMaxTemplateIndex}{Only available for templatized columns. This is the maximum index for the column for premium properties.}
#'
#' @seealso \code{\link{shiny_dimsmets}} \code{\link{get_ga}}
#'
#' @references
#' \href{https://developers.google.com/analytics/devguides/reporting/metadata/v3/}{Google Analytics Metadata API}
#'
#' \href{https://developers.google.com/analytics/devguides/reporting/core/dimsmets}{Core Reporting API - Dimensions & Metrics Reference}
#'
#' @examples
#' \dontrun{
#' ga_meta <- list_dimsmets("ga")
#' # a count of parameters types
#' table(ga_meta$type)
#' # parameters groups
#' table(ga_meta$group)
#' # get a deprecated parameters was replaced by
#' subset(ga_meta, status == "DEPRECATED", c(id, replacedBy))
#' # get a calculation metrics
#' subset(ga_meta, !is.na(calculation), c(id, calculation))
#' # get a not deprecated metrics from user group
#' subset(ga_meta, group == "User" & type == "METRIC" & status != "DEPRECATED", id)
#' # get parameters allowed in segments
#' subset(ga_meta, allowedInSegments, id)
#' }
#'
#' @importFrom httr GET
#' @include url.R
#' @include request.R
#' @export
list_dimsmets <- function(reportType = "ga") {
    url <- get_url(c("metadata", reportType, "columns"))
    response <- GET(url)
    json_content <- process_response(response)
    res <- json_content$items
    res$kind <- NULL
    names(res) <- gsub("attributes.", "", names(res), fixed = TRUE)
    res$allowedInSegments <- as.logical(res$allowedInSegments)
    tocenvert <- grep("TemplateIndex", names(res), fixed = TRUE)
    res[tocenvert] <- lapply(res[tocenvert], as.integer)
    class(res) <- c("tbl_df", "tbl", "data.frame")
    return(res)
}

#' @title The Shiny Dimensions & Metrics Explorer
#'
#' @description
#' The dimensions and metrics explorer lists and describes all the dimensions and metrics available through the Core Reporting API. This app deployed to the \url{https://artemklevtsov.shinyapps.io/ga-dimsmets}.
#'
#' @seealso \code{\link{list_dimsmets}} \code{\link{get_ga}}
#'
#' @export
shiny_dimsmets <- function() {
    appDir <- system.file("shiny-examples", "01-dimsmets", package = "RGA")
    if (appDir == "")
        stop("Could not find example directory. Try re-installing 'RGA' package.", call. = FALSE)
    shiny::runApp(appDir, display.mode = "normal")
}

#' @title Lists all columns for a Google Analytics core report type
#'
#' @usage
#' ga
#'
#' @description
#' This dataset represents all of the dimensions and metrics for the reporting API with their attributes. Attributes returned include UI name, description, segments support, etc.
#'
#' @format
#' A data frame with 436 rows and 14 variables containing the following columns:
#' \describe{
#' \item{id}{Parameter name.}
#' \item{type}{The type of column: \code{DIMENSION}, \code{METRIC}.}
#' \item{dataType}{The type of data this column represents: \code{STRING}, \code{INTEGER}, \code{PERCENT}, \code{TIME}, \code{CURRENCY}, \code{FLOAT}.}
#' \item{group}{The dimensions/metrics group the column belongs to.}
#' \item{status}{The status of the column: \code{PUBLIC}, \code{DEPRECATED}.}
#' \item{uiName}{The name/label of the column used in user interfaces (UI).}
#' \item{description}{The full description of the column.}
#' \item{allowedInSegments}{Indicates whether the column can be used in the segment query parameter.}
#' \item{addedInApiVersion}{API version with this param was added.}
#' \item{replacedBy}{The replacement column to use for a column with a \code{DEPRECATED} status.}
#' \item{calculation}{Only available for calculated metrics. This shows how the metric is calculated.}
#' \item{minTemplateIndex}{Only available for templatized columns. This is the minimum index for the column.}
#' \item{maxTemplateIndex}{Only available for templatized columns. This is the maximum index for the column.}
#' \item{premiumMinTemplateIndex}{Only available for templatized columns. This is the minimum index for the column for premium properties.}
#' \item{premiumMaxTemplateIndex}{Only available for templatized columns. This is the maximum index for the column for premium properties.}
#' }
#'
#' @source \url{https://www.googleapis.com/analytics/v3/metadata/ga/columns?pp=1}
#'
#' @references
#' \href{https://developers.google.com/analytics/devguides/reporting/metadata/v3/}{Google Analytics Metadata API}
#'
#' \href{https://developers.google.com/analytics/devguides/reporting/core/dimsmets}{Core Reporting API - Dimensions & Metrics Reference}
#'
#' @keywords data datasets
#' @docType data
#' @name ga
#'
#' @seealso \code{\link{get_ga}} \code{\link{list_dimsmets}}
#'
#' @examples
#' # a count of parameters types
#' table(ga$type)
#' # parameters groups
#' table(ga$group)
#' # get a deprecated parameters was replaced by
#' subset(ga, status == "DEPRECATED", c(id, replacedBy))
#' # get a calculation metrics
#' subset(ga, !is.na(calculation), c(id, calculation))
#' # get a not deprecated metrics from user group
#' subset(ga, group == "User" & type == "METRIC" & status != "DEPRECATED", id)
#' # get parameters allowed in segments
#' subset(ga, allowedInSegments, id)
NULL
