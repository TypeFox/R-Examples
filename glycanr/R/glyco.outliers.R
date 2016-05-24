# Fix the Non-standard evaluation usage for check()
if(getRversion() >= "2.15.1"){
    utils::globalVariables(c("value", "variable", "outlier", "desc"))
}

#' Discover outliers in glycan data
#'
#' Returns outliers within every glycan structure
#'
#' @author Ivo Ugrina
#' @export
#' @importFrom dplyr %>%
#' @importFrom grDevices boxplot.stats
#' @importFrom stats IQR kruskal.test p.adjust
#' @param data data frame in long format containing glycan measurements
#' @param group a possible grouping parameter on which
#'   stratification of \code{data} should be conducted. It should be
#'   a name of one of the columns in dataframe \code{data}
#'   and of type \code{factor}.
#' @param outlier.function a function that checks for outliers in
#'   a vector. Receives one parameter representing a vector and returns
#'   logical vector indicating outliers.
#' @param alpha If outlier.function parameter is set to NULL
#'   outliers are calculated as those points that are greater
#'   than upper quartile + alpha * IQR (interquartile range) or
#'   lower than lower quartile - alpha * IQR (interquartile range).
#'   If parameter outlier.function is not NULL parameter alpha is not used.
#' @details
#' Input data frame should have at least the following three columns: \cr
#'   - gid - representing a unique name of a sample \cr
#'   - glycan - representing glycan names \cr
#'   - value - representing measured values 
#' @return Returns a data.frame with outliers 
#' @examples
#' data(mpiu)
#' 
#' glyco.outliers(mpiu)
#'
#' # outliers per plate
#' glyco.outliers(mpiu, group="Plate")
glyco.outliers <- function(data,
                           group=NULL,
                           outlier.function=NULL,
                           alpha=1.5){
    warning("Version 0.3 of glycanr introduces a change in glyco.outliers function. It expects data frame input in long format now.",
            call. = FALSE)
    
    outf <- function(x) {
        lq <- boxplot.stats(x)$stats[2]
        uq <- boxplot.stats(x)$stats[4]

        iqr <- IQR(x, na.rm=TRUE)
        
        ifelse(x > uq + alpha * iqr | x < lq - alpha * iqr, TRUE, FALSE)
    }
    
    stopifnot(is.data.frame(data))
    if (!is.null(outlier.function)) {
        stopifnot(is.function(outlier.function))
        outf <- outlier.function
    }
    
    if (is.null(group)) {
        X.out <- data %>%
          dplyr::group_by(glycan) %>%
          dplyr::mutate(outlier = outf(value)) %>% 
          dplyr::ungroup()
    } else {
        g <- lapply(c("glycan", group), as.symbol)
        X.out <- data %>%
          dplyr::group_by_(.dots=g) %>%
          dplyr::mutate(outlier = outf(value)) %>% 
          dplyr::ungroup()
    }
    
    dplyr::select(dplyr::filter(X.out, outlier == TRUE), -outlier)
} 
