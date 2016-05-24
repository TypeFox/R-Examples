#' importPlatemap
#'
#' Imports a platemap configuration from a tab delimited file or dataframe.
#'
#' @param input Either a path to a text file or a data frame
#' @param control_cpd Specify the compound to use as baseline.  Defaults to DMSO
#'
#' @return IncucyteDRCPlateMap object
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' pm_file <- system.file(file='extdata/example_platemap.txt', package='IncucyteDRC')
#' test_pm <- importPlatemap(pm_file)
#' head(test_pm)
#' test_pm_df <- importPlatemap(as.data.frame(test_pm))
#' head(test_pm_df)
importPlatemap <- function(input, control_cpd='DMSO') {

    if(is.data.frame(input)) {
        message("Using data frame as an input")
        platemap_df <- input
    } else if(file.exists(input)) {
        message(sprintf("Importing platemap text file from %s",input))
        platemap_df <- read.table(input, header=TRUE, sep='\t', stringsAsFactors = FALSE)
    } else {
        stop('Require either a data frame or valid file path as input')
    }

    required_colnames <- c("row", "col", "sampleid", "conc", "samptype", "concunits",
                           "growthcondition", "celltype", "passage", "seedingdensity")

    if(!all(required_colnames %in% colnames(platemap_df))) {
        stop(sprintf("Input file should be tab delimited and contain all of the following columns: %s", paste(required_colnames, collapse=', ')))
    }

    if(!('wellid' %in% colnames(platemap_df))) {
        platemap_df <- platemap_df %>%
            dplyr::mutate(wellid = paste0(LETTERS[row], col)) %>%
            as.data.frame()
    }

    attr(platemap_df, 'IncucyteDRCPlatemap') <- TRUE
    #class(platemap_df) <- c('IncucyteDRCPlatemap', class(platemap_df))

    message('Plate map import successful!')
    return (platemap_df)

}
