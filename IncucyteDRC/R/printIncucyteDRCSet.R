#' print.IncucyteDRCSet
#'
#' Prints information on an IncucyteDRCSet object
#'
#' @param x IncucyteDRCSet object
#' @param ... Additional arguments
#' @export
#'
#' @examples
#' pm_file <- system.file(file='extdata/example.PlateMap', package='IncucyteDRC')
#' test_pm <- importPlatemapXML(pm_file)
#' data_file <- system.file(file='extdata/example_data.txt', package='IncucyteDRC')
#' test_data <- importIncucyteData(data_file, metric='pc')
#' test_list <- splitIncucyteDRCPlateData(test_pm, test_data, group_columns='growthcondition')
#' class(test_list)
#' class(test_list[[2]])
#' print(test_list[[2]])

print.IncucyteDRCSet <- function(x, ...) {

    cat("## This is an IncucyteDRCSet S3 object with the following elements\n")
    cat(names(x))

    cat('\n\n## Cut Time (cut_time)\n')
    cat(x$cut_time)
    cat('\n## Metadata (metadata)\n')
    print(x$metadata)

}
