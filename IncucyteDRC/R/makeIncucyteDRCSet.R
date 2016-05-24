#' makeIncucyteDRCSet
#'
#' Function to construct an IncucyteDRCSet object.  Contains data for just those wells that are from the
#'  same cell line background - ie cell line, passage number, cell number etc are all the same.  Give that the
#'  same cell line background is used, the control wells are common to all wells regardless of different
#'  compounds or concentrations of compound.
#'
#' @param platemap Platemap dataframe containing information for just those wells required
#' @param platedata IncucyteDRCPlateData object that will be filtered according to the wells present in
#'  the platemap
#' @param cut_time The time at which to extract the data for the dose response curve.  Default is NULL.
#' @param metadata A single row data frame containing any information specific to this IncucyteDRCSet object.
#' @param pm_warn Boolean to control platemap check warnings.  Default is TRUE.
#'
#' @return IncucyteDRCSet object
#' @export
#'
#' @examples
#' pm_file <- system.file(file='extdata/example.PlateMap', package='IncucyteDRC')
#' test_pm <- importPlatemapXML(pm_file)
#' data_file <- system.file(file='extdata/example_data.txt', package='IncucyteDRC')
#' test_data <- importIncucyteData(data_file, metric='pc')
#'
#' test_pm_filtered <-  dplyr::filter(test_pm,
#'                                    samptype %in% c('C', 'S') & growthcondition == '8 x 10e4/mL')
#' test_set <- makeIncucyteDRCSet(test_pm_filtered, test_data)
#'
#' print(test_set)


makeIncucyteDRCSet <- function(platemap, platedata, cut_time=NULL, metadata=NULL, pm_warn=TRUE) {

    stopifnot(inherits(platemap, 'data.frame'))
    stopifnot(inherits(platedata, 'IncucyteDRCPlateData'))

    if(is.null(attr(platemap, 'IncucyteDRCPlatemap')) & pm_warn==TRUE) {
        warning('Recommended that platemap data frames are parsed through importPlatemap function to check formatting')
    }

    platedata$data <- platedata$data %>% dplyr::filter(wellid %in% platemap$wellid) %>% as.data.frame()

    outdata <- list(platemap=as.data.frame(platemap),
                    platedata=platedata,
                    cut_time=cut_time,
                    metadata=NULL)

    class(outdata) <- 'IncucyteDRCSet'

    return(outdata)

}
