#' fitDoseResponseCurve
#'
#' Fits the dose response curve using drc for an IncucyteDRCSet object
#'
#' @param idrc_set IncucyteDRCSet object
#' @param include_control Whether to include control sample as zero conc control
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
#' test_list <- splitIncucyteDRCPlateData(test_pm, test_data, group_columns='growthcondition')
#'
#' print(test_list)
#'
#' test_idrc_set <- fitGrowthCurvesGrouped(test_list[[2]])
#' test_idrc_set <- fitGrowthCurvesIndividual(test_idrc_set)
#' test_idrc_set <- calculateDRCData(test_idrc_set, cut_time=100)
#' test_idrc_set <- fitDoseResponseCurve(test_idrc_set)
#' test_idrc_set <- calculateEC50(test_idrc_set)
#' exportEC50Data(test_idrc_set)
#'
fitDoseResponseCurve <- function(idrc_set, include_control=FALSE) {

    df <- exportDRCDataToDataFrame(idrc_set, include_control)

    #set up subroutine for fit
    drc_fit <- function(x) {
        m1 <- NULL
        if (length(unique(x$conc)) > 2) {
            try (m1 <- drc::drm(value~conc, fct=drc::LL.4(), data=x))
        }
        return(m1)
    }

    #fit models - derive a data frame
    drc_models <- df %>%
        dplyr::group_by(sampleid) %>%
        dplyr::do(drc_model=drc_fit(.))

    #sort out output
    output <- idrc_set
    output$drc_models <- drc_models

    return(output)

}
