#' calculateEC50
#'
#' Calculates the EC50s for the dose response curves in an IncucyteDRCSet object
#'
#' @param idrc_set IncucyteDRCSet object
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
calculateEC50 <- function(idrc_set) {

    #make sure that fitDoseResponseCurve has been run
    if(is.null(idrc_set$drc_models)) {
        stop("Need to run fitDoseResponseCurve first to fit dose response models")
    }

    #extract parameters from models filtering out NULLs
    drc_ec50 <- idrc_set$drc_models %>%
        dplyr::filter(!is.null(drc_model)) %>%
        dplyr::mutate(EC50 = drc::ED(drc_model,50, display=F)[1],
                      SE = drc::ED(drc_model,50, display=F)[2],
                      #R2 <- dr.calc.r2(dr),
                      b = coef(drc_model)[1],
                      c = coef(drc_model)[2],
                      d = coef(drc_model)[3],
                      e = coef(drc_model)[4]) %>%
        dplyr::select(-drc_model) %>%
        dplyr::ungroup()

    #deal with null models
    drc_ec50_nulls <- idrc_set$drc_models %>%
        dplyr::filter(is.null(drc_model)) %>%
        dplyr::mutate(EC50=NA, SE=NA, b=NA, c=NA, d=NA, e=NA) %>%
        dplyr::select(-drc_model) %>%
        dplyr::ungroup()

    #combine outputs
    drc_ec50 <- drc_ec50 %>%
        dplyr::bind_rows(drc_ec50_nulls) %>%
        as.data.frame()

    #sort out output
    output <- idrc_set
    output$drc_ec50 <- drc_ec50

    return(output)

}
