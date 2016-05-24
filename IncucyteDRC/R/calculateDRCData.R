#' calculateDRCData
#'
#' For a given cut time, calculate values from the growth curves in a IncucyteDRCSet object
#'
#' @param idrc_set IncucyteDRCSet object
#' @param cut_time Desired cut time.  If NULL will use the cut_time in the IncucyteDRCSet object (if set)
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
#' plotIncucyteDRCSet(test_idrc_set, grouped=FALSE)
#' plotIncucyteDRCSet(test_idrc_set, grouped=TRUE)
#' test_idrc_set <- calculateDRCData(test_idrc_set, cut_time=100)
#' print(test_idrc_set)
#' test_idrc_set$drc_data
#' plotIncucyteDRCSet(test_idrc_set)
#'
calculateDRCData <- function(idrc_set, cut_time=NULL) {

    if(is.null(idrc_set$cut_time) & is.null(cut_time)) {
        stop('Set cut_time in IncucyteDRCSet or function parameter')
    } else if (is.null(idrc_set$cut_time) & is.numeric(cut_time)) {
        message('Using cut time provided to function')
    } else if (is.numeric(idrc_set$cut_time) & is.null(cut_time)) {
        message('Using cut time from IncucyteDRCSet object')
        cut_time <- idrc_set$cut_time
    } else if (is.numeric(idrc_set$cut_time) & is.numeric(cut_time)) {
        warning('Using cut time provided to function and overriding cut time from IncucyteDRCSet object')
    } else {
        stop('Provide cut time as a numeric parameter')
    }

    if(is.null(idrc_set$fitted_models_indiv)) {
        stop('Need to fit growth curves first using fitGrowthCurvesIndividual')
    }

    #predict the value at the cut time using the splines
    drc_data <- idrc_set$fitted_models_indiv %>%
        dplyr::mutate(cut_val=predict(gc_model,cut_time),
               cut_time=cut_time) %>%
        dplyr::ungroup() %>%
        dplyr::inner_join(idrc_set$platemap, by='wellid') %>%
        dplyr::select(wellid, sampleid, conc, samptype, concunits, cut_time, cut_val) %>%
        as.data.frame()

    #construct the output object
    output <- idrc_set
    output$drc_data <- drc_data
    output$cut_time <- cut_time
    return(output)

}
