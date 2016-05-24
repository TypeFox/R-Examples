#' fitGrowthCurvesGrouped
#'
#' Function to fit loess function to the growth curve data in an IncucyteDRCSet object
#'
#' @param idrc_set IncucyteDRCSet object
#'
#' @return IncucyteDRCSet object
#' @importFrom stats loess loess.control predict coef
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
#'
fitGrowthCurvesGrouped <- function(idrc_set) {

    #combine the platemap and data
    data <- idrc_set$platemap %>% dplyr::inner_join(idrc_set$platedata$data, by='wellid')

    #platemap grouped
    platemap_grouped <- idrc_set$platemap %>%
        dplyr::select(sampleid, conc, samptype, concunits) %>%
        dplyr::distinct()

    #fit the splines
    fitted_models <- data %>%
        dplyr::group_by(sampleid, conc) %>%
        dplyr::do(gc_model=loess (value ~ elapsed  , ., span=0.3, control=loess.control(surface='direct')))

    #establish the data range
    data_range <- seq(from=min(data$elapsed, na.rm=TRUE), to=max(data$elapsed, na.rm = TRUE), by = 1)

    #generate the fitted data for plotting
    fitted_data <- fitted_models %>%
        dplyr::do(data.frame(value=predict(.$gc_model, data_range), elapsed=data_range,
                             sampleid=.$sampleid, conc=.$conc, stringsAsFactors=FALSE)) %>%
        dplyr::ungroup() %>%
        dplyr::inner_join(platemap_grouped, by=c('sampleid', 'conc')) %>%
        as.data.frame()

    #construct the output object
    output <- idrc_set
    output$fitted_data_grouped <- fitted_data
    output$fitted_models_grouped <- fitted_models

    return(output)


}
