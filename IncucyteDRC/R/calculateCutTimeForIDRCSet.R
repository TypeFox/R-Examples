#' calculateCutTimeForIDRCSet
#'
#' Uses the control growth curves in an IncucyteDRCSet to calclate the appropriate cut time for a specific
#' number of doublings
#'
#' @param idrc_set IncucyteDRCSet object
#' @param baseline_time The timepoint which forms the baseline for calculting number of doublings
#' @param no_doublings The number of doublings required
#' @param max_val The maximum allowable growth curve value
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
#' plotIncucyteDRCSet(test_idrc_set, grouped=TRUE)
#' test_idrc_set <- calculateCutTimeForIDRCSet(test_idrc_set)
#' plotIncucyteDRCSet(test_idrc_set, grouped=TRUE)
#' print(test_idrc_set$cut_plot)
#'
calculateCutTimeForIDRCSet <- function(idrc_set, baseline_time=24, no_doublings=4, max_val=80) {

    if(is.null(idrc_set$fitted_models_grouped)) {
        stop('Need to fit growth curves first using fitGrowthCurvesGrouped')
    }

    #baseline_time=24; no_doublings=4; max_val=80;

    #platemap grouped
    platemap_grouped <- idrc_set$platemap %>%
        dplyr::select(sampleid, conc, samptype, concunits, group_idx) %>%
        dplyr::distinct()

    #get the spline models for the control data
    control_drc_data <- idrc_set$fitted_models_grouped %>%
        dplyr::inner_join(platemap_grouped, by=c('sampleid', 'conc')) %>%
        dplyr::select(gc_model, sampleid, conc, samptype, concunits, group_idx) %>%
        dplyr::filter(samptype=='C')

    #should only be a single model
    if(nrow(control_drc_data) != 1) {
        stop('There should be a single combination of sampleid/conc marked as controls - please redefine your platemap')
    }

    gcm <- control_drc_data$gc_model[[1]]

    calc_cut_res <- calculateCutTimeForGrowthCurve(gcm, baseline_time, no_doublings, max_val)

    calculated_cut <- calc_cut_res$calculated_cut %>%
        dplyr::mutate(group_idx=control_drc_data$group_idx[1],
                      sampleid=control_drc_data$sampleid[1],
                      conc=control_drc_data$conc[1],
                      concunits=control_drc_data$concunits[1]) %>%
        as.data.frame()

    output <- idrc_set
    output$calculated_cut <- calculated_cut
    output$cut_time <- calc_cut_res$cut_time
    output$cut_plot <- calc_cut_res$cut_plot
    return(output)


}
