#' plotDoseResponseCurve
#'
#' Plots the dose response curve for a given sample id from an IncucyteDRCSet object
#'
#' @param idrc_set IncucyteDRCSet object
#' @param sampleid The sample id to plot
#' @param native deprecated
#'
#' @return a ggplot2 object (if native is FALSE) or NULL but draws to open graphics object (if native is TRUE)
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
#' plotDoseResponseCurve(test_idrc_set, 'PDD00017273', native=FALSE)
plotDoseResponseCurve <- function(idrc_set, sampleid, native=FALSE) {

    #make sure that fitDoseResponseCurve has been run
    if(is.null(idrc_set$drc_models)) {
        stop("Need to run fitDoseResponseCurve first to fit dose response models")
    }

    q_sampleid <- sampleid

    drc_models_filtered <- idrc_set$drc_models %>%
        dplyr::filter(sampleid == q_sampleid)

    drc_model <- drc_models_filtered$drc_model[[1]]
    if(is.null(drc_model)) return(NULL)

    EC50val <- drc::ED(drc_model,50, display=F)[1]

    #extract the generated function from the curve fit and make the curve
    drc_model_func <- drc_model$curve[[1]]
    conc_range <- unique(drc_model$dataList$dose)
    conc_range <- log10(conc_range[conc_range>0])
    conc_seq <- 10^seq(min(conc_range), max(conc_range), 0.01)
    drc_model_curve <- data.frame(conc=conc_seq,
                          value=drc_model_func(conc_seq))

    #do the plot
    p2 <- ggplot(drc_model$origData, aes(y=value, x=conc)) +
        geom_point(alpha=0.8, shape=21, size=rel(3)) +
        ggtitle(sprintf("%s", sampleid)) +
        scale_x_log10() +
        geom_line(data=drc_model_curve, colour='red') +
        geom_vline(xintercept = EC50val, color='red', alpha=0.5, linetype='dashed') +
        cowplot::theme_cowplot()

    return(p2)

}
