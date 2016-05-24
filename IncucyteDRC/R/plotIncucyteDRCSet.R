#' plotIncucyteDRCSet
#'
#' Plot the growth curve data from an IncucyteDRCSet object
#'
#' @param idrc_set IncucyteDRCSet object
#' @param grouped Boolean - whether the curves should be grouped or not.
#' @return ggplot object
#' @export
#' @import ggplot2
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
#'
plotIncucyteDRCSet <- function(idrc_set, grouped=FALSE) {

    if(grouped & is.null(idrc_set$fitted_data_grouped)) {
        stop('Need to fit growth curves first using fitGrowthCurvesGrouped')
    }

    if(!grouped & is.null(idrc_set$fitted_data_grouped)) {
        stop('Need to fit growth curves first using fitGrowthCurvesIndividual')
    }


    #combine the platemap and data
    data <- idrc_set$platemap %>% dplyr::inner_join(idrc_set$platedata$data, by='wellid')

    out_plot <- ggplot(data, aes(x=elapsed, y=value, colour=as.factor(round(conc, 3)))) +
                    geom_point(alpha=0.2) +
                    scale_colour_discrete(name='conc') +
                    facet_wrap(~sampleid) + theme_bw()

    if(grouped) {
        out_plot <- out_plot + geom_line(data=idrc_set$fitted_data_grouped, aes(group=paste0(sampleid, conc)))
    } else {
        out_plot <- out_plot + geom_line(data=idrc_set$fitted_data_indiv, aes(group=wellid))
    }

    if(is.numeric(idrc_set$cut_time)) {
        out_plot <- out_plot + geom_vline(xintercept=idrc_set$cut_time, colour='blue', linetype='dashed')
    }

    return(out_plot)

}
