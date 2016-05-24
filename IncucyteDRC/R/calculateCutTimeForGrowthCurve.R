#' calculateCutTimeForGrowthCurve
#'
#' For a given growth curve model, calculates the appropriate cut time for a specific
#' number of doublings
#'
#' @param gcm Growth curve model
#' @param baseline_time The timepoint which forms the baseline for calculting number of doublings
#' @param no_doublings The number of doublings required
#' @param max_val The maximum allowable growth curve value
#'
#' @return list object containing a single row data frame with the results and a plotted representation
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
#' test_idrc_set$fitted_models_grouped
#' test_gcm <- test_idrc_set$fitted_models_grouped$gc_model[[1]]
#' test_gcm
#' calculateCutTimeForGrowthCurve(test_gcm, baseline_time=24, no_doublings=4, max_val=80)
#'
calculateCutTimeForGrowthCurve <- function(gcm, baseline_time=24, no_doublings=4, max_val=80) {

    stopifnot(class(gcm) == 'loess')

    fitted_df <- data.frame(elapsed=seq(0, max(gcm$x), 1))
    fitted_df$value = predict(gcm, fitted_df, se=TRUE)$fit
    fitted_df$diff2 <- c(NA, NA,  diff(fitted_df$value, 1, 2) )  #generate second differences
    fitted_df$ma2 <- stats::filter(fitted_df$diff2,rep(1/30,30),sides=2) #moving average
    if (which.min(fitted_df$ma2) > which.max(fitted_df$ma2)) {  #if the minimum occurs after the maximum, ie typical situation where rate of growth fastest before tails off
        time_neg_inflex <- fitted_df [ which.min(fitted_df$ma2) , 'elapsed' ]  #get the timepoint where slowing of growth is happening fastest
        max_control <- min (predict(gcm , time_neg_inflex) ,  max(fitted_df$value) ) #dont go higher than absolute max control or value at time_neg_inflex
    } else {
        max_control <- max(fitted_df$value) #if the inflexion min occurs before the max this is a slow grower so don't apply inflection point criterion (just return max pc)
    }

    fitted_df <- fitted_df [ 1:which.max(fitted_df$value) , ] #get rid of datapoints after the max value
    val_at_baseline <- predict(gcm , baseline_time) #predict the value at baseline
    val_post_ndoublings <- min ( val_at_baseline * 2^no_doublings , max_val, max_control ) #work out the value post doublings, make sure it's not more than max spec val or max control val
    time_post_ndoublings <- fitted_df [ which.min(abs(fitted_df$value - val_post_ndoublings)) , 'elapsed' ]  #which timepoint closest to desired value
    cut_time <- round(time_post_ndoublings,0)
    actual_doublings <- log2 ( val_post_ndoublings / val_at_baseline )


    calculated_cut <- data.frame(baseline_time=baseline_time,
                                 no_doublings=no_doublings,
                                 max_val=max_val,
                                 val_at_baseline=val_at_baseline,
                                 val_post_ndoublings=val_post_ndoublings,
                                 time_post_ndoublings=time_post_ndoublings,
                                 cut_time=cut_time,
                                 actual_doublings=actual_doublings,
                                 stringsAsFactors=FALSE
    )


    p1 <- ggplot(fitted_df, aes(elapsed, value)) + geom_point(colour='red') +
        geom_vline(xintercept=cut_time, colour='blue', linetype='dashed') +
        geom_vline(xintercept=baseline_time, linetype='dotted') +
        geom_hline(yintercept=c(val_at_baseline, max_val), linetype='dotted') +
        geom_hline(yintercept=val_post_ndoublings, colour='blue', linetype='dashed') +
        cowplot::theme_cowplot()

    p2 <- ggplot(fitted_df, aes(elapsed, ma2)) + geom_point(colour='green') +
        geom_vline(xintercept=cut_time, colour='blue', linetype='dashed') +
        cowplot::theme_cowplot()
    cut_plot <- cowplot::plot_grid(p1,p2, align='v', nrow=2, rel_heights = c(2,1))
    cut_plot

    return(list(calculated_cut=calculated_cut,
                cut_plot=cut_plot,
                cut_time=cut_time
           ))


}


