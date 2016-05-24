utils::globalVariables(c("concentration", "lower_bound", "upper_bound",
                         "median_activity", "activity", "model_num", 
                         "median_slope"))

#bootstrap_metaregression function --------------------------------------
#' Perform Bootstrap Metaregression
#' 
#' 
#' Performs bootstrap metaregression on a concentration-response dataset.
#' 
#' This function performs bootstrap metaregression on a concentration-response
#' dataset. The dataset must be a data.frame with two columns: 
#' 1) Activity and 2) Concentration.
#' 
#' @param x an object of class \code{data.frame}.
#' 
#' @param dataset_size a numeric object with the size of an individual dataset
#' in \code{x}.
#' 
#' @param iterations the number of iterations to run; default is 1,000.
#' 
#' @return bmr_obj a \code{bmr} object that holds all of the bootstrap 
#' metaregression models produced.
#' 
#' @examples
#' bmr_obj <- bootstrap_metaregression(oxybenzone, 15, 100)
#' 
#' @import splines
#' @importFrom stats lm
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom stats predict
#' @importFrom stats quantile
#' 
#' @export
bootstrap_metaregression <- function(x, dataset_size, iterations=1000){
  model_fits <- list()
  bootstrap_smr_models <- list()
  for(i in 1:iterations){
    bootstrap_sample_rows <- sample(c(1:length(x$Activity)), size=dataset_size, replace=TRUE)
    bootstrap_data <- x[bootstrap_sample_rows, c(3,4)]
    fm1 <- lm(Activity ~ ns(Concentration), data=bootstrap_data)
    bootstrap_smr_models[[i]] <- fm1
    predict_fits <- predict(fm1, data.frame(Concentration=bootstrap_data$Concentration))
    names(predict_fits) <- bootstrap_data$Concentration
    model_fits[[i]] <- predict_fits
  }
  #Putting the models and model fits into the bootstrap metaregression class object
  bmr_obj <- bmr(models = bootstrap_smr_models, fits = model_fits)
  bmr_obj <- calculate_confidence_and_median(bmr_obj)
  return(bmr_obj)
}

#---------------------------------------------------------------------

#plot metaregression spaghetti plot function--------------------------------------
#' Make a spaghetti plot for the metaregression results
#' 
#' 
#' Plots a subset of the metaregression results.
#' 
#' This function plots the concentration-response curves for a subset of the 
#' metaregression models generated.
#' 
#' @param bootstrap_metaregression_obj the object that contains the bootstrap 
#' metaregression models as a \code{bmr} object.
#' 
#' @param number_to_plot the number of bootstrap metaregression 
#' concentration-response models to plot. The default is 100 models.
#' 
#' @examples
#' bmr_obj <- bootstrap_metaregression(oxybenzone, 15, 100)
#' plot_metaregression_spaghetti_plot(bmr_obj, number_to_plot=40)
#' 
#' @import splines
#' 
#' @export
plot_metaregression_spaghetti_plot <- function(bootstrap_metaregression_obj, number_to_plot=100){
  model_fits <- bootstrap_metaregression_obj@fits
  results_df <- NULL
  if(number_to_plot < length(model_fits)){
    random_sample <- sample(1:length(model_fits), size=number_to_plot, replace=FALSE)
    for(i in random_sample){
      temp_df <- data.frame(model_num = i, concentration = as.numeric(names(model_fits[[i]])), activity = model_fits[[i]])
      results_df <- rbind(results_df, temp_df)
    }
  }
  else{
    for(i in 1:length(model_fits)){
      temp_df <- data.frame(model_num = i, concentration = as.numeric(names(model_fits[[i]])), activity = model_fits[[i]])
      results_df <- rbind(results_df, temp_df)
    }
  }
  concentration_response_plot <- ggplot(results_df, aes(x=log10(concentration), y=activity, group=model_num))
  concentration_response_plot + 
    geom_point(color="black") + 
    geom_line(color="red")
}

#---------------------------------------------------------------------

#calculate confidence and median function--------------------------------------
#' Calculate the confidence envelope and the median for the bootstrap
#' metaregression concentration response models.
#' 
#' 
#' Calculates the 95% confidence envelope and the median for
#' concentration-response data.
#' 
#' This is an internal function that calculates the 95% confidence envelope and 
#' the median for concentration-response data. It is called by the 
#' bootstrap_metaregression function.
#' 
#' @param bootstrap_metaregression_obj the object that contains the bootstrap 
#' metaregression models as a \code{bmr} object.
#' 
#' @return bmr_obj a \code{bmr} object that holds all of the bootstrap 
#' metaregression models produced.
#' 
#' @examples
#' bmr_obj <- bootstrap_metaregression(oxybenzone, 15, 100)
#' 
#' @import splines
#' 
#' @export
calculate_confidence_and_median <- function(bootstrap_metaregression_obj){
  model_fits <- bootstrap_metaregression_obj@fits
  confidence_envelope <- NULL
  median_activity <- NULL
  temp_df <- NULL
  for(i in 1:length(model_fits)){
    temp_df_2 <- data.frame(model_num = i, concentration = as.numeric(names(model_fits[[i]])), activity = model_fits[[i]])
    temp_df <- rbind(temp_df, temp_df_2)
  }
  
  for(i in unique(temp_df$concentration)){
    temp_df_sub <- temp_df[which(temp_df$concentration == i), ]
    confidence_envelope_temp <- quantile(temp_df_sub$activity, probs=c(0.025, 0.975))
    median_activity <- rbind(median_activity, data.frame(concentration = i, median_activity = median(temp_df_sub$activity)))
    confidence_envelope_temp_df <- data.frame(concentration = i, lower_bound = confidence_envelope_temp[1], upper_bound = confidence_envelope_temp[2])
    confidence_envelope <- rbind(confidence_envelope, confidence_envelope_temp_df)
  }
  
  bootstrap_metaregression_obj@confidence_envelope <- confidence_envelope
  bootstrap_metaregression_obj@medians <- median_activity
  return(bootstrap_metaregression_obj)
}

#---------------------------------------------------------------------

#plot metaregression confidence envelope function--------------------------------------
#' Plot the metaregression confidence envelope and median results from the 
#' bootstrap metaregression models.
#' 
#' 
#' A function to plot the metaregression confidence envelope and median results from the 
#' bootstrap metaregression models.
#' 
#' @param bootstrap_metaregression_obj the object that contains the bootstrap 
#' metaregression models as a \code{bmr} object.
#' 
#' @param graph_pod a \code{boolean} that determines if the point of departure
#' will be displayed on the graph.
#' 
#' @param pod the chemical's point of departure as a \code{numeric} value
#' 
#' @param pod_threshold the threshold value used to calculate the chemical's
#' point of departure.
#' 
#' @param median_line_color the color for the median line, default is "orange".
#' 
#' @param pod_and_threshold_color the color of the POD and threshold
#' "crosshairs" on the plot. The default is "green".
#' 
#' @examples
#' bmr_obj <- bootstrap_metaregression(oxybenzone, 15, 100)
#' slope_pod <- slope_pod_analysis(bmr_obj, 0.0001, 10, 0.1)
#' pod_and_threshold <- pod_envelope_analysis(bmr_obj, slope_pod, 10, 
#'   min(oxybenzone$Concentration), max(oxybenzone$Concentration), 0.1)
#' plot_metaregression_confidence_envelope(bmr_obj, graph_pod = TRUE, 
#'   pod = pod_and_threshold$pod, pod_threshold=pod_and_threshold$threshold)
#' 
#' @import splines
#' @import ggplot2
#' 
#' @export
plot_metaregression_confidence_envelope <- function(bootstrap_metaregression_obj, graph_pod = FALSE, pod, pod_threshold, median_line_color = "orange", pod_and_threshold_color = "green"){
  if(!graph_pod){
    conf_env_plot <- ggplot(bootstrap_metaregression_obj@confidence_envelope, aes(x=log10(concentration)))
    conf_env_plot + 
      geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound)) +
      geom_line(data=bootstrap_metaregression_obj@medians, aes(x=log10(concentration), y=median_activity), color="orange")
  }
  else{
    conf_env_plot <- ggplot(bootstrap_metaregression_obj@confidence_envelope, aes(x=log10(concentration)))
    conf_env_plot + 
      geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound)) +
      geom_line(data=bootstrap_metaregression_obj@medians, aes(x=log10(concentration), y=median_activity), color=median_line_color) +
      geom_hline(yintercept = pod_threshold, color=pod_and_threshold_color) +
      geom_vline(xintercept = log10(pod), color = pod_and_threshold_color)
  }
}

#---------------------------------------------------------------------

#slope POD analysis function--------------------------------------
#' Slope-based POD analysis
#' 
#' This requires lower and upper limits to be specified. This is the function 
#' that calculates the slope as part of the basis for the POD. The slope is
#' used to identify the lower bound asymptote on the concentration-response
#' curve.
#' 
#' @param bootstrap_metaregression_obj the object that contains the bootstrap 
#' metaregression models as a \code{bmr} object.
#' 
#' @param lower_interpolation_range a \code{numeric} value where the 
#' interpolation should be bounded on the lower end.
#' 
#' @param upper_interpolation_range a \code{numeric} value where the
#' interpolation should be bounded on the upper end.
#' 
#' @param interval_size a \code{numeric} value that specifies how large the
#' interval should be between each value used for interpolation between
#' the lower and upper bounds.
#' 
#' @return slope_pod a two-column \code{data.frame} object that contains the
#' concentration (column 1) and the median slope.
#' 
#' @examples
#' bmr_obj <- bootstrap_metaregression(oxybenzone, 15, 100)
#' slope_pod <- slope_pod_analysis(bmr_obj, 0.0001, 10, 0.1)
#' 
#' @import splines
#' 
#' @export
slope_pod_analysis <- function(bootstrap_metaregression_obj, lower_interpolation_range, upper_interpolation_range, interval_size){
  all_slopes <- NULL
  median_slope <- NULL
  models <- bootstrap_metaregression_obj@models
  interval <- seq(from = lower_interpolation_range, to = upper_interpolation_range, by=interval_size)
  for(i in 1:length(models)){
    model_fits <- predict(models[[i]], data.frame(Concentration = interval))
    slopes <- diff(model_fits, lag=1)/diff(log10(interval), lag=1)
    slope_here <- cbind(Concentration = interval[-1], slope = slopes, model_group = rep(i, length(interval[-1])))
    all_slopes <- rbind(all_slopes, slope_here)
  }
  all_slopes <- as.data.frame(all_slopes, row.names=seq(1, length(all_slopes[,1])))
  
  for(i in unique(all_slopes$Concentration)){
    all_slopes_sub <- all_slopes[which(all_slopes$Concentration == i), ]
    median_slope <- rbind(median_slope, data.frame(concentration = i, median_slope = median(all_slopes_sub$slope)))
  }
  
  return(median_slope)
}

#---------------------------------------------------------------------

#plot slope analysis function--------------------------------------
#' Plot the median slope
#' 
#' This simply plots the slope as a function of concentration.
#' 
#' @param pod_slope_data the \code{data.frame} object that contains the 
#' concentration and slope data.
#' 
#' @param yaxis_limit a \code{boolean} value (default is FALSE) that identifies
#' if the user wants to specify y-axis limits. 
#' 
#' @param yaxis_limit_values a two-element \code{vector} that specifies the
#' y-axis limits. For instance \code{yaxis_limit_values = c(0, 20)}.
#' 
#' @examples
#' bmr_obj <- bootstrap_metaregression(oxybenzone, 15, 100)
#' slope_pod <- slope_pod_analysis(bmr_obj, 0.0001, 10, 0.1)
#' plot_slope_analysis(slope_pod, TRUE, c(0,30))
#' 
#' @import splines
#' @import ggplot2
#' 
#' @export
plot_slope_analysis <- function(pod_slope_data, yaxis_limit = FALSE, yaxis_limit_values){
  concentration_response_plot <- ggplot(pod_slope_data, aes(x=log10(concentration), y=median_slope))
  if(yaxis_limit){
    concentration_response_plot + 
      geom_point(color="black") + 
      geom_line(color="red") +
      ylim(yaxis_limit_values)
  }
  else{
    concentration_response_plot + 
      geom_point(color="black") + 
      geom_line(color="red")
  }
  
}

#---------------------------------------------------------------------

#POD slope analysis function--------------------------------------
#' This is intended as an internal function to facilitate the identification
#' of the concentration where the asymptote is likely to end based on a slope
#' threshold, so that the mean of the upper confidence limit for the
#' asymptote can be calculated.
#' 
#' An internal function for calculating the boundary of the lower asymptote.
#' 
#' @param pod_slope_data the \code{data.frame} object that contains the 
#' concentration and slope data.
#' 
#' @param slope_threshold the \code{numeric} object that sets the threshold
#' for the slope. This determines the upper bound on the concentration range
#' that is determined to be the asymptote. In other words, the asymptote is
#' defined as that region that has a slope less than the threshold at the lower
#' end of the concentration-response curve.
#' 
#' @import splines
#' 
#' @export
pod_slope_analysis <- function(pod_slope_data, slope_threshold = 10){
  pod_slope_data <- na.omit(pod_slope_data)
  thresholded_slope <- pod_slope_data[which(pod_slope_data$median_slope < slope_threshold), ]
  pod <- max(na.omit(thresholded_slope$concentration))
  return(pod)
}

#---------------------------------------------------------------------

#internal model fits function--------------------------------------
#' An internal function for calculating interpolation values.
#' 
#' An internal function for calculating interpolation values.
#' 
#' @param bmr_model a model of class \code{lm} based on a bootstrap natural 
#' splin metaregression of concentration-response data.
#' 
#' @param interval a \code{numeric} value that specifies how large the
#' interval should be between each value used for interpolation.
#' 
#' @return temp_df a \code{data.frame} consisting of concentration and 
#' activity columns that represent the interpolated model-based 
#' response/activity values.
#' 
#' @import splines
#' 
#' @export
internal_model_fits <- function(bmr_model, interval){
  model_fits <- predict(bmr_model, data.frame(Concentration = interval))
  temp_df <- data.frame(concentration = interval, activity = model_fits)
  return(temp_df)
}

#---------------------------------------------------------------------

#POD envelope analysis function--------------------------------------
#' This function calculates the chemical's point of departure.
#' 
#' This function calculates the chemical's point of departure based on the
#' concentration-response data.
#' 
#' @param bmr_obj a \code{bmr} object that holds all of the bootstrap 
#' metaregression models produced.
#' 
#' @param slope_data the \code{data.frame} object that contains the 
#' concentration and slope data.
#' 
#' @param slope_threshold the \code{numeric} object that sets the threshold
#' for the slope. This determines the upper bound on the concentration range
#' that is determined to be the asymptote. In other words, the asymptote is
#' defined as that region that has a slope less than the threshold at the lower
#' end of the concentration-response curve.
#' 
#' @param lower_interpolation_range a \code{numeric} value where the 
#' interpolation should be bounded on the lower end.
#' 
#' @param upper_interpolation_range a \code{numeric} value where the
#' interpolation should be bounded on the upper end.
#' 
#' @param interval_size a \code{numeric} value that specifies how large the
#' interval should be between each value used for interpolation between
#' the lower and upper bounds.
#' 
#' @param agonist_assay a \code{boolean} value that specifies if the assay
#' is an agonist or antagonist assay.
#' 
#' @return a two column \code{data.frame} that contains the chemical's point
#' of departure and the threshold value.
#' 
#' @examples
#' bmr_obj <- bootstrap_metaregression(oxybenzone, 15, 100)
#' slope_pod <- slope_pod_analysis(bmr_obj, 0.0001, 10, 0.1)
#' pod_and_threshold <- pod_envelope_analysis(bmr_obj, slope_pod, 
#'   slope_threshold = 10, min(oxybenzone$Concentration), 
#'   max(oxybenzone$Concentration), interval_size = 0.1)
#' 
#' @import splines
#' @importFrom plyr ldply
#' 
#' @export
pod_envelope_analysis <- function(bmr_obj, slope_data, slope_threshold = 1.0, lower_interpolation_range, upper_interpolation_range, interval_size, agonist_assay=TRUE){
  pod <- pod_slope_analysis(slope_data, slope_threshold)
  confidence_envelope <- bmr_obj@confidence_envelope[which(bmr_obj@confidence_envelope$concentration < pod),]
  if(agonist_assay){
    mean_upper_bound <- mean(confidence_envelope$upper_bound)
    interval <- seq(from = lower_interpolation_range, to = upper_interpolation_range, by=interval_size)
    model_fit_aggregates <- NULL
    median_activity <- NULL
    
    model_fit_aggregates <- ldply(bmr_obj@models, internal_model_fits, interval)
    
    for(i in unique(model_fit_aggregates$concentration)){
      temp_df_sub <- model_fit_aggregates[which(model_fit_aggregates$concentration == i), ]
      median_activity <- rbind(median_activity, data.frame(concentration = i, median_activity = median(temp_df_sub$activity)))
    }
    #Diagnostic to print the upper bound that needs to intersect with the median line to identify the POD
    #print(mean_upper_bound)
    #print(median_activity[which(median_activity$median_activity < mean_upper_bound),])
    return(data.frame(pod = median_activity[max(which(median_activity$median_activity < mean_upper_bound)),1], threshold = mean_upper_bound))
  }
  else{
    mean_lower_bound <- mean(confidence_envelope$lower_bound)
    interval <- seq(from = lower_interpolation_range, to = upper_interpolation_range, by=interval_size)
    model_fit_aggregates <- NULL
    median_activity <- NULL
    
    model_fit_aggregates <- ldply(bmr_obj@models, internal_model_fits, interval)
    
    for(i in unique(model_fit_aggregates$concentration)){
      temp_df_sub <- model_fit_aggregates[which(model_fit_aggregates$concentration == i), ]
      median_activity <- rbind(median_activity, data.frame(concentration = i, median_activity = median(temp_df_sub$activity)))
    }
    #Diagnostic to print the upper bound that needs to intersect with the median line to identify the POD
    #print(mean_upper_bound)
    #print(median_activity[which(median_activity$median_activity < mean_upper_bound),])
    return(data.frame(pod = median_activity[min(which(median_activity$median_activity < mean_lower_bound)),1], threshold = mean_lower_bound))
  }
  
}

#---------------------------------------------------------------------

#' High Throughput Screening Data (Tox21) for Assessing the Estrogenicity of
#' Oxybenzone
#' 
#' A dataset containing the concentration-response data for analyzing the
#' estrogenicity of oxybenzone from PubChem (Assay ID: 743075, Substance ID: 
#' 144209183, Chemical ID: 4632; Assay ID: 743079, Substance ID: 144203969
#' Chemical ID: 4632). 125 rows. 58 rows x 4 cols/variables.
#' 
#' AssayDataset DatasetReplicate Concentration  Activity
#' 
#' \itemize{
#'   \item AssayDataset. Coded value (1/2) that corresponds to an Assay ID  
#'   \item DatasetReplicate. This is the "biological" replicate within an 
#'   AssayDataset (1--3)
#'   \item Concentration. This is the concentration of the chemical in the assay (5.55694e-04--7.03500e+01) 
#'   \item Activity. This is the assay activity. Percent response for these assays. (1.646933e-03--1.140800e+02)
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name oxybenzone
#' @usage data(oxybenzone)
#' @format A data frame with 58 rows and 4 variables
NULL
