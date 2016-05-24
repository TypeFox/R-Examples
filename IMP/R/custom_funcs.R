
#' @title Model evaluation measures for Binary classification models
#' @description Generates & plots the following performance evaluation & validation measures for Binary Classification Models - Hosmer Lemeshow goodness of fit tests,
#'   Calibration plots, Lift index & gain charts & concordance-discordance measures
#' @param list_models A list of one (or more) dataframes for each model whose performance is to be evaluated. Each dataframe should comprise of 2 columns with the first column indicating the class labels (0 or 1)
#' and the second column providing the raw predicted probabilities
#' @param g The number of groups for binning. The predicted probabilities are binned as follows
#'
#' For Hosmer-Lemshow (HL) test: Predicted probabilities binned as per g unique quantiles i.e. cut_points = unique(quantile(predicted_prob,seq(0,1,1/g)))
#'
#' For Lift-Index & Gain charts: Same as HL test, however if g > unique(predicted_probability), the predicted probabilities
#' are used as such without binning
#'
#' For calibration plots, g equal sized intervals are created (of width 1/g each)
#'
#' @param sample_size_concord For computing concordance-discordance measures (and c-statistic) a random sample
#' is drawn from each dataset (if nrow(dataset) > 5000). Default sample size of 5000 can be adjusted by changing the value of this
#' argument
#'
#' @param perf_measures Select the required performance evaluation and validation measure/s, from the following
#' options - c('hosmer','calibration','lift','concord'). Default option is All
#'
#' @return A nested list with 2 components - a list of dataframes and a list of plots - containing
#' the outcomes of the different performance evaluations carried out.

#' @export
#' @import dplyr ggplot2
#' @importFrom tidyr gather
#' @importFrom stats as.formula pchisq quantile
#' @examples
#' model_1 <- glm(Species ~ Sepal.Length,data=iris,family=binomial)
#' model_2 <- glm(Species ~ Sepal.Width, data=iris, family = binomial)
#' df1 <- data.frame(model_1$y,fitted(model_1))
#' df2 <- data.frame(model_2$y,fitted(model_2))
#' staticPerfMeasures(list(df1,df2),g=10, perf_measures = c("hosmer","lift"))

staticPerfMeasures <- function(list_models, g,
                               perf_measures = c('hosmer','calibration','lift','concord'),sample_size_concord = 5000) {

    error_handler(list_models,g,method="Performance_Measure")

    perf_functions <- list(hosmer_df = hl_function, hosmer_results = hl_test, calibration = calib_function,
        lift = lift_function, concord = conc_disc)

    perf_functions <- perf_functions[grep(paste(perf_measures,collapse="|"),names(perf_functions))]

    df_out <- lapply(list_models, function(x) lapply(perf_functions, function(f) f(x, g, sample_size_concord)))
    df_out <- lapply(1:length(perf_functions), function(i) combine_df(df_out, i))
    names(df_out) <- names(perf_functions)

    if (length(grep("hosmer_results|concord",names(df_out))) != 0) {

      plot_dfs <- df_out[-grep("hosmer_results|concord",names(df_out))]

    } else {

      plot_dfs <- df_out
    }

    plot_functions <- list(hosmer = plot_HL,calibration = plot_calib,
                           lift = plot_lift)
    plot_functions <- plot_functions[grep(paste(perf_measures,collapse="|"),names(plot_functions))]

    plots_out <- Map(function(f,x) f(x), plot_functions, plot_dfs)

    if ("concord" %in% perf_measures ) {

      plots_out$concord <- plot_condis(list_models)
    }

    return(list(data = df_out, plots = plots_out))

}

################################################################ CONFUSION MATRIX #################################

#' @title Confusion Matrix for Binary Classification Models
#' @description Generates confusion matrix for a specified probability threshold. Also computes the
#' following metrics - Accuracy, True Positive Rate, False Positive Rate & Precision. Multiple models
#' can be passed as arguments to this function
#' @param list_models A list of one (or more) dataframes for each model whose performance is to be evaluated. Each dataframe should comprise of 2 columns with the first column indicating the class labels (0 or 1)
#' and the second column providing the raw predicted probabilities
#' @param t Probability threshold value
#' @param reps Performance measures derived from the confusion matrix (Accuracy, TPR, FPR & Precision) are
#' computed for a range of different probability thresholds. The "reps" argument controls the number of
#' different probability thresholds considered (threshold range given by the sequence - seq(0,1,1/reps))
#' @param reps.all.unique Logical; If set to True, Performance measures are computed for each unique Probability value
#'
#' @return If reps = NULL, the output will be a list with 2 components - a confusion matrix
#' dataframe and a dataframe with the values of the computed metrics (Accuracy,TPR,FPR,Precision). If reps argument is supplied, an
#' additional dataframe containing the metrics values for different probability thresholds is
#' included in the output
#' @export
#' @import dplyr ggplot2
#' @importFrom tidyr gather
#' @importFrom stats as.formula pchisq quantile
#'
#' @examples
#' model_1 <- glm(Species ~ Sepal.Length,data=iris,family=binomial)
#' model_2 <- glm(Species ~ Sepal.Width, data=iris, family = binomial)
#' df1 <- data.frame(model_1$y,fitted(model_1))
#' df2 <- data.frame(model_2$y,fitted(model_2))
#' staticConfMatrix(list(df1,df2),t=0.2)

staticConfMatrix <- function(list_models, t, reps = NULL, reps.all.unique = F) {

  error_handler(list_models,t,method="Confusion")

    if (!is.null(reps)) {
        out_metrics_range <- lapply(list_models, function(x) conf_range(x, reps))
        out_metrics_range <- do.call(rbind, out_metrics_range)
        Model = unlist(lapply(seq_along(list_models), function(i) rep(paste("Model", i), reps + 1)))
        out_metrics_range <- data.frame(Model, out_metrics_range)
    }

   if (reps.all.unique) {

     out_metrics_range <- lapply(list_models, function(x) conf_range(x, reps, all.unique = T))
     out_metrics_range <- do.call(rbind, out_metrics_range)
     Model = unlist(lapply(seq_along(list_models), function(i) rep(paste("Model", i), nrow(unique(list_models[[i]][2])))))
     out_metrics_range <- data.frame(Model, out_metrics_range)
   }

    out_conf <- lapply(list_models, function(x) conf_mat(x, t))
    out_conf <- as.data.frame(do.call(cbind, out_conf))
    out_conf <- data.frame(Pred = c("Pred-0", "Pred-1"), out_conf)

    colnames(out_conf)[-1] <- rep(c('Actuals-0','Actuals-1'),times=(ncol(out_conf)-1)/2)
    rownames(out_conf) <- NULL

    out_metrics <- lapply(list_models, function(x) conf_mat_metrics(x, t))
    out_metrics <- do.call(rbind, out_metrics)
    Model = sapply(seq_along(list_models), function(i) paste("Model", i))
    out_metrics <- data.frame(Model, out_metrics)

    if ((!is.null(reps)) || (reps.all.unique)) {
        return(list(metrics_range = out_metrics_range, conf = out_conf, metrics = out_metrics))
    } else {
        return(list(matrix = out_conf, measures = out_metrics))
    }

}

