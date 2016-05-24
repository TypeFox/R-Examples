#' Return the best VAR models found for a time series data set
#'
#' This function evaluates possible VAR models for the given time series data set and returns a sorted list of the best models found. The first item in this list is the "best model" found.
#'
#' AutovarCore evaluates eight kinds of models: models with and without log transforming the data, lag 1 and lag 2 models, and with and without day dummy variables. For each of these 8 model configurations, we evaluate all possible combinations for including outlier dummies (at 2.5x the standard deviation of the residuals) and retain the best model (the procedure for selecting the best model is described in more detail below).
#'
#' These eight models are further reduced to four models because AutovarCore determines whether adding daydummies improves the model fit (considering only significance bucket and the AIC/BIC score, NOT the number of outlier dummy columns), and only when the best model found with day dummies is a "better model" than the best model found without day dummies (other parameters kept the same), do we include the model with daydummies and discard the one without day dummies. Otherwise, we include only the model without daydummies and discard the one with daydummies.
#'
#' Thus, AutovarCore always returns four models (assuming that we find enough models that pass the Eigenvalue stability test: models that do not pass this test are immediately discarded). There are three points in the code where we determine the best model, which is done according to what we refer to as algorithm A or algorithm B, which we explain below.
#'
#' When evaluating all possible combinations of outlier dummies for otherwise identical model configurations, we use algorithm A to determine the best model. When comparing whether the best model found without day dummy columns is better than the best model found with day dummy columns, we use algorithm B. For sorting the two models with differing lag but both being either logtransformed or not, we again use algorithm A. Then at the end we merge the two lists of two models (with and without logtransform) to obtain the final list of four models. The sorting comparison here uses algorithm B.
#'
#' The reason for the different sorting algorithms is that in some cases we want to select the model with the fewest outlier dummy columns (i.e., the model that retains most of the original data), while in other cases we know that a certain operation (such as adding day dummies or logtransforming the data set) will affect the amount of day dummies in the model and so a fair comparison would exclude this property.
#'
#' Algorithm A applies the following rules for comparing two models in order: \enumerate{
#' \item We first consider the significance bucket. If the two models being compared are in different significance buckets, choose the one with the highest significance bucket, otherwise proceed to step 2.
#'
#' The significance buckets are formed between each of the (decreasingly sorted) specified \code{significance_levels} in the parameters to the autovar function call. For example, if the \code{signifance_levels} are \code{c(0.05, 0.01, 0.005)}, then the significance buckets are \code{(0.05 <= x), (0.01 <= x < 0.05), (0.005 <= x < 0.01),} and \code{(x < 0.005)}. The metric used to place a model into a bucket is the maximum p-level that can be chosen as cut-off for determining whether an outcome is statistically significant such that all residual tests will still pass ("pass" meaning not invalidating the assumption that the residuals are normally distributed). In other words: it is the minimum p-value of all three residual tests of all endogenous variables in the model.
#' \item If the two models being compared are in the same significance bucket, the number of outlier columns is most important. If the two models being compared have a different amount of outlier columns, choose the one with the least amount of outlier columns, otherwise proceed to step 4.
#'
#' For this count of outlier columns, the following rules apply: \itemize{
#' \item Day-part dummies do not add to the count. This is because when they are included, they are included for each model and thus never have any discriminatory power.
#' \item Day dummies count as one outlier column in total (so including day dummies will add one outlier column). This is because we do not necessarily want to punish models if the data happens to exhibit weekly cyclicity, but models that do not need the day dummies and are equally likely should be preferred.
#'
#' \item Outlier dummy variables are split up such that each time point that is considered an outlier has its own dummy outlier variable and adds one to the count of outlier columns. The outliers are, for each variable, the measurements at >2.5 times the standard deviation away from the mean of the residuals or of the squared residuals. Checks are in place to ensure that a time point identified as an outlier by multiple variables only adds a single dummy outlier column to the equation.
#' }
#' \item When the bucket and number of outlier columns for the two models being compared are the same, select the one with the lowest AIC/BIC score. Whether the AIC or BIC is used here depends on the  \code{criterion} option specified in the parameters to the autovar function call.
#' }
#'
#'
#' In the end, we should have one best logtransformed model and one best nonlogtransformed model. We then compare these two models in the same way as we have compared all other models up to this point with one exception: we do not compare the number of outlier columns. Comparing the number of outliers would have likely favored logtransformed models over models without logtransform, as logtransformations typically have the effect of reducing the outliers of a sample.
#'
#' Algorithm B is identical to algorithm A, except that we skip the step comparing the number of outlier dummy variables. Thus, we instead compare by bucket first and AIC/BIC score second. Notice that, if we may assume that the presence or absence of day dummies does not vary between the four models for any particular invocation of the autovar method (which is not an unreasonable assumption to make), that then the arbitrary choice of letting all daydummy columns together add one to the outlier count does not matter at all, since the only times where we are comparing the outlier dummy counts is when both models either both have or both do not have day dummy columns.
#'
#' We are able to compare the AIC/BIC scores of logtransformed and nonlogtransformed models fairly because we compensate the AIC/BIC scores to account for the effect of the logtransformation. We compensate for the logtransformation by adjusting the loglikelihood score of the logtransformed models in the calculation of their AIC/BIC scores (by subtracting the sum of the logtransformed data).
#' @param raw_dataframe The raw, unimputed data frame. This can include columns other than the \code{selected_column_names}, as those may be helpful for the imputation.
#' @param selected_column_names The endogenous variables in the models, specified as a vector of character strings. This argument is required. The selected column names should be a subset of the column names of \code{raw_dataframe}.
#' @param significance_levels A vector with descending p values that indicate cut-offs placing models in different buckets. If it is not specified, this parameter defaults to \code{c(0.05, 0.01, 0.005)}. For example, with the default configuration, a model whose worst (lowest) p-level for any test is 0.03 is always seen as a better model than one whose worst p-level for any test is 0.009, no matter the AIC/BIC score of that model. Also, the lowest significance level indicates the minimum p-level for any test of a valid model. Thus, if a test for a model has a lower p-level than the minimum specified significance level, it is considered invalid.
#' @param test_names The residual tests that should be performed, specified as a vector of character strings. If not specified, this parameter defaults to \code{c('portmanteau', 'portmanteau_squared', 'skewness')}. The possible tests are \code{c('portmanteau', 'portmanteau_squared', 'skewness', 'kurtosis', 'joint_sktest')}. In addition to the residual tests, please note that the Eigenvalue stability test is always performed.
#' @param criterion The information criterion used to sort the models. Valid options are 'AIC' (the default) or 'BIC'.
#' @param imputation_iterations The number of times we average over one Amelia call for imputing the data set. Since one Amelia call averages over five imputations on its own, the actual number of imputations is five times the number specified here. The default value for this parameter is \code{30}.
#' @param measurements_per_day The number of measurements per day in the time series data. The default value for this parameter is \code{1}. If this value is \code{0}, then daypart- and day-dummies variables are not included for any models.
#' @return A sorted list of "valid" models. A "model" is a list with the properties \code{logtransformed}, \code{lag}, \code{varest}, \code{model_score}, \code{bucket}, and \code{nr_dummy_variables}.
#' @examples
#' \dontrun{
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
#' while (sum(is.na(data_matrix)) == 0)
#'   data_matrix[as.logical(round(runif(ncol(data_matrix) * nrow(data_matrix), -0.3, 0.7)))] <- NA
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' dataframe <- as.data.frame(data_matrix)
#' autovar(dataframe, selected_column_names = c('rumination', 'happiness'),
#'                    significance_levels = c(0.05, 0.01, 0.005),
#'                    test_names = c('portmanteau',
#'                                   'portmanteau_squared',
#'                                   'skewness'),
#'                    criterion = 'AIC',
#'                    imputation_iterations = 30,
#'                    measurements_per_day = 1)
#' }
#' @export
autovar <- function(raw_dataframe, selected_column_names, significance_levels = c(0.05, 0.01, 0.005),
                    test_names = c('portmanteau', 'portmanteau_squared', 'skewness'),
                    criterion = 'AIC', imputation_iterations = 30, measurements_per_day = 1) {
  data_matrix <- validate_raw_dataframe(raw_dataframe)
  params <- validate_params(data_matrix, list(selected_column_names = selected_column_names,
                                              significance_levels = significance_levels,
                                              test_names = test_names,
                                              criterion = criterion,
                                              imputation_iterations = imputation_iterations,
                                              measurements_per_day = measurements_per_day))
  data_matrix <- impute_datamatrix(data_matrix,
                                   params$measurements_per_day,
                                   params$imputation_iterations)
  number_of_measurements <- nrow(data_matrix)
  ln_data_matrix <- apply_ln_transformation(data_matrix)
  daypart_dummy_data <- daypart_dummies(number_of_measurements,
                                        params$measurements_per_day)
  day_dummy_data <- day_dummies(number_of_measurements,
                                params$measurements_per_day)
  trend_column_matrix <- trend_columns(number_of_measurements)
  number_of_endo_vars <- length(params$selected_column_names)
  cluster <- makeCluster(detectCores(),
                         type = "PSOCK",
                         useXDR = FALSE,
                         methods = FALSE)
  all_outlier_masks <- 0:(2^(number_of_endo_vars) - 1)
  significance_buckets <- c(params$significance_levels, 0)
  returned_models <- list()
  for (use_logtransform in c(FALSE, TRUE)) {
    best_models <- list()
    if (use_logtransform)
      endo_matrix <- ln_data_matrix[, params$selected_column_names]
    else
      endo_matrix <- data_matrix[, params$selected_column_names]
    for (lag in 1:2) {
      if (needs_trend(endo_matrix, lag))
        trend_dummies <- cbind(daypart_dummy_data, trend_column_matrix)
      else
        trend_dummies <- daypart_dummy_data
      best_model_in_lag <- list(model_score = Inf, bucket = 0, nr_dummy_variables = Inf)
      for (use_daydummies in c(FALSE, TRUE)) {
        if (use_daydummies) {
          if (is.null(day_dummy_data)) next
          exo_matrix <- cbind(trend_dummies, day_dummy_data)
        } else {
          exo_matrix <- trend_dummies
        }
        outlier_dummies <- residual_outliers(residuals(run_var(endo_matrix,
                                                               exo_matrix,
                                                               lag)),
                                             number_of_measurements)
        outlier_masks <- select_valid_masks(all_outlier_masks,
                                            invalid_mask(outlier_dummies))
        model_vector <- clusterMap(cluster, evaluate_model, outlier_masks,
                                   MoreArgs = list(endo_matrix = endo_matrix,
                                                   exo_matrix = exo_matrix,
                                                   lag = lag,
                                                   outlier_dummies = outlier_dummies,
                                                   test_names = params$test_names,
                                                   criterion = params$criterion,
                                                   logtransformed = use_logtransform,
                                                   significance_buckets = significance_buckets),
                                   SIMPLIFY = FALSE, USE.NAMES = FALSE)
        best_model <- list(model_score = Inf, bucket = 0, nr_dummy_variables = Inf)
        for (model in model_vector) {
          if (is.null(model)) next
          # Compare with models with more or fewer outlier dummies. Prefer models with fewer outlier dummies.
          best_model <- compete(best_model, model, TRUE)
        }
        # Compare with models with/without day dummies, here we just look at what gives the better model.
        # Here we are not necessarily interested in which model needs more or fewer outliers, since that
        # is affected by the inclusion of day dummies.
        best_model_in_lag <- compete(best_model_in_lag, best_model, FALSE)
      }
      # Compare with models of other lags, here we are interested in which lag needs fewer dummies.
      # The assumption here is that if weekly cyclicity is present, it would be present on models
      # of all lags, and therefore what we are comparing here is merely which lag needs fewer outlier
      # dummies to find a model in the same significance bucket.
      if (!is.infinite(best_model_in_lag$model_score))
        best_models <- insert_model_into_list(best_model_in_lag, best_models, TRUE)
    }
    # Merge models with/without logtransform. Since logtransform affects outliers, don't compare those.
    # The reasoning follows from above when including day dummies.
    returned_models <- merge_model_lists(returned_models, best_models, FALSE)
  }
  stopCluster(cluster)
  returned_models
}

evaluate_model <- function(outlier_mask, endo_matrix, exo_matrix, lag, outlier_dummies,
                           test_names, criterion, logtransformed, significance_buckets) {
  if (outlier_mask != 0) {
    selected_column_indices <- selected_columns(outlier_mask)
    exploded_outlier_dummies <- explode_dummies(as.matrix(outlier_dummies[, selected_column_indices]))
    exo_matrix <- cbind(exo_matrix, exploded_outlier_dummies)
  }
  varest <- run_var(endo_matrix, exo_matrix, lag)
  if (!model_is_stable(varest))
    return(NULL)
  significance_p_values <- run_tests(varest, test_names)
  model_significance <- min(significance_p_values)
  score <- model_score(varest, criterion, logtransformed)
  significance_bucket <- 0
  for (bucket in significance_buckets) {
    if (model_significance < bucket) next
    significance_bucket <- bucket
    break
  }
  list(logtransformed = logtransformed,
       lag = lag,
       varest = varest,
       model_score = score,
       bucket = significance_bucket,
       nr_dummy_variables = nr_dummy_variables(varest))
}

nr_dummy_variables <- function(varest) {
  outlier_dummies <- length(grep("^outlier_[0-9]+$", colnames(varest$datamat)))
  day_dummies <- length(grep("^day_[0-9]+$", colnames(varest$datamat)))
  if (day_dummies > 0)
    day_dummies <- 1
  outlier_dummies + day_dummies
}

insert_model_into_list <- function(model, model_list, compare_outliers) {
  if (length(model_list) == 0) return(list(model))
  if (challenger_wins(model, model_list[[length(model_list)]], compare_outliers))
    return(append(model_list, list(model)))
  left <- 1
  right <- length(model_list)
  while (left != right) {
    middle <- (left + right)%/%2
    if (challenger_wins(model_list[[middle]], model, compare_outliers))
      right <- middle
    else
      left <- middle + 1
  }
  append(model_list, list(model), after = left - 1)
}

merge_model_lists <- function(list_a, list_b, compare_outliers) {
  pos_a <- 1
  pos_b <- 1
  result <- list()
  len_a <- length(list_a)
  len_b <- length(list_b)
  while (pos_a <= len_a && pos_b <= len_b) {
    if (challenger_wins(list_a[[pos_a]], list_b[[pos_b]], compare_outliers)) {
      result <- append(result, list(list_b[[pos_b]]))
      pos_b <- pos_b + 1
    } else {
      result <- append(result, list(list_a[[pos_a]]))
      pos_a <- pos_a + 1
    }
  }
  if (pos_a <= len_a)
    result <- append(result, list_a[pos_a:len_a])
  if (pos_b <= len_b)
    result <- append(result, list_b[pos_b:len_b])
  result
}
