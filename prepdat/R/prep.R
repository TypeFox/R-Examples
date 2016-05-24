#' Creates One Finalized Table Ready for Statistical Analysis
#'
#'@description \code{prep()} aggregates a single dataset in a long format
#'   according to any number of grouping variables. This makes \code{prep()}
#'   suitable for aggregating data from various types of experimental designs
#'   such as between-subjects, within-subjects (i.e., repeated measures), and
#'   mixed designs (i.e., experimental designs that include both between- and
#'   within- subjects independent variables). \code{prep()} returns a data
#'   frame with a number of dependent measures for further analysis for each
#'   aggregated cell (i.e., experimental cell) according to the  provided
#'   grouping variables (i.e., independent variables). Dependent measures for
#'   each experimental cell include among others means before and after
#'   rejecting all observations according to a flexible standard deviation
#'   criteria, number of rejected observations according to the flexible
#'   standard deviation criteria, proportions of rejected observations
#'   according to the flexible standard deviation criteria, number of
#'   observations before rejection, means after rejecting observations
#'   according to procedures described in Van Selst & Jolicoeur (1994;
#'   suitable when measuring reaction-times), standard deviations, medians,
#'   means according to any percentile (e.g., 0.05, 0.25, 0.75, 0.95) and
#'   harmonic means. The data frame \code{prep()} returns can also be exported
#'   as a txt file to be used for statistical analysis in other statistical
#'   programs.
#'
#' @usage prep(
#'    dataset = NULL
#'    , file_name = NULL
#'    , id = NULL
#'    , within_vars = c()
#'    , between_vars = c()
#'    , dvc = NULL
#'    , dvd = NULL
#'    , keep_trials = NULL
#'    , drop_vars = c()
#'    , keep_trials_dvc = NULL
#'    , keep_trials_dvd = NULL
#'    , id_properties = c()
#'    , sd_criterion = c(1, 1.5, 2)
#'    , percentiles = c(0.05, 0.25, 0.75, 0.95)
#'    , outlier_removal = NULL
#'    , keep_trials_outlier = NULL
#'    , decimal_places = 4
#'    , notification = TRUE
#'    , dm = c()
#'    , save_results = TRUE
#'    , results_name = "results.txt"
#'    , save_summary = TRUE
#' )
#' @param dataset Name of the data frame in R that contains the long format
#'   table after merging the individual data files using
#'   \code{file_merge()}. Either \code{dataset} or \code{file_name} must be
#'   provided. Default is \code{NULL}.
#' @param file_name A string with the name of a txt or csv file (including the
#'   file extension, e.g. \code{"my_data.txt"}) with the merged
#'   dataset in case the user already merged the individual data files. Either
#'   \code{dataset} or \code{file_name} must be provided. Default is
#'   \code{NULL}.
#' @param id A string with the name of the column in \code{file_name} or in
#'   \code{dataset} that contains the variable specifying the case identifier
#'   (i.e., the variable upon which the measurement took place; e.g.,
#'   \code{"subject_number"}). This should be a unique value per case. Values
#'   in this column must be numeric. Argument must be provided. Default is
#'   \code{NULL}.
#' @param within_vars A vector with names of the grouping variables in
#'   \code{file_name} or in \code{dataset} that contain independent variables
#'   manipulated (or observed) within-ids (i.e., within-subjects, repeated
#'   measures). Single or multiple values must be specified as a string (e.g.,
#'   \code{c("SOA", "condition")}) according to the hierarchical order you
#'   wish. Note that the order of the names in \code{within_vars()} is
#'   important because \code{prep()} aggregates the data for the dependent
#'   measures by first dividing them to the levels of the first grouping
#'   variable in \code{witin_vars()}, and then within each of those levels
#'   \code{prep()} divides the data according to the next variable in
#'   \code{within_vars()} and so forth. Values in these columns must be
#'   numeric. Either \code{within_vars} or \code{between_vars} or both
#'   arguments must be provided. Default is \code{c()}.
#' @param between_vars A vector with names of the grouping variables in
#'   \code{file_name} or in \code{dataset} that contain independent variables
#'   manipulated (or observed) between-ids (i.e., between-subjects). Single
#'   or multiple values must be specified as a string (e.g., \code{c("order")}).
#'   Order of the names in \code{between_vars()} does not matter. Values in
#'   this column must be numeric. Either \code{between_vars} or
#'   \code{within_vars} or both arguments must be provided. Default is
#'   \code{c()}.
#' @param dvc A string with the name of the column in \code{file_name} or in
#'   \code{dataset} that contains the dependent variable (e.g., "rt" for
#'   reaction-time as a dependent variable). Values in this column must be in
#'   an interval or ratio scale. Either \code{dvc} or \code{dvd} or both
#'   arguments must be provided. Default is \code{NULL}.
#' @param dvd A string with the name of the column in \code{file_name} or in
#'   \code{dataset} that contains the dependent variable (e.g., \code{"ac"}
#'   for accuracy as a dependent variable). Values in this column must be
#'   numeric and discrete (e.g., 0 and 1). Either \code{dvc} or \code{dvd} or
#'   both arguments must be provided. Default is \code{NULL}.
#' @param keep_trials A string. Allows deleting unnecessary observations and keeping
#'   necessary observations in \code{file_name} or in \code{dataset} according
#'   to logical conditions specified as a string. For example, if the dataset
#'   contains practice trials for each subject, these trials should not be
#'   included in the aggregation. The user should remove these trials by
#'   specifying how they were coded in the raw data (i.e., data before
#'   aggregation). For example, if practice trials are all the ones for which
#'   the "block" column in the raw data tables equals to zero, the
#'   \code{keep_trials} argument should be \code{"raw_data$block !== 0"}.
#'   \code{raw_data} is the internal object in \code{prep()} representing the
#'   merged table. All logical conditions in \code{keep_trials} should be put
#'   in the same string and be concatenated by \code{&} or \code{|}. Logical
#'   conditions for this argument can relate to different columns in the merged
#'   table. Note that all further arguments of \code{prep()} will relate to the
#'   remaining observations in the merged table. Default is \code{NULL}.
#' @param drop_vars A vector with names of columns to delete in \code{file_name}
#'   or in \code{dataset}. Single or multiple values must be specified as a
#'   string (e.g., \code{c("font_size")}). Order of the names in
#'   \code{drop_vars} does not matter. Note that all further arguments of
#'   \code{prep()} will relate to the remaining variables in the merged table.
#'   Default is \code{c()}.
#' @param keep_trials_dvc A string. Allows deleting unnecessary observations
#'   and keeping necessary observations in \code{file_name} or in \code{dataset}
#'   for calculations and aggregation of the dependent variable in \code{dvc}
#'   according to logical conditions specified as a string. Logical conditions
#'   should be specified as a string as in the \code{keep_trials} argument
#'   (e.g., \code{"raw_data$rt > 100 & raw_data$rt < 3000 & raw_dada$ac == 1"}).
#'   All dependent measures for \code{dvc} except for those specified in
#'   \code{outlier_removal} will be calculated on the remaining observations.
#'   Defalut is \code{NULL}.
#' @param keep_trials_dvd A string. Allows deleting unnecessary observations
#'   and keeping necessary observations in \code{file_name} or in \code{dataset}
#'   for calculations and aggregation of the dependent variable in \code{dvd}
#'   according to logical conditions specified as a string. Logical conditions
#'   should be specified as a string as in the \code{keep_trials argument}
#'   (e.g., \code{raw_data$rt > 100 & raw_data$rt < 3000}). All dependent
#'   measures for \code{dvd} (i.e., \code{"mdvd"} and \code{"merr"}) will be
#'   calculated on the remaining observations. Default is \code{NULL}.
#' @param id_properties A vector with names of columns in \code{dataset} or in
#'   \code{file_name} that describe the ids (e.g., subjects) in the data and
#'   are not manipulated (or observed) within-or between-ids. For example, in
#'   case the user logged for each observation and for each id in an experiment
#'   also the age and the gender of the subject, this argument will be
#'   \code{c("age", "gender")}. Order of the names in \code{id_properties} does
#'   not matter. Single or multiple values must be specified as a string.
#'   Values in these columns must be numeric. Default is \code{c()}.
#' @param sd_criterion A vector specifying a number of standard deviation
#'   criteria for which \code{prep()} will calculate the mean \code{dvc} for
#'   each cell in the finalized table after rejecting observations that did not
#'   meet the criterion (e.g., rejecting all observations that were more than 2
#'   standard deviations above or below the mean of that cell). Values in this
#'   vector must be numeric. Default is \code{c(1, 1.5, 2)}.
#' @param percentiles A vector containing wanted percentiles for \code{dvc}.
#'   Values in this vector must be decimal numbers between 0 to 1. Percentiles
#'   are calculated according to \code{type = 7} (see
#'   \code{\link[stats]{quantile}} for more information). Default is
#'   \code{c(0.05, 0.25, 0.75, 0.95)}.
#' @param outlier_removal Numeric. Specifies which outlier removal procedure
#'   with moving criterion to calculate for \code{dvc} according to procedures
#'   described by Van Selst & Jolicoeur (1994). If \code{1} then non-recursive
#'   procedure is calculated, if \code{2} then modified recursive procedure is
#'   calculated, if \code{3} then hybrid recursive procedure is calculated.
#'   Moving criterion is according to Table 4 in Van Selst & Jolicoeur (1994).
#'   If experimental cell has 4 trials or less it will result in \code{NA}.
#'   Default is \code{NULL}.
#' @param keep_trials_outlier A string. Allows deleting unnecessary
#'   observations and keeping necessary observations in \code{file_name} or in
#'   \code{dataset} for calculations and aggregation of the outlier removal
#'   procedures by Van Selst & Jolicoeur (1994). Logical conditions should be
#'   specified as a string as in the \code{keep_trials} argument (e.g.,
#'   \code{"raw_data$ac == 1"}). \code{outlier_removal} procedure will be
#'   calculated on the remaining observations. Defalut is \code{NULL}.
#' @param decimal_places Numeric. Specifies number of decimals to be written
#'   in \code{results_name} for each value of the dependent measures for
#'   \code{dvc}. Value must be numeric. Default is \code{4}.
#' @param notification Logical. If TRUE, prints messages about the progress of
#'   the function. Default is \code{TRUE}.
#' @param dm a Vector with names of dependent measures the function returns. If
#'   empty (i.e., \code{c()}) the function returns a data frame with all
#'   possible dependent measures in \code{prep()}. Values in this vector must
#'   be strings from the following list: "mdvc", "sdvc", "meddvc", "tdvc",
#'   "ntr", "ndvc", "ptr", "prt", "rminv", "mdvd", "merr". Default is
#'   \code{c()}. See return for more details.
#' @param save_results Logical. If TRUE, the function creates a txt file
#'   containing the returned data frame. Default is \code{TRUE}.
#' @param results_name A string of the name of the data frame the function
#'   returns in case \code{save_results} is TRUE. Default is
#'   \code{"results.txt"}.
#' @param save_summary Logical. if TRUE, creates a summary txt file. Default is
#'   \code{TRUE}.
#' @references Grange, J.A. (2015). trimr: An implementation of common response
#'  time trimming methods. R Package Version 1.0.1.
#'  \url{https://CRAN.R-project.org/package=trimr}
#'
#' Van Selst, M., & Jolicoeur, P. (1994). A solution to the effect of sample
#' size on outlier elimination. \emph{The quarterly journal of experimental
#' psychology, 47}(3), 631-650.
#'
#'
#' @return A data frame with dependent measures for the dependent variables in
#'   \code{dvc} and \code{dvd} by \code{id} and grouping variables:
#'
#'      \code{mdvc}: mean \code{dvc}.
#'
#'      \code{sdvc}: SD for \code{dvc}.
#'
#'      \code{meddvc}: median \code{dvc}.
#'
#'      \code{tdvc}: mean \code{dvc} after rejecting observations above
#'      standard deviation criteria specified in \code{sd_criterion}.
#'
#'      \code{ntr}: number of observations rejected for each standard deviation
#'      criterion specified in \code{sd_criterion}.
#'
#'      \code{ndvc}: number of observations before rejection.
#'
#'      \code{ptr}: proportion of observations rejected for each standard
#'      deviation criterion specified in \code{sd_criterion}.
#'
#'      \code{rminv}: harmonic mean of \code{dvc}.
#'
#'      \code{prt}: \code{dvc} according to each of the percentiles specified
#'      in \code{percentiles}.
#'
#'      \code{mdvd}: mean \code{dvd}.
#'
#'      \code{merr}: mean error.
#'
#'      \code{nrmc}: mean \code{dvc} according to non-recursive procedure with
#'      moving criterion.
#'
#'      \code{nnrmc}: number of observations rejected for \code{dvc} according
#'      to non-recursive procedure with moving criterion.
#'
#'      \code{pnrmc}: percent of observations rejected for \code{dvc} according
#'      to non-recursive procedure with moving criterion.
#'
#'      \code{tnrmc}: total number of observations upon which the non-recursive
#'      procedure with moving criterion was applied.
#'
#'      \code{mrmc}: mean \code{dvc} according to modified-recursive procedure
#'      with moving criterion.
#'
#'      \code{nmrmc}: number of observations rejected for \code{dvc} according
#'      to modified-recursive procedure with moving criterion.
#'
#'      \code{pmrmc}: percent of observations rejected for \code{dvc} according
#'      to modified-recursive procedure with moving criterion.
#'
#'      \code{tmrmc}: total number of observations upon which the
#'      modified-recursive procedure with moving criterion was applied.
#'
#'      \code{hrmc}: mean \code{dvc} according to hybrid-recursive procedure
#'      with moving criterion.
#'
#'      \code{nhrmc}: number of observations rejected for \code{dvc} according
#'      to hybrid-recursive procedure with moving criterion.
#'
#'      \code{thrmc}: total number of observations upon which the
#'      hybrid-recursive procedure with moving criterion was applied.
#' @export
#' @examples
#' data(stroopdata)
#' finalized_data <- prep(
#'          dataset = stroopdata
#'          , file_name = NULL
#'          , id = "subject"
#'          , within_vars = c("block", "target_type")
#'          , between_vars = c("order")
#'          , dvc = "rt"
#'          , dvd = "ac"
#'          , keep_trials = NULL
#'          , drop_vars = c()
#'          , keep_trials_dvc = "raw_data$rt > 100 & raw_data$rt < 3000 & raw_data$ac == 1"
#'          , keep_trials_dvd = "raw_data$rt > 100 & raw_data$rt < 3000"
#'          , id_properties = c()
#'          , sd_criterion = c(1, 1.5, 2)
#'          , percentiles = c(0.05, 0.25, 0.75, 0.95)
#'          , outlier_removal = 2
#'          , keep_trials_outlier = "raw_data$ac == 1"
#'          , decimal_places = 0
#'          , notification = TRUE
#'          , dm = c()
#'          , save_results = TRUE
#'          , results_name = "results.txt"
#'          , save_summary = TRUE
#'       )
#'
prep <- function(dataset = NULL, file_name = NULL, id = NULL,
                 within_vars =  c(), between_vars = c(), dvc = NULL,
                 dvd = NULL, keep_trials = NULL, drop_vars = c(),
                 keep_trials_dvc = NULL, keep_trials_dvd = NULL,
                 id_properties = c(), sd_criterion = c(1, 1.5, 2),
                 percentiles = c(0.05, 0.25, 0.75, 0.95),
                 outlier_removal = NULL, keep_trials_outlier = NULL,
                 decimal_places = 4, notification = TRUE, dm = c(),
                 save_results = TRUE, results_name = "results.txt",
                 save_summary = TRUE) {
  ## Error handling
  # Checks if large dataset was provided
  if (is.null(dataset) & is.null(file_name)) {
    stop("Oops! Did not find dataset or file_name. Either dataset or file_name must be provided")
  }

  # Get results_name extension
  extension <- substr(results_name, nchar(results_name) - 3,
                      nchar(results_name))
  # Check if extension is correct
  if (extension != ".txt") {
    stop("Oops! results_name must end with txt extension")
  }

  # Check if id was provided
  if (is.null(id)) {
    stop("Oops! id is missing. Please provide name of id column")
  }

  # Check if independent variables were provided
  if (length(within_vars) == 0 & length(between_vars) == 0) {
    stop("Oops! Did not find independent variables. You must enter an independent varible to either within_vars,
          between_vars (or to both)")
  }

  # Check if dependent variables were provided
  if (is.null(dvc) & is.null(dvd)) {
    stop("Oops! Did not find dvc or dvd. You must enter at least one dependent variable")
  }

  ## Number of decimal places for dvd
  decimal_dvd <- 3

  ## Create name for summary file
  sum_file_name <- paste(substr(results_name, 0, nchar(results_name) - 4),
                         "_summary.txt", sep = "")

  ## Read table
  if (!is.null(file_name)) {
    # Call read_data()
    raw_data <- read_data(file_name, notification)
    # For summary txt file (to be used later in the code)
    dataname <- file_name
  } else {
    # Take table from dataset
    raw_data <- dataset
    # For summary txt file (to be used laster in the code)
    dataname <- "dataset"
  }

  ## For summary txt file (to be used later in the code)
  if (!is.null(outlier_removal)) {
    outlier_name <- c("non-recursive", "modified-recursive", "hybrid-recursive")
    outlier_name <- paste(outlier_name[outlier_removal], "procedure with moving criterion", sep = " ")
  }

  ## Print dimensions of the large dataset to console
  if (notification == TRUE) {
    # Message
    message(paste("raw_data has", dim(raw_data)[1], "observations and", dim(raw_data)[2], "variables"))
    print(head(raw_data))
  }

  ## Save dimensions of the large dataset before doing anything
  dim_raw_data1 <- dim(raw_data)

  ## Delete unnecessary trials in case keep_trials is not NULL
  # All further calculations will be done on the remaining trials
  if (!is.null(keep_trials)) {
    if (notification == TRUE) {
      # Message
      message("Keeping trials according to keep_trials and deleting unnecessary trials in raw_data")
      message("All further calculations will be done on the remaining trials")
    }
    # Save keep_trials as a string before parsing in order to later use when
    # writing Summary file
    keep_trials_sum <- keep_trials
    # Subset trials
    keep_trials <- eval(parse(text = keep_trials))
    raw_data <- raw_data[keep_trials, ]
    if (notification == TRUE) {
      # Message
      message(paste("raw_data has", dim(raw_data)[1], "observations and", dim(raw_data)[2], "variables"))
      print(head(raw_data))
    }
  } else {
    keep_trials_sum <- class(keep_trials)
  }
  # End of !is.null(keep_trials)

  ## Delete unnecessary variables in case drop_vars is not NULL
  # All further calculations will be done on the remaining variables
  if (length(drop_vars) > 0) {
    if (notification == TRUE) {
      # Message
      message("Dropping variables in raw_data according to drop_vars")
      message("All further calculations will be done on the remaining variables")
    }
    # Save drop_vars as a string before parsing in order to later use when
    # writing Summary file
    drop_col_sum <- drop_vars

    # Subset variables by keeping all columns except the ones in drop_vars
    raw_data <- raw_data[, !(colnames(raw_data) %in% drop_vars)]
    if (notification == TRUE) {
      # Message
      message(paste("raw_data has", dim(raw_data)[1], "observations and", dim(raw_data)[2], "variables"))
      print(head(raw_data))
    }
  } else {
    drop_col_sum <- "NULL"
  }
  # End of if(length(drop_vars) > 0)

  ## Save dimensions of raw_data after removing unnecessary trials and
  # variables according to keep_trials and drop_vars
  dim_raw_data2 <- dim(raw_data)

  ## Creates an array in the size of number of ids in raw_data
  id_col <- tapply(raw_data[[id]], list(raw_data[[id]]), mean)

  ## Calculates independent between-subjects variables in case they exists
  if (length(between_vars) > 0) {
    # Found between-subjects variables
    if(notification == TRUE) {
      # Message
      if (length(between_vars) == 1) {
        message(paste("Found", length(between_vars), "between-subjects independent variable:", between_vars))
      } else {
        message(paste("Found", length(between_vars), "between-subjects independent variables:", between_vars))
      }
    }
    # Creates a data frame for between-subject variables according to the
    # number of ids
    between_vars_df <- data.frame(id_col)
    if (notification == TRUE) {
      # Message
      message("Calculating between-subjects independent variables")
    }
    # Calculates between-subjects variables
    # Reset i to 1
    i <- 1
    while (i <= length(between_vars)) {
      between_vars_df[i] <- tapply(raw_data[[between_vars[i]]],
                                   list(raw_data[[id]]), mean)
      i <- i + 1
    }
    if (notification == TRUE) {
      # Message
      message("Giving names for between-subjects variables")
    }
    # Give names for between-subjects variables
    colnames(between_vars_df) <- between_vars
    if (notification == TRUE) {
      # Message
      message("Finished between-subjects independent variables")
    }
  }

  ## Calculates id_properties in case they exists
  if (length(id_properties) > 0) {
    # Found id properties
    if (notification == TRUE) {
      # Message
      message(paste("Found", length(id_properties), "id properties"))
    }
    # Creates a data frame for id_properties according to the number of subjects
    id_properties_df <- data.frame(id_col)
    if (notification == TRUE) {
      # Message
      message("Calculating id_properties")
    }
    # Calculates id properties
    # Reset i to 1
    i <- 1
    while (i <= length(id_properties)) {
      id_properties_df[i] <- tapply(raw_data[[id_properties[i]]],
                                    list(raw_data[[id]]), mean)
      i <- i + 1
    }
    if (notification == TRUE) {
      # Message
      message("Giving names for id_properties variables")
    }
    # Give names for id properties variables
    colnames(id_properties_df) <- id_properties
    if (notification == TRUE) {
      # Message
      message("Finished id_properties")
    }
  }

  ## Creates within_vars column in case within-subject indpendent variables
  # exists
  if (length(within_vars) > 0) {
    # Creates within_col: a data frame for within_vars variables
    if (notification == TRUE) {
      # Message
      message(paste("Found", length(within_vars), "within-subject independent variables"))
    }
    # within_col is in the same length as number of trials in raw_data
    within_col <- data.frame(1:length(raw_data[, 1]))
    # Reset l to 1
    l <- 1
    # Get all within_vars variables in raw_data into within_col
    while(l <= length(within_vars)) {
      within_col[l] <- raw_data[[within_vars[l]]]
      l <- l + 1
    }
    # Give names to within_col
    colnames(within_col) <- within_vars
    # Creates a column for within-subject independent variables
    within_condition <- do.call(interaction, c(within_col, sep = "_",
                                               lex.order = TRUE))
    # Append within_condition to raw_data
    raw_data$within_condition <- factor(within_condition)
    if (notification == TRUE) {
      # Message
      message("within_condition column was added to raw_data")
      print(head(raw_data))
      message("Finished within_condition column")
    }
  }

  ## Creates a data frame for dvc data and delete unnecessary trials for dvc in
  # case dvc exists
  if (!is.null(dvc)) {
    # Found dvc data
    if (notification == TRUE) {
      # Message
      message("Found dvc data")
      message("Creats a data frame for dvc data")
    }
    # Delete unnecessary trials in case keep_trials_dvc is not NULL
    # All further calculations on dvc except outlier removal procedures will
    # be done on the remaining trials
    if (!is.null(keep_trials_dvc)) {
      if (notification == TRUE) {
        # Message
        message("Keeping trials according to keep_trials_dvc and deleting unnecessary trials in raw_data_dvc")
        message("All further calculations on dvc except outlier removal procedures will be done on the remaining trials")
      }
      # Save keep_trials_dvc as a string before parsing in order to later use
      # when writing Summary file
      keep_trials_dvc_sum <- keep_trials_dvc
      # Subsetting trials
      keep_trials_dvc <- eval(parse(text = keep_trials_dvc))
      raw_data_dvc <- raw_data[keep_trials_dvc, ]
    } else {
      # Save keep_trials_dvc as "NULL" in order to later use when writing
      # Summary file
      keep_trials_dvc_sum <- class(keep_trials_dvc)
      # Save raw_data as raw_data_dvc
      raw_data_dvc <- raw_data
    }
    if (notification == TRUE) {
      # Message
      message(paste("raw_data_dvc has", dim(raw_data_dvc)[1], "observations and", dim(raw_data_dvc)[2], "variables"))
      print(head(raw_data_dvc))
    }
    # Get dimensions of raw_data_dvc for summary file
    dim_raw_data_dvc <- dim(raw_data_dvc)
  }

  ## Creates a data frame for dvd data and delete unnecessary trials for dvd
  # in case dvd exists
  if (!is.null(dvd)) {
    # Found dvd data
    if (notification == TRUE) {
      # Message
      message("Found dvd data")
      message("Creats a data frame for dvd data")
    }
    # Delete unnecessary trials in case keep_trials_dvd is not NULL
    # All further calculations on dvd will be done on the remaining trials
    if (!is.null(keep_trials_dvd)) {
      # Delete unnecessary trials
      if (notification == TRUE) {
        # Message
        message("Keeping trials according to keep_trials_dvd and deleting unnecessary trials in raw_data_dvd")
        message("All further calculations on dvd will be done on the remaining trials")
      }
      # Save keep_trials_dvd as a string before parsing in order to later use
      # when writing Summary file
      keep_trials_dvd_sum <- keep_trials_dvd
      # Subsetting trials
      keep_trials_dvd <- eval(parse(text = keep_trials_dvd))
      # Save raw_data after subsetting
      raw_data_dvd <- raw_data[keep_trials_dvd, ]
    } else {
      # Save keep_trials_dvd as "NULL" in order to later use when writing
      # Summary file
      keep_trials_dvd_sum <- class(keep_trials_dvd)
      # Save raw_data as raw_data_dvd
      raw_data_dvd <- raw_data
    }
    if (notification == TRUE) {
      # Message
      message(paste("raw_data_dvd has", dim(raw_data_dvd)[1], "observations and", dim(raw_data_dvd)[2], "variables"))
      print(head(raw_data_dvd))
    }
    # Get dimensions of raw_data_dvd for summary file
    dim_raw_data_dvd <- dim(raw_data_dvd)
  }

  # Creates a data frame for dependent measures
  dm_df <- data.frame(1:length(id_col))

  # Reset counter for dependent measures to 1
  count <- 1

  ## Calculates dependent measures
  # Checks if the design is a between-subjects design, within-subjects design
  # or a mixed design
  if (length(between_vars) > 0 & length(within_vars) == 0) {
    # Between-subjects design: all dependent measures will be calculated by id
    if (notification == TRUE) {
      message("Your design is a between-subjects design")
    }

    # dvc
    # Calculates dependent measures for dvc in case dvc exists
    if (!is.null(dvc)) {

      # mdvc: mean dvc
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "mdvc" among others (i.e., "mdvc" %in% dm equals to TRUE)
      if (length(dm) == 0 | "mdvc" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating mean dvc")
        }
        dm_df[count] <- tapply(raw_data_dvc[[dvc]], list(raw_data_dvc[[id]]),
                               mean)
        colnames(dm_df)[count] <- "mdvc"
        # Round mdvc according to decimal_places
        dm_df[count] <- round(dm_df[count], digits = decimal_places)
        count <- count + 1
      }

      # sdvc: standard deviation (SD) for dvc using denominator n
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "sdvc" among others (i.e., "sdvc" %in% dm equals to TRUE)
      if (length(dm) == 0 | "sdvc" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating SD for dvc using denominator n")
        }
        # SD using denominator n - 1
        dvc_sd_r <- tapply(raw_data_dvc[[dvc]], list(raw_data_dvc[[id]]), sd)
        # Length of each cell in the raw_data_dvc
        l_dvc <- tapply(raw_data_dvc[[dvc]], list(raw_data_dvc[[id]]), length)
        # Calculating sdvc using denominator n
        dm_df[count] <- dvc_sd_r * sqrt((l_dvc - 1) / (l_dvc))
        colnames(dm_df)[count] <- "sdvc"
        # Round sdvc according to decimal_places
        dm_df[count] <- round(dm_df[count], digits = decimal_places)
        count <- count + 1
      }

      # meddvc: median dvc
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "medrt" among others (i.e., "meddvc" %in% dm equals to TRUE)
      if (length(dm) == 0 | "meddvc" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating median dvc")
        }
        dm_df[count] <- tapply(raw_data_dvc[[dvc]], list(raw_data_dvc[[id]]),
                               median)
        colnames(dm_df)[count] <- "meddvc"
        # Round meddvc according to decimal_places
        dm_df[count] <- round(dm_df[count], digits = decimal_places)
        count <- count + 1
      }

      # Z scores (to be used later when removing values according to
      # sd_criterion)
      # deviation from the mean
      dev_from_mean  <- ave(raw_data_dvc[[dvc]], raw_data_dvc[[id]],
                            FUN = function(x) x - mean(x))
      raw_data_dvc$dev_from_mean <- dev_from_mean
      # sdvc for z scores (again using denominator n)
      sdvc_for_zscores <- ave(raw_data_dvc[[dvc]], raw_data_dvc[[id]],
                              FUN = function(x) sd(x) * sqrt((length(x) - 1 ) / (length(x))))
      raw_data_dvc$sdvc_for_zscores <- sdvc_for_zscores
      # z scores
      z_score <- raw_data_dvc$dev_from_mean / raw_data_dvc$sdvc_for_zscores
      raw_data_dvc$z_score <- z_score

      # tdvc: trimmed dvc according to SD in sd_criterion
      # Create an empty vector for names of sd_criterion coulmns
      sd_criterion_names <- c()
      # Reset l to 1
      l <- 1
      # Fill empty vector with names of sd_criterion columns
      while (l <= length(sd_criterion)) {
        sd_criterion_names[l] <- paste("t", sd_criterion[l], "dvc" ,sep = "")
        l <- l + 1
      }
      # Change names of id and dvc columns because later I use "dplyr" which
      # makes it easier
      # Change name of id column to "id"
      names(raw_data_dvc)[names(raw_data_dvc) == id] <- "id"
      # Change name of dvc column to "dvc"
      names(raw_data_dvc)[names(raw_data_dvc) == dvc] <- "dvc"
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "tdvc" among others (i.e., "tdvc" %in% dm equals to TRUE)
      if (length(dm) == 0 | "tdvc" %in% dm) {
        # Reset j to 1
        j <- 1
        # Calculating mean dvc after rejecting values above SD criterions
        # according to sd_criterion
        while (j <= length(sd_criterion)) {
          if (notification == TRUE) {
            # Message
            message(paste("Calculating mean dvc after rejecting values above", sd_criterion[j], "SD",
                    sep = " "))
          }
          temp <- raw_data_dvc %>%
                  dplyr::group_by(id) %>%
                  dplyr::summarise(n = mean(dvc[abs(z_score) < sd_criterion[j]])) %>%
                  reshape2::dcast(id ~ sd_criterion_names[j], value.var = "n", fun = mean)
          dm_df[count] <- temp[-1]
          # Round according to decimal_places
          dm_df[count] <- round(dm_df[count], digits = decimal_places)
          count <- count + 1
          j <- j + 1
        }
      }

      # ntr: number of trimmed values for each SD in sd_criterion
      # Checks if the user does not want ntr in results (i.e., !("ntr" %in% dm)
      # equals to FALSE)
      if (!("ntr" %in% dm)) {
        # In case the user does not want ntr in the filnal table,
        # store ntr in temp_dm to be used later when calculating proportion
        temp_dm <- data.frame(1:length(id_col))
        count_temp <- 1
      }
      # Create an empty vector for names of ntr columns
      ntr_names <- c()
      # Rest l to 1
      l <- 1
      # Fill empty vector with names of ntr columns
      while (l <= length(sd_criterion)) {
        ntr_names[l] <- paste("n", sd_criterion[l], "tr" ,sep = "")
        l <- l + 1
      }
      # Reset j to 1 for next while loop
      j <- 1
      # Counting number of rejected values for each sd_criterion
      while (j <= length(sd_criterion)) {
        if (notification == TRUE) {
        # Messgae
        message(paste("Counting number of rejected values for", sd_criterion[j], "SD", sep = " "))
        }
        temp <- raw_data_dvc %>%
                dplyr::group_by(id) %>%
                dplyr::summarise(n = sum(abs(z_score) > sd_criterion[j])) %>%
                reshape2::dcast(id ~ ntr_names[j], value.var = "n", fun = sum)
        if (length(dm) == 0 | "ntr" %in% dm) {
          # Store it in the dependent measures df
          dm_df[count] <- temp[-1]
          count <- count + 1
        } else {
          # Stroe it in a temporary df
          temp_dm[count_temp] <- temp[-1]
          count_temp <- count_temp + 1
        }
        j <- j + 1
      }
      if (exists("temp_dm")) {
        temp_dm <- temp_dm[-1]
      }
      # End of ntr

      # ndvc: number of values before rejecting values according to sd_criterion
      # Calculate ndvc anyway
      ndvc <- tapply(raw_data_dvc$dvc, list(raw_data_dvc$id), length)
      # Checks if to store ndvc in dm_df ( in case the user wants all dpendent
      # measures (i.e., length(dm) == 0) or he wants "ndvc" among others
      # (i.e., "ndvc" %in% dm equals to TRUE)
      if (length(dm) == 0 | "ndvc" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Counting number of values before rejecting values according to sd_criterion")
        }
        dm_df[count] <- ndvc
        colnames(dm_df)[count] <- "ndvc"
        count <- count + 1
      }
      # End of ndvc

      # ptr: proportion of rejected values according to sd_criterion
      # Create an empty vector for names of ptr columns
      ptr_names <- c()
      # Rest l to 1
      l <- 1
      # Fill empty vector with names of ptr columns
      while (l <= length(sd_criterion)) {
        ptr_names[l] <- paste("p", sd_criterion[l], "tr" ,sep = "")
        l <- l + 1
      }
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "ptr" among others (i.e., "ptr" %in% dm equals to TRUE)
      if (length(dm) == 0 | "ptr" %in% dm) {
        # Reset j to 1 for next while loop
        j <- 1
        while (j <= length(sd_criterion)) {
          if (notification == TRUE) {
            # Message
            message(paste("Calculating proportion of values rejected for", sd_criterion[j], "SD", sep = " "))
          }
          # Get ntr from temp_dm or from dm_df
          if (!(length(dm) == 0 | "ntr" %in% dm)) {
            temp <- temp_dm[[ntr_names[j]]] / ndvc
          } else {
            temp <- dm_df[[ntr_names[j]]] / ndvc
          }
          dm_df[count] <- temp
          colnames(dm_df)[count] <- ptr_names[j]
          # Round according to digits = 3
          dm_df[count] <- round(dm_df[count], digits = 3)
          count <- count + 1
          j <- j + 1
        }
      }
      # End of ptr

      # rminv: harmonic mean
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "rminv" among others (i.e., "rminv" %in% dm equals to TRUE)
      if (length(dm) == 0 | "rminv" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating harmonic mean for dvc")
        }
        dm_df[count] <- tapply(raw_data_dvc$dvc, list(raw_data_dvc$id),
                               psych::harmonic.mean)
        colnames(dm_df)[count] <- "rminv"
        # Round rminv according to decimal_places
        dm_df[count] <- round(dm_df[count], digits = decimal_places)
        count <- count + 1
      }
      # End of rminv

      # Percentiles: dvc according to requested percentiles
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "pdvc" among others (i.e., "pdvc" %in% dm equals to TRUE)
      if (length(dm) == 0 | "pdvc" %in% dm) {
        # Create an empty vector for names of pdvc columns
        percentails_names <- c()
        # Rest l to 1
        l <- 1
        # Fill vector
        while (l <= length(percentiles)) {
          percentails_names[l] <- paste("p", percentiles[l], "dvc" ,sep = "")
          l <- l + 1
        }
        # Reset j to 1 before next while loop
        j <- 1
        # Calculating percentiles
        while (j <= length(percentiles)) {
          if (notification == TRUE) {
            # Message
            message(paste("Calculating the", percentiles[j], "percentail for dvc" ), sep = "")
          }
          dm_df[count] <- tapply(raw_data_dvc$dvc, list(raw_data_dvc$id),
                                 quantile, probs = percentiles[j])
          colnames(dm_df)[count] <- percentails_names[j]
          # Round according to decimal_places
          dm_df[count] <- round(dm_df[count], digits = decimal_places)
          count <- count + 1
          j <- j + 1
        }
      }
      # End of pdvc

      # Outlier removal procedures according to Van Selst & Jolicoeur (1994)
      # Checks if the user wants an outlier removal procedure
      if (!is.null(outlier_removal)) {
        # User wants an outlier removal procedure

        # Delete unnecessary trials in case keep_trials_outlier is not NULL
        # Outlier removal procedure will be calculated on the remaining trials
        if (!is.null(keep_trials_outlier)) {
          # Save keep_trials_outlier as a string before parsing in order to
          # later use when writing Summary file
          keep_trials_outlier_sum <- keep_trials_outlier
          # Delete unnecessary trials
          if (notification == TRUE) {
            # Message
            message("Keeping trials according to keep_trials_outlier and deleting unnecessary trials in raw_data")
            message("Outlier removal procedure will be calculated on the remaining trials")
          }
          # Subsetting trials
          keep_trials_outlier <- eval(parse(text = keep_trials_outlier))
          # Save raw_data after subsetting
          raw_data <- raw_data[keep_trials_outlier, ]
        } else {
          # Save keep_trials_outlier as "NULL" in order to later use when
          # writing Summary file
          keep_trials_outlier_sum <- class(keep_trials_outlier)
          # Save raw_data
          raw_data <- raw_data
        }

        if (notification == TRUE) {
          # Message
          message("Calculating mean dvc by id according to selected outlier removal procedure")
        }

        # Change name of id column to "id"
        names(raw_data)[names(raw_data) == id] <- "id"
        # Change name of dvc column to "dvc"
        names(raw_data)[names(raw_data) == dvc] <- "dvc"

        # Create final array that will contain trimmed means per id
        final_data <- matrix(0, nrow = length(id_col), ncol = 1)
        # Create final array that will contain number of trials trimmed per id
        numtrials_data <- matrix(0, nrow = length(id_col), ncol = 1)
        # Create final array that will contain percent of trials removed per id
        percent_data <- matrix(0, nrow = length(id_col), ncol = 1)
        # Create final array that will contain total number of trials per id
        totaltrials_data <- matrix(0, nrow = length(id_col), ncol = 1)

        # Give names to columns according to the procedure
        if (outlier_removal == 1) {
          # Non recursive moving criterion procedure
          # nrmc: means for dvc per id
          colnames(final_data) <- "nrmc"
          # nnrmc: number of trials trimmed per id
          colnames(numtrials_data) <- "nnrmc"
          # pnrmc: percent of trials removed per id
          colnames(percent_data) <- "pnrmc"
          # tnrmc: total number of trials per id
          colnames(totaltrials_data) <- "tnrmc"
        } else if (outlier_removal == 2) {
          # Modified recursive moving criterion procedure
          # mrmc: means for dvc per id
          colnames(final_data) <- "mrmc"
          # nmrmc: number of trials trimmed per id
          colnames(numtrials_data) <- "nmrmc"
          # pmrmc: percent of trials removed per id
          colnames(percent_data) <- "pmrmc"
          # tmrmc: total number of trials per id
          colnames(totaltrials_data) <- "tmrmc"
        } else if (outlier_removal == 3) {
          # Hybrid recursive moving criterion procedure
          # hrmc: means for dvc per id
          colnames(final_data) <- "hrmc"
          # nhrmc: number of trials trimmed per id
          colnames(numtrials_data) <- "nhrmc"
          # phrmc: percent of trials removed per id
          colnames(percent_data) <- "phrmc"
          # thrmc: total number of trials per id
          colnames(totaltrials_data) <- "thrmc"
        }

        # Do outlier removal procedure
        for (i in 1:length(id_col)) {
          # For each id (i.e., subject)
          # Isolate the current participant's data
          temp_data <- raw_data[raw_data$id == id_col[i], ]

          if (outlier_removal == 1) {
            # Non recursive procedure with moving criterion
            # Do current id trimming
            current_id <- non_recursive_mc(temp_data$dvc)
            final_data[i, 1] <- current_id[1]
            percent_data[i, 1] <- current_id[2]
            numtrials_data[i, 1] <- current_id[3]
            totaltrials_data[i, 1] <- current_id[4]
          } else if (outlier_removal == 2) {
            # Modified recursive procedure with moving criterion
            # Do current id trimming
            current_id <- modified_recursive_mc(temp_data$dvc)
            final_data[i, 1] <- current_id[1]
            percent_data[i, 1] <- current_id[2]
            numtrials_data[i, 1] <- current_id[3]
            totaltrials_data[i, 1] <- current_id[4]
          } else if (outlier_removal == 3) {
            # Hybrid recursive procedure with moving criterion
            # Do current id trimming
            current_id <- hybrid_recursive_mc(temp_data$dvc)
            final_data[i, 1] <- current_id[1]
            percent_data[i, 1] <- current_id[2]
            totaltrials_data[i, 1] <- current_id[3]  # Should be 3. Do not change
          }
        }
        # End of id for loop

        # Round final_data according to decimal_places
        final_data <- round(final_data, digits = decimal_places)
        # Round percent_data according to digits = 3
        percent_data <- round(percent_data, digits = 3)

        if (notification == TRUE) {
          # Messgae
          if (outlier_removal == 1) {
            message("Finished calculating means for dvc by id according to non-recursive procedure with moving criterion ")
          } else if (outlier_removal == 2) {
            message("Finished calculating means for dvc by id according to modified-recursive procedure with moving criterion")
          } else if(outlier_removal == 3) {
            message("Finished calculating means for dvc by id according to hybrid procedure with moving criterion")
          }
        }
      }
      # End of outlier removal

      # Finished dvc
      if (notification == TRUE) {
        # Message
        message("Finished calculating dependent measures for dvc")
      }
    }  # End of !is.null(dvc)

    # dvd
    # Check if dvd exists
    if(!is.null(dvd)) {
      # Change name of id column to "id"
      names(raw_data_dvd)[names(raw_data_dvd) == id] <- "id"

      # mdvd: mean dvd
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "mdvd" among others (i.e., "mdvd" %in% dm equals to TRUE)
      if (length(dm) == 0 | "mdvd" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating mean dvd")
        }
        dm_df[count] <- tapply(raw_data_dvd[[dvd]], list(raw_data_dvd$id), mean)
        colnames(dm_df)[count] <- "mdvd"
        # Round mdvd according to decimal_dvd
        dm_df[count] <- round(dm_df[count], digits = decimal_dvd)
        count <- count + 1
      }

      # merr: mean error
      # Checks if the user wants all dpendent measures (i.e., length(dm) == 0)
      # or he wants "merr" among others (i.e., "merr" %in% dm equals to TRUE)
      if (length(dm) == 0 | "merr" %in% dm) {
        if (notification == TRUE) {
         # Message
         message("Calculating mean error")
        }
       dm_df[count] <- 1 - dm_df[["mdvd"]]
        colnames(dm_df)[count] <- "merr"
        # Round merr according to decimal_dvd
        dm_df[count] <- round(dm_df[count], digits = decimal_dvd)
        if (notification == TRUE) {
          # Message
          message("Finished calculating dependent measures for dvd")
        }
      }
    } # End of !is.null(dvd)

    # Creates results file
    if (notification == TRUE) {
      # Message
      message("Creating results")
    }

    # Check if to bind id_properties_df
    if (length(id_properties) > 0) {
      results <- cbind(id_col, id_properties_df, between_vars_df, dm_df)
    } else {
      results <- cbind(id_col, between_vars_df, dm_df)
    }

    # Check if to bind results of outlier removal procedure
    if (!is.null(outlier_removal)) {
      if (outlier_removal == 1 | outlier_removal == 2) {
        results <- cbind(results, final_data, percent_data, numtrials_data,
                         totaltrials_data)
      } else if (outlier_removal == 3) {
        results <- cbind(results, final_data, percent_data, totaltrials_data)
      }
    }

    # Change back name of id to original
    names(results)[names(results) == "id_col"] <- id
    if (notification == TRUE) {
      # Message
      message("Printing head results to console")
      print(head(results))
    }

    ## Summary file
    if (save_summary == TRUE) {
      # Write a summary txt file
      if (notification == TRUE) {
        # Message
        message("Creating summary file")
      }
      # Dim of raw_data, raw_data_dvc and raw_data_dvd
      cat("====================================================================", "\n", file = sum_file_name)
      cat(paste("Summary:", dataname), date(), file = sum_file_name, sep = "\n", append = TRUE)
      cat("====================================================================", "\n", file = sum_file_name, append = TRUE)
      cat(paste(dataname, "has", dim_raw_data1[1], "observations and", dim_raw_data1[2], "variables"),
          "\n", file = sum_file_name, append = TRUE)
      if (!is.null(keep_trials) | length(drop_vars) > 0) {
        cat(paste("* keep_trials:", keep_trials_sum, "drop_vars:", drop_col_sum), "\n", file = sum_file_name, append = TRUE)
        cat(paste("* After deleting unnecessary trials and variables according to keep_trials and drop_vars", dataname, "has", dim_raw_data2[1],
                  "observations and", dim_raw_data2[2], "variables"), "\n", file = sum_file_name, append = TRUE)

      }
      if (!is.null(dvc)) {
        cat(paste("* keep_trials_dvc:", keep_trials_dvc_sum), "\n", file = sum_file_name, append = TRUE)
        cat(paste("* After deleting unnecessary trials according to keep_trials_dvc", dvc, "data has", dim_raw_data_dvc[1], "observations and",
                  dim_raw_data_dvc[2], "variables"), "\n", file = sum_file_name, append = TRUE)
        if(!is.null(outlier_removal)) {
          cat(paste("* keep_trials_outlier:", keep_trials_outlier_sum), "\n", file = sum_file_name, append = TRUE)
          cat(paste("* After deleting unnecessary trials and variables according to keep_trials_outlier", dvc, "data for outlier removal procedures has",
                    dim(raw_data)[1], "observations and", dim(raw_data)[2], "variables"), "\n", file = sum_file_name, append = TRUE)
          cat(paste("*", outlier_name, "was calculated on", dim(raw_data)[1], "observations and", dim(raw_data)[2], "variables"), "\n",
              file = sum_file_name, append = TRUE)
        }
      } else {
        cat("* No dvc was found", "\n", file = sum_file_name, append = TRUE)
      }
      if(!is.null(dvd)) {
        cat(paste("* keep_trials_dvd:", keep_trials_dvd_sum), "\n", file = sum_file_name, append = TRUE)
        cat(paste("* After deleting unnecessary trials", dvd, "data has", dim_raw_data_dvd[1],
                        "observations and", dim_raw_data_dvd[2], "variables"), "\n", file = sum_file_name, append = TRUE)
      } else {
        cat("* No dvd was found", "\n", file = sum_file_name, append = TRUE)
      }
      # Dim of results_name
      cat(paste("*", results_name, "has", dim(results)[1], "ids and", dim(results)[2], "variables"), "\n",
          file = sum_file_name, append = TRUE)
      # Write ids
      cat("* ids:", levels(factor(results[[id]])), "\n", file = sum_file_name, append = TRUE)
      # Write levels of between-subjects indpendent variables
      if (length(between_vars) > 0) {
        cat("* Between-subject independent variables:", "\n",file = sum_file_name, append = TRUE)
        for (lev in 1:length(between_vars)) {
          cat(paste("  ", between_vars[[lev]], ": ", sep = ""), file = sum_file_name, append = TRUE)
          cat(levels(factor(raw_data[[between_vars[[lev]]]])), "\n", file = sum_file_name, append = TRUE)
        }
      } else {
        cat("* No between-subjects variables were found", "\n",file = sum_file_name, append = TRUE)
      }
      if (length(id_properties) > 0) {
        # Write names of id_properties
        cat("* id properties:", id_properties, "\n", file = sum_file_name, append = TRUE)
      } else {
        cat("* No id properties were found", "\n", file = sum_file_name, append = TRUE)
      }
      if (length(within_vars) == 0) {
        # Write to Summary file
        cat("* No within-subjects variables were found", "\n",file = sum_file_name, append = TRUE)
      }
      if (notification == TRUE) {
        message("Summary file finished")
      }
    }

    # Save results
    if (save_results == TRUE) {
      write.table(results, row.names = FALSE, file = results_name)
    }
    # Message
    if (save_results == TRUE) {
      message(paste(results_name, "has", dim(results)[1], "observations and", dim(results)[2], "variables"))
    }
    message("prep() returned a data frame to console")
    message("Hip Hip Hooray! prep() finished. Have a great day and may all your results be significant!")

    # Return
    return(results)

  # Checks if there are within-subject independent variables
  } else if (length(within_vars) > 0) {
    if (notification == TRUE) {
      if (length(between_vars) > 0) {
        # Message
        message("Your design is a mixed design (i.e., includes both between-subjects and within-subjects independent variables")
      } else {
        # Message
        message("Your design is a within-subject design")
      }
    }

    # Creats a temporay df to hold all dependent measures
    temp_dm <- cbind(id_col)

    ## All dependent measures will be calculated by id by within_condition
    # dvc
    # Calculates dependent measures for dvc in case dvc exists
    if (!is.null(dvc)) {

      # Change name of id column to "id"
      names(raw_data_dvc)[names(raw_data_dvc) == id] <- "id"
      # Change name of dvc column to "dvc"
      names(raw_data_dvc)[names(raw_data_dvc) == dvc] <- "dvc"

      # mdvc: mean dvc
      if (length(dm) == 0 | "mdvc" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating mean dvc by id by within_condition")
        }
        mdvc <- reshape2::dcast(raw_data_dvc, id ~ within_condition, mean,
                                value.var = "dvc")
        mdvc <- mdvc[-1]
        colnames(mdvc)[1:length(mdvc)] <- paste("mdvc", 1:length(mdvc), sep = "")
        # Round mdvc according to decimal_places
        mdvc <- round(mdvc, digits = decimal_places)
        # Bind
        temp_dm <- cbind(temp_dm, mdvc)
      }

      # sdvc: SD according to denominator n - 1
      if (length(dm) == 0 | "sdvc" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating SD for dvc by id by within_condition using denominator n")
        }
        # SD using denominator n - 1
        dvc_sd_r <- reshape2::dcast(raw_data_dvc, id ~ within_condition, sd,
                                    value.var = "dvc")
        dvc_sd_r <- dvc_sd_r[-1]
        # Length of each cell in the raw_data_dvc
        l_dvc <- reshape2::dcast(raw_data_dvc, id ~ within_condition, length,
                                 value.var = "dvc")
        l_dvc <- l_dvc[-1]
        # Calculating SD using denominator n
        sdvc <- dvc_sd_r * sqrt((l_dvc - 1) / (l_dvc))
        colnames(sdvc) <- paste("sdvc", 1:dim(sdvc)[2], sep = "")
        # Round sdvc according to decimal_places
        sdvc <- round(sdvc, digits = decimal_places)
        # Bind
        temp_dm <- cbind(temp_dm, sdvc)
      }

      # meddvc: median dvc
      if (length(dm) == 0 | "meddvc" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating median dvc by id by within_condition")
        }
        meddvc <- reshape2::dcast(raw_data_dvc, id ~ within_condition, median,
                                  value.var = "dvc", fill = NaN)
        meddvc <- meddvc[-1]
        colnames(meddvc) <- paste("meddvc", 1:dim(meddvc)[2], sep = "")
        # Round meddvc according to decimal_places
        meddvc <- round(meddvc, digits = decimal_places)
        # Bind
        temp_dm <- cbind(temp_dm, meddvc)
      }

      # Deviation from the mean for z scores
      dev_from_mean <- ave(raw_data_dvc$dvc, raw_data_dvc$id,
                           raw_data_dvc$within_condition,
                           FUN = function(x) x - mean(x))
      # sdvc for z scores
      sdvc_for_zscores <- ave(raw_data_dvc$dvc, raw_data_dvc$id,
                              raw_data_dvc$within_condition,
                              FUN = function(x) sd(x) * sqrt((length(x) - 1) / length(x)))
      # Z scores
      z_score <- dev_from_mean / sdvc_for_zscores
      raw_data_dvc$z_score <- z_score

      # tdvc: trimmed dvc according SD in sd_criterion
      # Create names for sd_criterion coulmn
      sd_criterion_names <- c()
      # Reset l to 1
      l <- 1
      while (l <= length(sd_criterion)) {
        sd_criterion_names[l] <- paste("t", sd_criterion[l], "dvc" ,sep = "")
        l <- l + 1
      }
      if (length(dm) == 0 | "tdvc" %in% dm) {
        # Create a data frame for sd_criterion
        tdvc_df <- data.frame(1:length(id_col))
        # Reset j to 1
        j <- 1
        # Calculating mean dvc by id by within_condition after rejecting values
        # above SD specified in sd_criterions according to sd_criterion
        while (j <= length(sd_criterion)) {
          if (notification == TRUE) {
            # Message
            message(paste("Calculating mean dvc by id by within_condition after rejecting values above", sd_criterion[j], "SD", sep = " "))
          }
          temp <- raw_data_dvc %>%
                  dplyr::group_by(id, within_condition) %>%
                  dplyr::summarise(n = mean(dvc[abs(z_score) < sd_criterion[j]])) %>%
                  reshape2::dcast(id ~ within_condition, value.var = "n")
          temp <- temp[-1]
          colnames(temp) <- paste(sd_criterion_names[j], 1:dim(temp)[2], sep = "")
          if (j == 1) {
            tdvc_df[2:(dim(temp)[2] + 1)] <- temp
          } else {
            tdvc_df[(dim(tdvc_df)[2] + 1):(dim(tdvc_df)[2] + dim(temp)[2])] <- temp
          }
          j <- j + 1
        }
        tdvc_df <- tdvc_df[-1]
        # Round tdvc_df according to decimal_places
        tdvc_df <- round(tdvc_df, digits = decimal_places)
        # Bind
        temp_dm <- cbind(temp_dm, tdvc_df)
      }

      # ntr: number of trimmed values for each SD in sd_criterion
      ntr_names <- c()
      # Reset l to 1
      l <- 1
      while (l <= length(sd_criterion)) {
        ntr_names[l] <- paste("n", sd_criterion[l], "tr" ,sep = "")
        l <- l + 1
      }
      # Create a data frame for ntr
      ntr_df <- data.frame(1:length(id_col))
      # Reset j to 1 for next while loop
      j <- 1
      # Counting number of rejected values for each sd_criterion
      while (j <= length(sd_criterion)) {
        if(length(dm) == 0 | "ntr" %in% dm) {
          if (notification == TRUE) {
            # Messgae
            message(paste("Counting number of rejected values for", sd_criterion[j], "SD", sep = " "))
          }
        }
        temp <- raw_data_dvc %>%
                dplyr::group_by(id, within_condition) %>%
                dplyr::summarise(n = length(dvc[abs(z_score) > sd_criterion[j]])) %>%
                reshape2::dcast(id ~ within_condition, value.var = "n")
        temp <- temp[-1]
        colnames(temp) <- paste(ntr_names[j], 1:dim(temp)[2], sep = "")
        if (j == 1) {
          ntr_df[2:(dim(temp)[2] + 1)] <- temp
        } else {
          ntr_df[(dim(ntr_df)[2] + 1):(dim(ntr_df)[2] + dim(temp)[2])] <- temp
        }
        j <- j + 1
      }
      ntr_df <- ntr_df[-1]
      if (length(dm) == 0 | "ntr" %in% dm) {
        # Bind
        temp_dm <- cbind(temp_dm, ntr_df)
      }

      # ndvc: number of values before removing values according to SD in
      # sd_criterion
      if (notification == TRUE) {
        # Message
        message("Counting number of values for dvc by id by within_condition before rejecting values according to sd_criterion")
      }
      ndvc <- reshape2::dcast(raw_data_dvc, id ~ within_condition, length,
                              value.var = "dvc")
      ndvc <- ndvc[-1]
      colnames(ndvc) <- paste("ndvc", 1:dim(ndvc)[2], sep = "")
      if (length(dm) == 0 | "ndvc" %in% dm) {
        # Bind
        temp_dm <- cbind(temp_dm, ndvc)
      }

      # ptr: proportion of rejected values according to SD in sd_criterion
      if (length(dm) == 0 | "ptr" %in% dm) {
        # Create names for ptr
        ptr_names <- c()
        # Reset l to 1
        l <- 1
        while (l <= length(sd_criterion)) {
          ptr_names[l] <- paste("p", sd_criterion[l], "tr" ,sep = "")
          l <- l + 1
        }
        # Create a data frame for ptr
        ptr_df <- data.frame(1:length(id_col))
        # Reset counter j for the next while loop
        j <- 1
        # Reset i to 1
        i <- 1
        while (j <= length(sd_criterion)) {
          if (notification == TRUE) {
            # Message
            message(paste("Calculating proportion of rejected values for dvc by id by within_condition for", sd_criterion[j], "SD", sep = " "))
          }
          if (j == 1) {
            temp <- ntr_df[i:dim(ndvc)[2]]
          } else {
            temp <- ntr_df[i:(i + dim(ndvc)[2] - 1)]
          }
          ptr <- temp / ndvc
          colnames(ptr) <- paste(ptr_names[j], 1:dim(ndvc)[2], sep = "")
          if (j == 1) {
            ptr_df[2:(dim(ndvc)[2] + 1)] <- ptr
            ptr_df <- ptr_df[-1]
          } else {
            ptr_df[(dim(ptr_df)[2] + 1):(dim(ptr_df)[2] + dim(ptr)[2])] <- ptr
          }
          i <- i + dim(ndvc)[2]
          j <- j + 1
        }
        # Round ptr_df according to digits = 3
        ptr_df <- round(ptr_df, digits = 3)
        # Bind
        temp_dm <- cbind(temp_dm, ptr_df)
      }

      # rminv: harmonic mean
      if (length(dm) == 0 | "rminv" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating harmonic mean for dvc by id by condition")
        }
        rminv <- reshape2::dcast(raw_data_dvc, id ~ within_condition,
                                 psych::harmonic.mean, value.var = "dvc")
        rminv <- rminv[-1]
        colnames(rminv) <- paste("rminv", 1:dim(rminv)[2], sep = "")
        # Round rminv according to decimal_places
        rminv <- round(rminv, digits = decimal_places)
        # Bind
        temp_dm <- cbind(temp_dm, rminv)
      }

      # pdvc: Percentiles
      if (length(dm) == 0 | "pdvc" %in% dm) {
        # Create names for pdvc
        percentails_names <- c()
        # Rest counter l
        l <- 1
        while (l <= length(percentiles)) {
          percentails_names[l] <- paste("p", percentiles[l], "dvc" ,sep = "")
          l <- l + 1
        }
        # Create a data frame for percentiles
        percentails_df <- data.frame(1:length(id_col))
        # Reset j to 1 before next while loop
        j <- 1
        while (j <= length(percentiles)) {
          if (notification == TRUE) {
            # Message
            message(paste("Calculating the", percentiles[j], "percentail for dvc by id by within_condition" ), sep = "")
          }
          temp <- reshape2::dcast(raw_data_dvc, id ~ within_condition,
                                  quantile, probs = percentiles[j],
                                  value.var = "dvc")
          temp <- temp[-1]
          colnames(temp) <- paste(percentails_names[j], 1:dim(temp)[2], sep = "")
          if (j == 1) {
            percentails_df[2:(dim(temp)[2] + 1)] <- temp
          } else {
            percentails_df[(dim(percentails_df)[2] + 1):(dim(percentails_df)[2] + dim(temp)[2])] <- temp
          }
          j <- j + 1
        }
        percentails_df <- percentails_df[-1]
        # Round percentails_df according to decimal_places
        percentails_df <- round(percentails_df, digits = decimal_places)
        # Bind
        temp_dm <- cbind(temp_dm, percentails_df)
      }

      ## Outlier removal procedures according to Van Selst & Jolicoeur (1994)
      # Check to see if the user wants an outlier removal procedure
      if (!is.null(outlier_removal)) {
        # User wants an outlier removal procedure

        # Save keep_trials_outlier as a string before parsing in order to later
        # use when writing Summary file
        keep_trials_outlier_sum <- keep_trials_outlier

        if (!is.null(keep_trials_outlier)) {
          # Delete unnecessary trials
          if (notification == TRUE) {
            # Message
            message("Keeping trials according to keep_trials_outlier and deleting unnecessary trials in raw_data")
            message("Outlier removal procedure will be calculated on the remaining trials")
          }
          # Subsetting trials
          keep_trials_outlier <- eval(parse(text = keep_trials_outlier))
          # Save raw_data after subsetting
          raw_data <- raw_data[keep_trials_outlier, ]
        } else {
          # Save keep_trials_outlier as "NULL" in order to later use when
          # writing Summary file
          keep_trials_outlier_sum <- class(keep_trials_outlier)
          # Save raw_data as raw_data
          raw_data <- raw_data
        }

        if (notification == TRUE) {
          # Message
          message("Calculating mean dvc by id by within_condition according to selected outlier removal procedure")
        }

        # Change name of id column to "id"
        names(raw_data)[names(raw_data) == id] <- "id"
        # Change name of dvc column to "dvc"
        names(raw_data)[names(raw_data) == dvc] <- "dvc"

        # Create final array that will contain trimmed means per id per within_condition
        final_data <- matrix(0, nrow = length(id_col), ncol = length(levels(raw_data$within_condition)))
        # Create final array that will contain number of trials trimmed per id per within_condition
        numtrials_data <- matrix(0, nrow = length(id_col), ncol = length(levels(raw_data$within_condition)))
        # Create final array that will contain percent of trials removed per id per within_condition
        percent_data <- matrix(0, nrow = length(id_col), ncol = length(levels(raw_data$within_condition)))
        # Create final array that will contain total number of trials per id per within_condition before removal
        totaltrials_data <- matrix(0, nrow = length(id_col), ncol = length(levels(raw_data$within_condition)))

        # Give name to columns according to the procedure
          if (outlier_removal == 1) {
            # Non recursive moving criterion procedure
            # nrmc: means for dvc per id per within_condition
            colnames(final_data) <- paste("nrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # nnrmc: number of trials trimmed per id per within_condition
            colnames(numtrials_data) <- paste("nnrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # pnrmc: percent of trials removed per id per within_condition
            colnames(percent_data) <- paste("pnrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # tnrmc: total number of trials per id per within_condition before removal
            colnames(totaltrials_data) <- paste("tnrmc", 1:length(levels(raw_data$within_condition)), sep = "")
          } else if (outlier_removal == 2) {
            # Modified recursive moving criterion procedure
            # mrmc: means for dvc per id per within_condition
            colnames(final_data) <- paste("mrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # nmrmc: number of trials trimmed per id per within_condition
            colnames(numtrials_data) <- paste("nmrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # pmrmc: percent of trials removed per id per within_condition
            colnames(percent_data) <- paste("pmrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # tmrmc: total number of trials per id per within_condition before removal
            colnames(totaltrials_data) <- paste("tmrmc", 1:length(levels(raw_data$within_condition)), sep = "")
          } else if (outlier_removal == 3) {
            # Hybrid recursive moving criterion procedure
            # hrmc: means for dvc per id per within_condition
            colnames(final_data) <- paste("hrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # nhrmc: number of trials trimmed per id per within_condition
            colnames(numtrials_data) <- paste("nhrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # phrmc: percent of trials removed per id per within_condition
            colnames(percent_data) <- paste("phrmc", 1:length(levels(raw_data$within_condition)), sep = "")
            # thrmc: total number of trials per id per within_condition before removal
            colnames(totaltrials_data) <- paste("thrmc", 1:length(levels(raw_data$within_condition)), sep = "")
          }

        for (i in 1:length(id_col)) {
          # For each id (i.e., subject)

          # Reset j to 1
          j <- 1

          # Do outlier removal procedure
          for (cond in levels(raw_data$within_condition)) {
            # For each within_condition

            # Isolate the current participant & condition's data
            temp_data <- raw_data[raw_data$id == id_col[i] & raw_data$within_condition == levels(raw_data$within_condition)[j], ]

            if (outlier_removal == 1) {
              # Non recursive procedure with moving criterion
              # Do current id trimming
              current_id <- non_recursive_mc(temp_data$dvc)
              final_data[i, j] <- current_id[1]
              percent_data[i, j] <- current_id[2]
              numtrials_data[i, j] <- current_id[3]
              totaltrials_data[i, j] <- current_id[4]
            } else if (outlier_removal == 2) {
              # Modified recursive procedure with moving criterion
              # Do current id trimming
              current_id <- modified_recursive_mc(temp_data$dvc)
              final_data[i, j] <- current_id[1]
              percent_data[i, j] <- current_id[2]
              numtrials_data[i, j] <- current_id[3]
              totaltrials_data[i, j] <- current_id[4]
            } else if (outlier_removal == 3) {
              # Hybrid recursive procedure with moving criterion
              # Do current id trimming
              current_id <- hybrid_recursive_mc(temp_data$dvc)
              final_data[i, j] <- current_id[1]
              percent_data[i, j] <- current_id[2]
              totaltrials_data[i, j] <- current_id[3]  # Should be 3. Do not change
            }
            j <- j + 1
          }
          # End of condition for loop
        }
        # End of id for loop

        # Round final data according to decimal_places
        final_data <- round(final_data, digits = decimal_places)
        # Round percent_data according to digits = 3
        percent_data <- round(percent_data, digits = 3)

        if (notification == TRUE) {
          # Messgae
          if (outlier_removal == 1) {
            message("Finished calculating means for dvc by id by within_condition according to non-recursive procedure with moving criterion ")
          } else if (outlier_removal == 2) {
            message("Finished calculating means for dvc by id by within_condition according to modified-recursive procedure with moving criterion")
          } else if(outlier_removal == 3) {
            message("Finished calculating means for dvc by id by within_condition according to hybrid procedure with moving criterion")
          }
        }
      }
      # End of outlier removal procedures

      # Finished dvc
      if (notification == TRUE) {
          # Message
          message("Finished calculating dependent measures for dvc")
      }
    }  # End of !is.null(dvc)

    # dvd
    # Check if dvd exists
    if (!is.null(dvd)) {

      # Change name of id column to "id"
      names(raw_data_dvd)[names(raw_data_dvd) == id] <- "id"
      # Change name of dvd column to "ac"
      names(raw_data_dvd)[names(raw_data_dvd) == dvd] <- "dvd"

      # mdvd: mean dvd
      if (length(dm) == 0 | "mdvd" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating mean dvd by id by within_condition")
        }
        mdvd <- reshape2::dcast(raw_data_dvd, id ~ within_condition, mean, value.var = "dvd")
        mdvd <- mdvd[-1]
        colnames(mdvd) <- paste("mdvd", 1:dim(mdvd)[2], sep = "")
        # Round mdvd according to decimal_dvd
        mdvd <- round(mdvd, digits = decimal_dvd)
        # Bind
        temp_dm <- cbind(temp_dm, mdvd)
      }

      # merr: mean error
      if (length(dm) == 0 | "merr" %in% dm) {
        if (notification == TRUE) {
          # Message
          message("Calculating mean error by id by within_condition")
        }
        merr <- 1 - mdvd
        colnames(merr) <- paste("merr", 1:dim(merr)[2], sep = "")
        # Round merr according to decimal_dvd
        merr <- round(merr, digits = decimal_dvd)
        # Bind
        temp_dm <- cbind(temp_dm, merr)
      }
      if (notification == TRUE) {
        # Message
        message("Finished calculating dvd data")
      }
    } # End of !is.null(dvd)

    # Remove the first column of the temporary df
    temp_dm <- temp_dm[-1]

    ## Creates results file
    if (notification == TRUE) {
      # Message
      message("Creating results")
    }
    if (length(id_properties) > 0) {
      # Bind also id_properties_df
      if (length(between_vars) > 0) {
        # Bind also between_bars_df
        results <- cbind(id_col, id_properties_df, between_vars_df, temp_dm)
      } else {
        results <- cbind(id_col, id_properties_df, temp_dm)
      }
    } else {
      if (length(between_vars) > 0) {
        # Bind also between_vars_df
        results <- cbind(id_col, between_vars_df, temp_dm)
      } else {
        results <- cbind(id_col, temp_dm)
      }
    }

    # In case of an outlier procedure, bind also results for that
    if (!is.null(outlier_removal)) {
      if (outlier_removal == 1 | outlier_removal == 2) {
        results <- cbind(results, final_data, percent_data, numtrials_data, totaltrials_data)
      } else if (outlier_removal == 3) {
        results <- cbind(results, final_data, percent_data, totaltrials_data)
      }
    }

    names(results)[names(results) == "id_col"] <- id
    if (notification == TRUE) {
      # Message
      message("Printing head results")
      print(head(results))
    }

    ## Summary
    if (save_summary == TRUE) {
      # Write a summary txt file
      if (notification == TRUE) {
        # Message
        message("Creating summary file")
      }
      # Dim of raw_data, raw_data_dvc and raw_data_dvd
      cat("====================================================================", "\n", file = sum_file_name)
      cat(paste("Summary:", dataname), date(), file = sum_file_name, sep = "\n", append = TRUE)
      cat("====================================================================", "\n", file = sum_file_name, append = TRUE)
      cat(paste(file_name, "has", dim_raw_data1[1], "observations and", dim_raw_data1[2], "variables"),
          "\n", file = sum_file_name, append = TRUE)
      if (!is.null(keep_trials) | length(drop_vars) >0) {
        cat(paste("* keep_trials:", keep_trials_sum, "drop_vars:", drop_col_sum), "\n", file = sum_file_name, append = TRUE)
        cat(paste("* After deleting unnecessary trials and variables", file_name, "has", dim_raw_data2[1],
                  "observations and", dim_raw_data2[2], "variables"), "\n", file = sum_file_name, append = TRUE)
      }
      if (!is.null(dvc)) {
        cat(paste("* keep_trials_dvc:", keep_trials_dvc_sum), "\n", file = sum_file_name, append = TRUE)
        cat(paste("* After deleting unnecessary trials and variables", dvc, "data has", dim_raw_data_dvc[1],
                        "observations and", dim_raw_data_dvc[2], "variables"), "\n", file = sum_file_name, append = TRUE)
        if (!is.null(outlier_removal)) {
          cat(paste("* keep_trials_outlier:", keep_trials_outlier_sum), "\n", file = sum_file_name, append = TRUE)
          cat(paste("* After deleting unnecessary trials and variables according to keep_trials_outlier", dvc, "data for outlier removal procedures has",
                    dim(raw_data)[1], "observations and", dim(raw_data)[2], "variables"), "\n", file = sum_file_name, append = TRUE)
          cat(paste("*", outlier_name, "was calculated on", dim(raw_data)[1], "observations and", dim(raw_data)[2], "variables"), "\n",
              file = sum_file_name, append = TRUE)
        }
      } else {
        cat("* No dvc was found", "\n", file = sum_file_name, append = TRUE)
      }
      if(!is.null(dvd)) {
        cat(paste("* keep_trials_dvd:", keep_trials_dvd_sum), "\n", file = sum_file_name, append = TRUE)
        cat(paste("* After deleting unnecessary trials and variables", dvd, "data has", dim_raw_data_dvd[1],
                  "observations and", dim_raw_data_dvd[2], "variables"), "\n", file = sum_file_name, append = TRUE)
      } else {
        cat("* No dvd was found", "\n", file = sum_file_name, append = TRUE)
      }
      # Dim of results_name
      cat(paste("*", results_name, "has", dim(results)[1], "ids and", dim(results)[2], "variables"), "\n",
          file = sum_file_name, append = TRUE)
      # Write ids
      cat( "* ids:", levels(factor(results[[id]])), "\n", file = sum_file_name, append = TRUE)
      # Write levels of between-subjects indpendent variables
      if (length(between_vars) > 0) {
        cat("* Between-subject independent variables:", "\n",file = sum_file_name, append = TRUE)
        for (lev in 1:length(between_vars)) {
          cat(paste("  ", between_vars[[lev]], ": ", sep = ""), file = sum_file_name, append = TRUE)
          cat(levels(factor(raw_data[[between_vars[[lev]]]])), file = sum_file_name, append = TRUE)
          cat("\n", file = sum_file_name, append = TRUE)
        }
      } else {
        cat("* No between-subjects variables were found", "\n",file = sum_file_name, append = TRUE)
      }
      if (length(id_properties) > 0) {
        # Write names of id_properties
        cat("* id properties:", id_properties, "\n", file = sum_file_name, append = TRUE)
      } else {
        cat("*", "No id properties were found", "\n", file = sum_file_name, append = TRUE)
      }
      if (notification == TRUE) {
        message("Summary file finished")
      }
      # Write levels of within-subject indpendent variables
      cat("* Within-subject independent variables:", "\n", file = sum_file_name, append = TRUE)
      for (lev in 1:length(within_vars)) {
        cat(paste("  ", within_vars[[lev]], ": ", sep = ""), file = sum_file_name, append = TRUE)
        cat(levels(factor(raw_data[[within_vars[[lev]]]])), file = sum_file_name, append = TRUE)
        cat("\n", file = sum_file_name, append = TRUE)
      }
    }

    # Save results
    write.table(results, row.names = FALSE, file = results_name)
    # Message
    if (save_results == TRUE) {
      message(paste(results_name, "has", dim(results)[1], "observations and", dim(results)[2], "variables"))
    }
    message("prep() returned a data frame to console")
    message("Hip Hip Hooray! prep() finished. Have a great day and may all your results be significant!")
    # Return
    return(results)
  }  # End of if (length(between_vars) > 0 & length(within_vars) == 0)
}  # End of prep()
