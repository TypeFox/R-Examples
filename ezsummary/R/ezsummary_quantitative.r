#' Simple Summary for Quantitative/binary variables
#'
#' @description Function ezsummary_quantitative provides simple summary (Mean and
#' standard deviation with/without N) for quantitative data while function ezsummary_binary
#' provides simple summary (freq and percentage with/without total counts) for binary data.
#' These two function are simply wrappers outside of a summarise_each function. If
#' we just want to know the most basic statistical summary, this function can save us some
#' typing time. It also provide the option to include number of subjects inside the analyses.
#'
#' @param tbl The input matrix of data you would like to analyze.
#' @param n n is a True/False switch that controls whether counts(N) should be included in
#' the output
#' @param mean a T/F switch to control whether mean should be calculated
#' @param sd a T/F switch to control whether standard deviation should be calculated
#' @param sem a T/F switch to control whether standard error of the mean should be calculated
#' @param median a T/F switch to control whether median should be calculated
#' @param quantile a T/F switch to control whether 0%, 25%, 50%, 75% and 100% quantile should
#' be calculated
#' @param round.N Rounding Number
#' @param flavor Flavor has two possible inputs: "long" and "wide". "Long" is the default
#' setting which will put grouping information on the left side of the table. It is more
#' machine readable and is good to be passed into the next analytical stage if needed.
#' "Wide" is more print ready (except for column names, which you can fix in the next step,
#' or fix in LaTex or packages like \code{htmlTable}). In the "wide" mode, the analyzed
#' variable will be the only "ID" variable and all the stats values will be presented ogranized
#' by the grouping variables (if any). If there is no grouping, the outputs of "wide" and
#' "long" will be the same.
#' @param unit_markup When unit_markup is not NULL, it will call the ezmarkup function and
#' perform column combination here. To make everyone's life easier, I'm using the term "unit"
#' here. Each unit mean each group of statistical summary results. If you want to
#' know mean and stand deviation, these two values are your units so you can put something
#' like "[. (.)]" there
#'
#' @return It will return in the same format as a summarise_each function does
#'
#' @examples
#' library(dplyr)
#' mtcars %>% group_by(am) %>% select(mpg, wt, qsec) %>% ezsummary_quantitative()
#'
#' @export
ezsummary_quantitative <- function(tbl, n = FALSE, mean = TRUE, sd = TRUE, sem = FALSE, median = FALSE, quantile = FALSE, round.N=3, flavor = "long", unit_markup = NULL){
  # Define the following variable to avoid NOTE on RMD check
  variable = NULL
  # If the input tbl is a vector, convert it to a 1-D data.frame and set it as a 'tbl' (dplyr).
  if(is.vector(tbl)){
    tbl <- as.tbl(as.data.frame(tbl))
    attributes(tbl)$names <- "unknown"
    warning("ezsummary cannot detect the naming information from an atomic vector. Please try to use something like 'select(mtcars, gear)' to replace mtcars$gear in your code.")
  }

  if(flavor != "long" & flavor !="wide"){warning('The value of flavor has to be either "long" or "wide". Now the input is evalued as if you entered "long" by default. Please revise your function inputs!')}

  # Try to obtain grouping and variable information from the input tbl
  group.name <- attributes(tbl)$vars
  var.name <- attributes(tbl)$names
  if (!is.null(group.name)){
    group.name <- as.character(group.name)
    var.name <- var.name[!var.name %in% group.name]
    }
  n.group <- length(group.name)
  n.var <- length(var.name)

  # Set up the summarise_each formula based on the input.
  options <- c("N = length(na.omit(.))", "mean = round(mean(na.omit(.)), round.N)", "sd = round(sd(na.omit(.)), round.N)", "sem = round(sd(na.omit(.)) / sqrt(length(na.omit(.))), round.N)", "median = median(na.omit(.))", "q0 = quantile(.,names = F)[1], q25 = quantile(.,names = F)[2], q50 = quantile(.,names = F)[3], q75 = quantile(.,names = F)[4], q100 = quantile(.,names = F)[5]")
  option_names <- c("N", "mean", "sd", "sem", "median", "q0", "q25", "q50", "q75", "q100")
  option_switches <- c(n, mean, sd, sem, median, quantile)
  option_name_switches <- c(n, mean, sd, sem, median, quantile, quantile, quantile, quantile, quantile)
  calculation_formula <- paste0("summarise_each(tbl, funs(",
                                paste0(options[option_switches],collapse = ", "),
                                "))")

  # Perform the summarise_each calculation and pass the results to table_raw
  table_raw <- eval(parse(text=calculation_formula))

  # Standardize names when either n.var or the number of statistical analyses is small
    # When there is only one analysis (except(quantile))
    if (sum(option_name_switches) == 1){names(table_raw)[(n.group+1):ncol(table_raw)] <- paste(names(table_raw)[(n.group+1):ncol(table_raw)], option_names[option_name_switches], sep="_")}
    # When there is only one variable but the number of analyses is not 1
    if (n.var == 1 & sum(option_name_switches) != 1){names(table_raw)[(n.group+1):ncol(table_raw)] <- paste(var.name, names(table_raw)[(n.group+1):ncol(table_raw)], sep="_")}

  # Transform the raw table to a organized tabular output
  table_export <- melt(table_raw, id.vars = group.name) %>%
    separate(variable, into = c("variable", "statistics"), sep="_(?=[^_]*$)")
  # I have to stop in the middle and provide sort information to 'statistics' as
  # the dplyr::dcast function will automaticall sort all the values
  table_export$statistics <- factor(table_export$statistics, levels = c("N", "mean", "sd", "sem", "median", "q0", "q25", "q50", "q75", "q100"))

  # Use dplyr::dcast to reformat the table
  dcast_formula <- as.formula(paste0(c(group.name, "variable ~ statistics"), collapse = " + "))
  table_export <- table_export %>% dcast(dcast_formula)

  # Ezmarkup
  if(!is.null(unit_markup)){
    ezmarkup_formula <- paste0(paste0(rep(".", n.group), collapse = ""), ".", unit_markup)
    table_export <- ezmarkup(table_export, ezmarkup_formula)
  }

  # Turn the table from long to wide if needed
  if(flavor == "wide"){
    for(i in 1:n.group){
      table_export[,group.name[i]] <- paste(group.name[i], unlist(table_export[,group.name[i]]), sep=".")
    }
    table_export <- table_export %>% melt(id.var = c(group.name, "variable"), variable.name = "stats.var")
    dcast_formula <- paste0("dcast(table_export, variable ~ ", paste0(c(group.name, "stats.var"), collapse = " + "), ")")
    table_export <- eval(parse(text = dcast_formula))
  }

  # Fix the sorting of variable
  table_export$variable <- factor(table_export$variable, levels = var.name)
  table_export <- arrange(table_export, variable)

  attributes(table_export)$vars <- group.name
  attributes(table_export)$n.group <- n.group
  attributes(table_export)$n.var <- n.var
  attr(table_export, "class") <- c("tbl_df", "tbl", "data.frame")
  attr(table_export, "group_sizes") <- NULL
  attr(table_export, "biggest_group_size") <- NULL
  return(table_export)
}

