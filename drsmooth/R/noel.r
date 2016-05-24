#' @name noel
#' @aliases noel
#' @title No/Lowest Observed Effect Levels
#' @usage
#' noel(dosecolumn = "", 
#'       targetcolumn = "", 
#'       data_type = "continuous",
#'       tests = c("all"),
#'       alternatives = c("greater", "two.sided"),
#'       alpha = .05,
#'       control = "",
#'       tot.obs = "",
#'       data = NA)
#' @description
#' This function calculates and displays the results of the requested no/lowest observed effect level tests.
#' @details
#' Dosecolumn should be assigned the name of the dose column in the input dataframe.
#' Targetcolumn should be assigned the name of the response column in the input dataframe.
#' Tests is defaulted to run all tests available, given the data_type defined.
#' 
#' If the data_type is defined as "continuous", Dunnett's, Dunn's, and Dunnett's T3 tests
#' are available and are all executed (default), unless tests is defined as a subset of the 
#' list: c("dunnetts", "dunns", "dunnettst3"). The input data frame is expected to include
#' a numeric column containing the dose and a numeric column containing a continuous
#' response variable.
#' 
#' When data_type is defined as "dichotomous" for a dichotomous response variable,
#' the Fisher's exact test is available and executed, with the tests parameter is either left
#' at the c("all") default or specified as c("fishers.exact").  The control dose is defaulted
#' to the lowest observed dose, unless a different control dose is provided as a string to
#' the control parameter.  The input data frame is expected to include summarized data of
#' dichotomous outcome tests:  one row for each dose, the total number tested at that dose,
#' and the total number of events observed at that dose.
#' 
#' The alternative parameter specifies the direction(s) of the alternative hypothesis.
#' All alternatives listed will be tested.
#' 
#' The alpha level determining significance can be specified.
#' 
#' @param dosecolumn  Character string, name of dose column to be tested.
#' @param targetcolumn  Character string, name of response column to be tested.
#' @param data_type  Allowed values "continuous" (default) or "dichotomous".
#' @param tests  Available tests depend on data_type.  See details.  
#' @param alternatives  Character string(s) specifying the direction of the alternative hypothesis.
#' Must be one or more of "greater", "two.sided", or "less".
#' @param alpha  Significance level (numeric) to be used.
#' @param control Level of dose to be used as the control for dichotomous data.
#' @param tot.obs Character string, column with number tested at each dose for dichotomous outcomes.
#' @param data  Input dataframe.  See details for expected formats.
#' @return
#' Tables are printed giving the comparisons of the active dose levels to the zero dose control
#' along with indications of significance specific to each type of test.
#' @examples
#' # Prints all available tests of no/lowest observed effect levels at default alpha=.05:
#' data(DRdata)
#' noel("dose", "MF_Log", data=DRdata) 
#'
#' # Dunnett's T3 tests at user-specified alpha of .01:
#' data(DRdata)
#' noel("dose", "MF_Log", tests=c("dunnettst3"), alpha=.01, data=DRdata) 
#' 
#' # Fisher's exact test for dichotomous outcome data:
#' data(DIdata)
#' noel(dosecolumn = "Dose", 
#'      targetcolumn = "Tumor", 
#'      data_type = "dichotomous", 
#'      tot.obs = "n", 
#'      data = DIdata)
#' @export

noel <- function (dosecolumn = "", 
                  targetcolumn = "", 
                  data_type = "continuous", 
                  tests = c("all"), 
                  alternatives=c("greater", "two.sided"), 
                  alpha=0.05, 
                  control = "",
                  tot.obs = "",
                  data=NA) {
  
  # future: call validate function from here, parameter from = noel.
  # future: add test for event that user enters wrong data type, given actual data entered
  
  if (data_type == "continuous") {
    allowed_tests <- c("dunnetts", "dunns", "dunnettst3")} else {
      if (data_type == "dichotomous") {
        allowed_tests <- c("fishers.exact")} else {
          stop("data_type must be specified as 'continuous' or 'dichotomous'")
        }
    }
  
  if (tests == "all") {tests <- allowed_tests} else {
    if (tests %in% allowed_tests) {} else {stop("test is not allowed for given data type")}
  }
  
  f <- get("dosefactor", envir = environment(drsmooth))
  data <- f(dosecolumn, data)
  
  for (i in 1:length(tests)) {
    test_name <- tests[i]
    if ((test_name == "dunns" | test_name == "dunnettst3") & "less" %in% alternatives) {
      alternatives <- gsub('less', NA, alternatives)
      alternatives <- alternatives[!is.na(alternatives)]
    }
    for (j in 1:length(alternatives)) {
      alternative <- alternatives[j]
      if (alternative == "greater") {label <- "Direction-Greater, One-Tailed"} else {
        if (alternative == "two.sided") {label <- "Two-Tailed"} else {
          if (alternative == "less") {label <- "Direction-Less, One-Tailed"} else {
            next}
        }
      }
      f <- get(test_name, envir = environment(drsmooth))
      result <- f(targetcolumn, alternative, alpha, control, tot.obs, label, data)
      if (is.null(result)) {next} else {
        suppressWarnings(for (k in 1:length(result)) {
          if (class(result[[k]])=="matrix") {drsmooth.print(result[[k]])} else {print(result[[k]])}
        }
        )
      }
    }
  }
}
