#' Summary statistics for sample surveys
#' @description Wraps functions for summary statistics from survey package.
#' @param design an output form \code{\link{DesignSurvey}} function.
#' @param variables \code{\link{character}} \code{\link{vector}} with the type of estimate for each variable contained in \code{design} (see details).
#' @param conf.level the confidence level required.
#' @param rnd the number of decimal places (round) or significant digits (signif) to be used. If \code{NA}, scientific notation is used.
#' @return Matrix with survey summaries.
#' @details The length of \code{variables} must be equal to the length of \code{names(design$variables)} (see examples).
#' @references Lumley, T. (2011). Complex surveys: A guide to analysis using R (Vol. 565). Wiley.
#' 
#' \url{http://oswaldosantos.github.io/capm}
#' @export
#' @examples
#' # Load data.
#' data(psu.ssu)
#' data(survey.data)
#' 
#' ##########################################
#' ## Example 1 (two-stage cluster design) ##
#' ## General estimates                    ##
#' ##########################################
#' 
#' # Specify the two-stage cluster design.
#' design <- DesignSurvey(sample = survey.data, psu.ssu = psu.ssu,
#'                        psu.col = 2, ssu.col = 1, psu.2cd = 20)
#' 
#' # Look at the variables contained in the survey design
#' names(design$variables)
#' 
#' # Specify the type of estimate for each variable
#' variables <- c("total", "prop", "mean", rep("prop", 2),
#'                "total", rep("prop", 8))
#' 
#' # Make sure you specify the correct type of estimate for each variable
#' cbind(names(design$variables), variables)
#' 
#' # Calculate the summary statistics for the survey.
#' # Uncomment the following two lines (will take some seconds).
#' # estimates <- SummarySurvey(design, variables = variables, rnd = 3)
#' 
#' ##########################################
#' ## Example 2 (two-stage cluster design) ##
#' ## Sex-specific estimates               ##
#' ##########################################
#' 
#' # Make a copy of the dataset and select some 
#' # variables of interest.
#' sample1 <- survey.data[, c(1:4, 6:7, 11)]
#' 
#' # Transform to numeric the "sterilized" variable in order
#' # to estimate its total.
#' sample1[, 5] <- as.character(sample1[, 5])
#' sample1[which(sample1$sterilized == "yes"), 5] <- 1
#' sample1[which(sample1[, 5] == "no"), 5] <- 0
#' sample1[, 5] <- as.numeric(sample1[, 5])
#' 
#' # Define a survey design for each sex.
#' design.sex <- DesignSurvey(sample = sample1, psu.ssu = psu.ssu,
#'                            psu.col = 2, ssu.col = 1, psu.2cd = 20)
#' design.f <- subset(design.sex, sex == 'Female')
#' design.m <- subset(design.sex, sex == 'Male')
#' 
#' # Look at the variables contained in the survey design
#' names(design.sex$variables)
#'
#' # Specify the type of estimate for each variable
#' variables.sex <- c("total", "", "total", "prop", "prop")
#'
#' # Make sure you specify the correct type of 
#' # estimate for each variable
#' cbind(names(design.sex$variables), variables.sex)
#'
#' # Calculate the summary statistics for the surveys.
#' # Uncomment the following two lines (will take some seconds).
#' # estimates.f <- SummarySurvey(design.f, variables.sex, rnd = 3)
#' # estimates.m <- SummarySurvey(design.m, variables.sex, rnd = 3)
#' 
SummarySurvey <- function(design = NULL, variables = NULL, conf.level = 0.95, rnd = 3) {
  if (length(variables) != length(names(design$variables))) {
    stop('The length of variables argument must be equal to the length of names(design$variables)')
  }
  match1 <- names(design$variables)
  match2 <- c('psu.id', 'ssu.id', 'pop.size', 'psu.size')
  matches <- which(!is.na(match(match1, match2)))
  variables[matches] <- ''
  z <- abs(round(qnorm((1 - conf.level) / 2, 0, 1), 2))
  vrs <- design$variables
  out <- NULL
  for (i in 1:length(variables)) {
    if (variables[i] == 'total') {
      tmp <- svytotal(~ vrs[, i], design, na.rm = T, deff = T)
      tmp1 <- as.matrix(cbind(tmp, SE(tmp), confint(tmp), 
                              deff(tmp), cv(tmp) * z * 100), nr = 1)
      ci <- attributes(confint(tmp, level = conf.level))$dimnames[[2]]
      rownames(tmp1) <- paste0('Total.', names(vrs)[i])
      out <- rbind(out, tmp1)
    }
    if (variables[i] == 'mean') {
      tmp <- svymean(~ vrs[, i], design, na.rm = T, deff = T)
      tmp1 <- as.matrix(cbind(tmp, SE(tmp), confint(tmp), 
                              deff(tmp), cv(tmp) * z * 100), nr = 1)
      ci <- attributes(confint(tmp, level = conf.level))$dimnames[[2]]
      rownames(tmp1) <- paste0('Mean.', names(vrs)[i])
      out <- rbind(out, tmp1)
    }
    if (variables[i] == 'prop') {
      tmp <- svymean(~ vrs[, i], design, na.rm = T, deff = T)
      tmp1 <- as.matrix(cbind(tmp, SE(tmp), confint(tmp), 
                              deff(tmp), cv(tmp) * z * 100), nr = 1)
      ci <- attributes(confint(tmp, level = conf.level))$dimnames[[2]]
      rownames(tmp1) <- paste0('Prop.', rownames(tmp1))
      rownames(tmp1) <- gsub('vrs\\[, i\\]', paste0(names(vrs)[i], '.'), rownames(tmp1))
      out <- rbind(out, tmp1)
    }
  }
  colnames(out) <- c('Estimate', 'SE', ci[1], ci[2], 'Deff', 'Error (%)')
  if ('simple' %in% names(design)) {
    out <- out[ , -5]
  }
  ifelse (is.na(rnd), return(out), return(round(out, rnd)))
}
