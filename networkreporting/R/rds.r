#####################################################
## rds.r
##
## limited implementation of RDS estimator

#####################################################
##' rdsII.estimator
##'
##' compute an estimate for the prevalence of a trait
##' from an RDS sample, using the estimator described
##' in TODO [Volz + Heckathorn '08]
##'
##' NOTE: we have no weights for now, right? RDS doesn't
##' get used with weights?
##'
##' @param survey.data the dataframe with RDS survey results
##' @param d.hat.vals the name or index of the column that contains
##'                  each respondent's estimated degree
##' @param y.vals the name or index of the column that contains
##'               the quantity of interest. if this is a
##'               dichotomous trait, it should be 0 / 1
##' @param missing if "ignore", then proceed with the analysis without
##'                doing anything about missing values. if "complete.obs"
##'                then only use rows that have no missingness for the
##'                computations (listwise deletion). care
##'                must be taken in using this second option
##' @param verbose if TRUE, print messages to the screen
##' @return the RDS-II estimate of the average of the quantity of interest
##' @export
rdsII.estimator <- function(survey.data,
                            d.hat.vals,
                            y.vals,
                            missing="ignore",
                            verbose=FALSE)
{

  if (! missing %in% c("ignore", "complete.obs")) {
    stop("error in specifying procedure for handling missing values in rdsII.estimator. invalid option.\n")
  }

  d.hat.vals <- surveybootstrap:::get.var(survey.data, d.hat.vals)

  y.vals <- surveybootstrap:::get.var(survey.data, y.vals)

  #### compute the actual estimates
  num <- y.vals/d.hat.vals
  denom <- 1/d.hat.vals

  ## figure out if we have to only use non-missing entries
  touse.idx <- 1:length(y.vals)
  if (missing == "complete.obs") {
    touse.idx <- which( (! is.na(y.vals) & ! is.na(d.hat.vals)))

    notused <- length(y.vals) - length(touse.idx)

    if (verbose & notused > 0) {
      cat(str_c("missing=='complete.obs', so dropping ",
                notused, " rows with missing data\n"))
    }

  }

  ## NB: for now, this will return NA if either the degrees or
  ##     the y_i's has any NAs and missing == "ignore"
  res <- sum(num[touse.idx])/sum(denom[touse.idx])

  return(res)

}

