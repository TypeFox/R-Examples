######################################################################################################################

# Function: Design.
# Argument: enroll.period,  enroll.dist, enroll.dist.par, followup.period, study.duration, dropout.dist,
#           dropout.dist.par
# Description: This function is used to create an object of class Design.
#' @export
Design = function(enroll.period = NULL,  enroll.dist = NULL, enroll.dist.par = NULL, followup.period = NULL, study.duration = NULL, dropout.dist = NULL, dropout.dist.par = NULL) {

  # Error checks
  if (!is.null(enroll.period) & !is.numeric(enroll.period)) stop("Design: enrollment period must be numeric.")
  if (!is.null(enroll.dist) & !is.character(enroll.dist)) stop("Design: enrollment distribution must be character.")
  if (!is.null(enroll.dist.par) & !is.list(enroll.dist.par)) stop("Design: enrollment distribution parameters must be provided in a list.")
  if (!is.null(followup.period) & !is.numeric(followup.period)) stop("Design: follow-up period must be provided in a list.")
  if (!is.null(study.duration) & !is.numeric(study.duration)) stop("Design: study duration must be provided in a list.")
  if (!is.null(dropout.dist) & !is.character(dropout.dist)) stop("Design: dropout distribution must be character.")
  if (!is.null(dropout.dist.par) & !is.list(dropout.dist.par)) stop("Design: enrollment distribution parameters must be provided in a list.")
  if (is.null(followup.period) & is.null(study.duration)) stop("Design: follow-up period or study duration must be defined")
  if (!is.null(followup.period) & !is.null(study.duration)) stop("Design: either follow-up period or study duration must be defined")
  if (is.null(enroll.dist) & !is.null(dropout.dist)) stop("Design: Dropout parameters cannot be specified without enrollment parameters.")

  design = list(enroll.period = enroll.period,
                enroll.dist = enroll.dist,
                enroll.dist.par = enroll.dist.par,
                followup.period = followup.period,
                study.duration = study.duration,
                dropout.dist = dropout.dist,
                dropout.dist.par = dropout.dist.par)


  class(design) = "Design"
  return(design)
  invisible(design)

}
