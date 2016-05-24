#' Create landmark data set
#' 
#' Create landmark data set from original data, which can be either in wide or
#' long format, see details.
#' 
#' For a given landmark time point \code{LM}, patients who have reached the
#' event of interest (outcome) or are censored before or at \code{LM} are
#' removed. Administrative censoring is applied at the time horizon.
#' Time-varying covariates are evaluated at the landmark time point \code{LM}.
#' Time-varying covariates can be specified in the \code{varying} item of the
#' \code{covs} argument, in two ways. In the first way (data in long format)
#' different values of time-dependent covariate(s) are stored different rows of
#' the data, with \code{id} identifying which values belong to the same
#' subject; the column specified through \code{rtime} then contains the time
#' points at which the value of the covariate changes value; with
#' \code{right=TRUE} (default), it is assumed that the covariate changes value
#' at the time point specified in \code{rtime} (and hence is not used for
#' prediction of an event at \code{rtime}), while with \code{right=FALSE}, it
#' is assumed that the covariate changes value just before the time point
#' specified in \code{rtime}. The second way (data in wide format) can only be
#' used for a specific type of time-varying covariates, often used to model
#' whether some other event has occurred or not, namely those that change value
#' from 0 (event not yet occurred) to 1 (event has occurred).
#' 
#' @param data Data frame from which to construct landmark dataset
#' @param outcome List with items \code{time} and \code{status}, containing
#' character strings identifying the names of time and status variables,
#' respectively, of the survival outcome
#' @param LM Scalar, the value of the landmark time point
#' @param horizon Scalar, the value of the horizon. Administrative censoring is
#' applied at \code{horizon}.
#' @param covs List with items \code{fixed} and \code{varying}, containing
#' character strings specifying column names in the data containing time-fixed
#' and time-varying covariates, respectively
#' @param format Character string specifying whether the original data are in
#' wide (default) or in long format
#' @param id Character string specifying the column name in \code{data}
#' containing the subject id; only needed if \code{format="long"}
#' @param rtime Character string specifying the column name in \code{data}
#' containing the (running) time variable associated with the time-varying
#' covariate(s); only needed if \code{format="long"}
#' @param right Boolean (default=\code{TRUE}), indicating if the intervals for
#' the time-varying covariates are closed on the right (and open on the left)
#' or vice versa, see \code{\link{cut}}
#' @return A landmark data set, containing the outcome and the values of
#' time-fixed and time-varying covariates taken at the landmark time points.
#' The value of the landmark time point is stored in column \code{LM}.
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#' Clinical Survival Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' test0 <- data.frame(id=c(1,1,1,2,2,2),survyrs=c(2.3,2.3,2.3,2.7,2.7,2.7),
#'   survstat=c(1,1,1,0,0,0),age=c(76,76,76,68,68,68),gender=c(1,1,1,2,2,2),
#'   bp=c(80,84,88,92,90,89),bptime=c(1,2,2.2,0,1,2))
#' cutLM(test0, outcome=list(time="survyrs", status="survstat"),
#'   LM=1, horizon=2.5, covs=list(fixed=c("age","gender"),varying="bp"),
#'   format="long", id="id", rtime="bptime")
#' # Note how the previous example does not use the value of the time-varying
#' # covariate AT time=LM, only just before (if available). This is in line
#' # with the time-varying covariates being predictable.
#' # If you want the value of the time-varying covariate at time=LM if it
#' # changes value at LM, then use right=FALSE
#' cutLM(test0, outcome=list(time="survyrs", status="survstat"),
#'   LM=1, horizon=2.5, covs=list(fixed=c("age","gender"),varying="bp"),
#'   format="long", id="id", rtime="bptime", right=FALSE)
#' 
#' # An example of a time-varying covariate in wide format; recyrs and recstat
#' # are time and status of a (cancer) recurrence. Here it is assumed that the
#' # value of the time-varying covariate is 0 and changes value to 1 at recyrs.
#' # The status variable is not used!
#' test1 <- data.frame(id=1:4,survyrs=c(7.6,8.4,5.3,2.6),survstat=c(0,1,1,0),
#'   age=c(48,52,76,18),gender=c(1,2,2,1),recyrs=c(7.6,5.2,0.8,2.6),
#'   recstat=c(0,1,1,0))
#' cutLM(test1, outcome=list(time="survyrs", status="survstat"),
#'   LM=3, horizon=8, covs=list(fixed=c("id","age","gender"),varying="recyrs"))
#' 
#' # The same example in long format, similar to (but not the same as) the way
#' # one would use a time-varying covariate in long format.
#' test2 <- data.frame(id=c(1,2,2,3,3,4),survyrs=c(7.6,8.4,8.4,5.3,5.3,2.6),
#'   survstat=c(0,1,1,1,1,0),age=c(48,52,52,76,76,18),gender=c(1,2,2,2,2,1),
#'   rec=c(0,0,1,0,1,0),rectime=c(0,0,5.2,0,0.8,0))
#' cutLM(test2, outcome=list(time="survyrs", status="survstat"),
#'   LM=3, horizon=8, covs=list(fixed=c("age","gender"),varying="rec"),
#'   format="long", id="id", rtime="rectime")
#' 
#' @export cutLM
cutLM <- function(data, outcome, LM, horizon, covs,
  format = c("wide","long"), id, rtime, right=TRUE)
{
    format <- match.arg(format)
    if (format=="wide") {
      LMdata <- data
      if (!is.null(covs$varying))
        LMdata[[covs$varying]] <- 1 - as.numeric(LMdata[[covs$varying]] > LM)
    } else {
      if (missing(id))
        stop("argument 'id' should be specified for long format data")
      if (missing(rtime))
        stop("argument 'rtime' should be specified for long format data")
      ord <- order(data[[id]],data[[rtime]])
      data <- data[ord,]
      ids <- unique(data[[id]])
      n <- length(ids)
      # initialize LMdata; copy first row of each subject
      LMdata <- data[which(!duplicated(data[[id]])),]
      for (i in 1:n) {
          wh <- which(data[[id]]==ids[i])
          di <- data[wh,]
          idx <- cut(LM,c(data[[rtime]][wh],Inf),right=right,labels=FALSE)
          if (!is.na(idx)) LMdata[i,] <- di[idx,]
          else {
            LMdata[i,] <- di[1,]
            LMdata[[covs$varying]][i] <- NA
            LMdata[[rtime]][i] <- NA
          } 
      }
    }
    LMdata <- LMdata[LMdata[[outcome$time]] > LM,]
    if (format=="long") LMdata <- LMdata[!is.na(LMdata[[id]]),]
    # apply administrative censoring at horizon
    LMdata[outcome$status] <- LMdata[[outcome$status]] *
      as.numeric(LMdata[[outcome$time]] <= horizon)
    LMdata[outcome$time] <- pmin(as.vector(LMdata[[outcome$time]]),horizon)
    LMdata$LM <- LM
    if (format=="long")
      cols <- match(c(id,outcome$time,outcome$status,covs$fixed,covs$varying,rtime,"LM"),
        names(LMdata))
    else
      cols <- match(c(outcome$time,outcome$status,covs$fixed,covs$varying,"LM"),
        names(LMdata))
    return(LMdata[,cols])
}
