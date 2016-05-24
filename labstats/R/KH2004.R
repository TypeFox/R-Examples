#' @title Repeated measures data on rat muscles
#'
#' @description Repeated measures data set from Kristensen and Hansen
#' (2004) testing the effect of pinacidil on the force of muscle
#' contraction.
#'
#' @details Seven rats were euthanized and both soleus muscles (lower
#' leg) were extracted from each rat. One muscle from each pair was
#' randomised to the pinacidil condition and the other to the control
#' condition. Fifteen measurements of the force of muscle contraction
#' were taken every thirty seconds on each muscle.
#'
#' @format A data frame with 210 rows and 4 variables:
#' \describe{
#'   \item{time:}{Observation number, from 0 to 14. Observations were
#' taken every 30 seconds.}
#'   \item{cond:}{Condition, either Placebo or Pinacidil.}
#'   \item{rat:}{Rat identification number.}
#'   \item{values:}{Force of muscle contraction, normalised to the
#' first time point.}  }
#'
#' @references Kristensen M, Hansen T (2004). Statistical analyses of
#' repeated measures in physiological research: a tutorial. \emph{Adv
#' Physiol Educ.} 28(1-4):2-14.
#' 
#' @docType data
#' @name KH2004
NULL
