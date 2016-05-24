# This page contains Utility functions


#' Utility Input Switch 
#'
#' Utility functions 
#'
#' @param pux The expected average criterion score of selected applicants.
#' @param uxs The average predicter score of those selected
#' @param rxy The correlation between the predictor composite and the criterion.
#' @param sr A selection ratio or a vector of selection ratios.
#' @return The expected average criterion score of selected applicants.
#' @note Many utility functions can except either (1) pux, (2) uxs and rxy, or
#'       (3) sr and rxy.
#' @author Allen Goebl and Jeff Jones
#' @keywords internal
#' @rdname internal.utilitySwitch
.utilitySwitch <- function(pux=NULL, uxs=NULL, rxy=NULL, sr=NULL){
    if(exists("pux", inherits=FALSE, mode="numeric")) {
        return(pux) 
    } else if(exists("uxs", inherits=FALSE, mode="numeric")) {
        return(uxs * rxy) 
    } else if(exists("sr", inherits=FALSE, mode="numeric")) {
        return(ux(sr) * rxy)
    } else { stop("Invalid or Missing Arguments: The expected average criterion
                  score could not be determined") }
}

#' Brogeden-Cronbach-Gleser Utility Model.
#'
#' Estimates the utility of an employee selection system.
#'
#' @param n The size of the applicant pool
#' @param sdy The standard deviation of performance in monetary units.
#' @param rxy The correlation between the predictor composite and the criterion.
#' @param uxs The average predicter score of those selected. If the uxs is 
#'  unknown, the sr argument can used instead.
#' @param sr A selection ratio or a vector of selection ratios.
#' @param pux The expected average criterion score of selected applicants
#' @param cost The cost per applicant of a selection system.
#' @param period The anticipated tenure of selected employees.
#' @return Estimated gain in utility.
#' @note This functions can except either (1) pux, (2) uxs and rxy, or (3) sr and rxy.
#' @author Allen Goebl and Jeff Jones
#' @references Cronbach, L. J., & Gleser, G. C. (1965). \emph{Psychological  
#' tests and personnel decisions.}, 37-40.
#' @examples
#' utilityBcg(sdy=10000, rxy=.50, sr=.30)
#' @export
utilityBcg <- function(n=1, sdy, rxy=NULL, uxs=NULL, sr=NULL, pux=NULL, cost=0, period=1) {
    pux <- .utilitySwitch(pux, uxs, rxy, sr)
    return((n * period * sdy * pux) - cost)
}

#' Boudreau Utility Model.
#'
#' This utility model extends the BCG model with additional financial variables.
#'
#' @param n The size of the applicant pool
#' @param sdy The standard deviation of performance in monetary units.
#' @param rxy the correlation between the predictor composite and the criterion.
#' @param uxs The average predicter score of those selected. If the uxs is 
#'  unknown, the sr argument can used instead.
#' @param sr A selection ratio or a vector of selection ratios.
#' @param pux The expected average criterion score of selected applicants.
#' @param cost The cost per applicant of a selection system.
#' @param period The anticipated tenure of selected employees.
#' @param v The proportion of new costs to new revenue (i.e. sc/sv).
#' @param tax The marginal tax rate.
#' @param i Discount rate.
#' @return Estimated gain in utility.
#' @note This functions can except either (1) pux, (2) uxs and rxy, or (3) sr and rxy.
#' @author Allen Goebl and Jeff Jones
#' @references Boudreau, J.W. (1983). Economic considerations in estimating
#'   the utility of human resource productivity improvement programs. 
#'   \emph{Personnel Psychology}, 36, 551-576.
#' @examples
#' utilityB(sdy=10000, rxy=.50, sr=.30, period=4, v=.5, tax=.1, i=.02)
#' @export
utilityB <- function(n=1, sdy, rxy=NULL, uxs=NULL, sr=NULL, pux=NULL, cost=0, period=1,
                     v=0, tax=0, i=0) {
    pux <- .utilitySwitch(pux, uxs, rxy, sr)
    fn <- function(x) ((sdy * (1-v) * (1-tax) * pux) / (1+i)^x)
    pv <- sum(fn(1:period))
    return((n * pv) - (cost * (1-tax)))
}

#' Schmidt-Hunter-Pearlman Utility Model.
#'
#' This model calculates the utility of an intervention accepting d rather
#' than rxy as an argument. 
#'
#' @param n The number of employees involved in the intervention.
#' @param sdy The standard deviation of performance in monetary units.
#' @param d The difference in job performance between the group recieving a 
#'  treatment and the group not recieving a treatment, expressed in
#'  standard deviation units. 
#' @param cost The cost of the intervention per participant.
#' @param period The anticipate duration of the training effect.
#' @return Estimated gain in utility.
#' @author Allen Goebl and Jeff Jones
#' @references Schmidt, F. L., Hunter, J. E., & Pearlman, K. (1982). Assessing
#'   the economic impact of personnel programs on workforce productivity.
#'   \emph{Personnel Psychology}, 35(2), 333-347.
#' @examples
#' utilityShp(sdy=10000, d=.50, period=4)
#' @export
utilityShp <- function(n=1, sdy, d, cost=0, period=1) {
    return((n * period * sdy * d) - cost)
}

#' Raju-Burke-Normand Utility Model
#'
#' This utility model uses SD of job performance ratings rather than the SD of
#' job performance in monetary units.
#'
#' @param n The size of the applicant pool.
#' @param sdr The standard deviation of ratings of job performance.
#' @param a The average total compensation.
#' @param rxy The correlation between the predictor composite and the criterion.
#' @param uxs The average predicter score of those selected. If the uxs is 
#'  unknown, the sr argument can used instead.
#' @param sr A selection ratio or a vector of selection ratios.
#' @param pux The expected average criterion score of selected applicants.
#' @param cost The cost per applicant of a selection system.
#' @param period The anticipated tenure of selected employees.
#' @return Estimated gain in utility.
#' @note This functions can except either (1) pux, (2) uxs and rxy, or (3) sr and rxy.
#' @author Allen Goebl and Jeff Jones
#' @references Raju, N.S., Burke, M.J. and Normand, J. (1990). A new approach
#' for utility analysis. \emph{Journal of Applied Psychology}, 75, 3-12.
#' @examples
#' utilityRbn(sdr=10000, a=90000, rxy=.50, sr=.30)
#' @export
utilityRbn <- function(n=1, sdr, a, rxy, uxs=NULL, sr=NULL, pux=NULL, cost=0, period=1) {
    pux <- .utilitySwitch(pux, uxs, rxy, sr)
    return((n * period * a * sdr * pux) - cost)
}

#' Raju-Cabrera-Lezotte Utility Model
#'
#' @keywords internal
#' @rdname todo
.rclUtility <- function() {}

#' Taylor-Russell Ratio
#'
#' Computes the Taylor Russel ratio
#'
#' @param rxy The correaltion between the predictor composite and the criterion.
#' @param sr The selection ratio.
#' @param br The base rate of the criterion. The cutoff point indicating
#'        success or failure.
#' @return The success ratio.
#' @author Allen Goebl and Jeff Jones
#' @references Taylor, H. C., & Russell, J. T. (1939). The relationship of
#' validity coefficients to the practical effectiveness of tests in selection:
#' Discussion and tables. \emph{Journal of Applied Psychology}, 25(5), 565.
#' @examples
#' trModel(rxy=.5, sr=.5, br=.6)
#' @export
trModel <- function(rxy, sr, br) {
    #Terms
    r_mat <- matrix(c(1, rxy, rxy, 1), nrow=2, ncol=2)
    qa_lower <- c(qnorm(1-sr), qnorm(1-br))
    qb_lower <- c(-Inf, qnorm(1-sr))
    qb_upper <- c(qnorm(1-br), Inf)
    #Density
    qa <- mvtnorm::pmvnorm(lower=qa_lower, mean=c(0,0), corr=r_mat)[1]
    qb <- mvtnorm::pmvnorm(lower=qb_lower, upper=qb_upper, mean=c(0,0), corr=r_mat)[1]
    return(qa / (qa + qb))
}

#' Taylor-Russell Utility
#'
#' The Taylor Russel Model can be used to estimate the utility of selecting
#' for a dichotomous criterion.
#'
#' @param n The size of the applicant pool.
#' @param rxy The correaltion between the predictor composite and the criterion.
#' @param sr The selection ratio.
#' @param br The criterion ratio. The cutoff point indicating success or failure.
#' @param dbr The monetary value of a 1 percent change in the basis rate per applicant.
#' @param cost The cost per applicant of a selection system.
#' @param period The anticipated tenure of selected employees
#' @return Estimated gain in utility.
#' @author Allen Goebl and Jeff Jones
#' @references Roomsburg (1989). \emph{Utility as a function of selection ratio
#' and base rate: An empirical investigation of military aviation selection.}
#' (Doctoral dissertation).
#' @examples
#' #trUtility(rxy=.5, sr=.5, br=.6, dbr=1000)
#' @keywords internal
#' @rdname internal.fuseRma
trUtility <- function(n=1, rxy, sr, br, dbr, cost=0, period=1) {
    br_change <- ((trModel(rxy=rxy, sr=sr, br=br) - br) * 100)
    return((n * period * dbr * br_change) - cost)
}
