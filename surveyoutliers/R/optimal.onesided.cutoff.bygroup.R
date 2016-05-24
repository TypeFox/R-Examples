#' Optimal one-sided winsorization for survey outliers by group
#'
#' This function calculates optimal tuning parameter, cutoffs, and
#' winsorized values for one-sided winsorization, by group.
#'
#' @details
#' This function calculates optimal one-sided cutoffs for winsorization
#' where regression residuals are truncated at Q / (weight-1)
#' and Q satisfies the optimality result in Kokic and Bell (1994) and Clark (1995).
#' @param formula The regression formula (e.g. income ~ employment + old.turnover if income is survey variable and employment and old.turnover are auxiliary variables).
#' @param surveydata A data frame of the survey data including the variables in formula, piwt (inverse probability of selection) and gregwt (generalized regression estimator weight).
#' @param historical.reweight A set of reweighting factors for use when a historical dataset is being used.
#' It reweights from the historical sample to the sample of interest.
#' The default value of 1 should be used if the sample being used for optimising Q is the
#' same sample (or at least the same design) as the sample to which the winsorizing cutoffs are to be applied.
#' @param groupname The variable of this name in surveydata defines the groups for which Q is to be optimised.
#' If groupname is missing, it is assumed that cutoffs are to be optimised for the overall mean or total.
#' @param estimated.means.name The variable of this name in surveydata should contain an estimator of the expected values for each sample value of the variable of interest.
#'  If set to "", the regression model is estimated using IRLS.
#' @param stop Set to T to open a browser window (for debugging purposes)
#' @return A list consisting of Q.opt (the optimal Q), rlm.coef (the robust regression coefficients),
#' cutoffs (the winsorizing cutoffs for each unit in sample), y (the values of the variable of interest),
#' win1.values (the type 1 winsorized values of interest, i.e. the minimums of the cutoff and y)
#' win2.values (the type 2 winsorized values of interest, so that sum(surveydata$gregwt*win2.values) is the winsorized estimator)
#' @references
#' Clark, R. G. (1995), "Winsorisation methods in sample surveys," Masters thesis, Australian National University, http://hdl.handle.net/10440/1031.
#'
#' Kokic, P. and Bell, P. (1994), "Optimal winsorizing cutoffs for a stratified finite population estimator," J. Off. Stat., 10, 419-435.
#' @examples
#' test <- optimal.onesided.cutoff.bygroup(formula=y~x1+x2,
#' surveydata=survdat.example,groupname="industry")
#' plot(test$windata$y,test$windata$win1.values)

optimal.onesided.cutoff.bygroup <- function(formula,surveydata,historical.reweight=1,groupname,estimated.means.name="",stop=F){
  if(stop) browser()
  surveydata$ID <- c(1:nrow(surveydata))
  if(missing(groupname)){
    groupname <- "dummy"
    surveydata$dummy <- 1
  }
  windata <- NULL
  Qopt <- NULL
  for(g in sort(unique(surveydata[,groupname]))){
    surveydata.g <- surveydata[surveydata[,groupname]==g,]
    winresults.g <- optimal.onesided.cutoff( formula=formula , surveydata=surveydata.g ,
                                             historical.reweight=historical.reweight , estimated.means.name=estimated.means.name )
    Qopt <- c(Qopt,winresults.g)
    windata <- rbind(windata,winresults.g$windata)
  }
  names(Qopt) <- as.character(sort(unique(surveydata[,groupname])))
  windata <- windata[order(windata$ID),]
  list(Qopt=Qopt,windata=windata)
}
