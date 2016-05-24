WCE.data.frame <- function(data, analysis, nknots, cutoff, constrained = FALSE, int.knots = NULL, aic = FALSE, MatchedSet = NULL, id='Id', event = 'Event',  start='Start', stop='Stop', expos ='dose', covariates = NULL, controls = NULL, ...) {
   if (constrained == 'right') constrained = 'Right'
   if (constrained == 'left') constrained = 'Left'
   if (is.data.frame(data) == F)  stop("ERROR: data must be a data frame")
   if (is.null(covariates) == F & sum(covariates %in% names(data))!= length(covariates)) stop("ERROR: At least one covariate does not belong to the data set supplied") 
   if (analysis == 'Cox' | analysis == 'cox') {WCE.cox(data, nknots, cutoff, constrained, int.knots = NULL, aic, id, event, start, stop, expos, covariates, controls)} else 
   if (analysis == 'NCC' | analysis == 'ncc') {stop("Methods for nested case control designs are not implemented yet.")} else 
   if  (analysis == 'CC' | analysis == 'cc') {stop("Methods for case control designs are not implemented yet")} else stop("ERROR: Requested analysis is invalid")
 }