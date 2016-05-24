#' likelihood_phenology estimate likelihood for a set of parameters.
#' @title Estimate the likelihood of timeseries based on a set of parameters.
#' @author Marc Girondot
#' @return The likelihood of the data with the parameters
#' @param data Dataset generated with add_format
#' @param parametersfixed Set of fixed parameters
#' @param parametersfit Set of parameters to be fitted
#' @param method_incertitude 2 [default] is the correct one from a statistical point of view; \cr
#'                           0 is an aproximate method more rapid; \cr
#'                           1 is an alternative more rapid but biased.
#' @param zero_counts example c(TRUE, TRUE, FALSE) indicates whether the zeros have 
#'                    been recorder for each of these timeseries. Defaut is TRUE for all.
#' @param result An object obtained after fit_phenology()
#' @description This function is used to estimate the likelihood based on a set of parameters.
#' @examples
#' \dontrun{
#' # Read a file with data
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' data(Gratiot)
#' # Generate a formated list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Estimate likelihood with this initial set of parameters
#' likelihood_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL)
#' # Or directly from a result object
#' likelihood_phenology(result=result_Gratiot)
#' }
#' @export

likelihood_phenology <-
function(data=NULL, parametersfit=NULL, parametersfixed=NULL, zero_counts=NULL, method_incertitude=NULL, result=NULL) {

# if result est donné, on prend les données dedans et on remplace celles introduites en plus

if (!is.null(result)) {
  if (class(result) != "phenology") {
    stop("The object result must be the result of a fit_phenology()")
  }

  if (is.null(data)) {data <- result$data}
  if (is.null(parametersfit)) {parametersfit <- result$par}
  if (is.null(parametersfixed)) {parametersfixed <- result$parametersfixed}
  if (is.null(zero_counts)) {zero_counts <- result$zero_counts}
  if (is.null(method_incertitude)) {method_incertitude <- result$method_incertitude}

}

# if (!is.null(data)) {
#   if (class(data) != "phenologydata") {
#     stop("The data object must be the result of a add_format()")
#   }
# }


# if (is.null(parametersfixed)) {parametersfixed <- NA}
if (is.null(method_incertitude)) {method_incertitude <- 0}
if (is.null(zero_counts)) {zero_counts <- TRUE}
	
	if (length(zero_counts)==1) {zero_counts <- rep(zero_counts, length(data))}
	if (length(zero_counts)!=length(data)) {
		stop("zero_counts parameter must be TRUE (the zeros are used for all timeseries) or FALSE (the zeros are not used for all timeseries) or possess the same number of logical values than the number of series analyzed.")
	}
	
LnL <- getFromNamespace(".Lnegbin", ns="phenology")(parametersfit, pt=list(data=data, fixed=parametersfixed, incertitude=method_incertitude, zerocounts=zero_counts))
	
	return(LnL)

}
