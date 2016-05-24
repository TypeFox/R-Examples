#' Result of the mcmc using the nest database
#' @title Result of the mcmc using the nest database
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name result_mcmc_6p
#' @description Fit using the nest database
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage result_mcmc_6p
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(115.758929130522, 428.649022170996, 503.687251738993, 
#' 12.2621455821612, 306.308841227278, 116.35048615105), .Names = c("DHA", 
#' "DHH", "DHL", "DT", "T12L", "Rho25"))
#' pfixed <- c(rK=2.093313)
#' resultNest_6p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_6p)
#' pMCMC <- TRN_MHmcmc_p(resultNest_6p, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_6p <- GRTRN_MHmcmc(result=resultNest_6p,  
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_6p)
#' plot(result_mcmc_6p, parameters="T12L", main="", xlim=c(290, 320), bty="n")
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with 6 parameters and Gompertz model of growth
NULL
