#' Result of the mcmc using the nest database
#' @title Result of the mcmc using the nest database
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name result_mcmc_4p
#' @description Fit using the nest database
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage result_mcmc_4p
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(118.431040984352, 498.205702157603, 306.056280989839, 
#' 118.189669472381), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_4p <- GRTRN_MHmcmc(result=resultNest_4p,  
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_4p)
#' plot(result_mcmc_4p, parameters="T12H", main="", xlim=c(290, 320), bty="n")
#' plotR(resultNest_4p, SE=result_mcmc_4p$SD, ylim=c(0,0.4), las=1)
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with 4 parameters and Gompertz model of growth
NULL
