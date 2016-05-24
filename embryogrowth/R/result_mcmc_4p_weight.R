#' Result of the mcmc using the nest database with weight
#' @title Result of the mcmc using the nest database
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name result_mcmc_4p_weight
#' @description Fit using the nest database
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage result_mcmc_4p_weight
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' w <- weightmaxentropy(formated, control_plot=list(xlim=c(20,36)))
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' # K or rK are not used for dydt.linear or dydt.exponential
#' resultNest_4p_weight <- searchR(parameters=x,  
#' 	fixed.parameters=pfixed, temperatures=formated,  
#' 	derivate=dydt.Gompertz, M0=1.7, test=c(Mean=39.33, SD=1.92),  
#' 	method = "BFGS", weight=w)
#' data(resultNest_4p_weight)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_weight, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_4p_weight <- GRTRN_MHmcmc(result=resultNest_4p_weight,  
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_4p_weight)
#' plot(result_mcmc_4p_weight, parameter="T12H", main="", xlim=c(290, 320), bty="n")
#' plotR(resultNest_4p_weight, SE=result_mcmc_4p_weight$SD, 
#'  ylim=c(0,0.3), las=1)
#' data(resultNest_4p)
#' data(result_mcmc_4p)
#' par(xpd=TRUE)
#' plotR(list(resultNest_4p_weight, resultNest_4p), 
#'  SE=list(result_mcmc_4p_weight$SD, result_mcmc_4p$SD), 
#'  ylim=c(0,0.4), las=1, col=list("red", "black"), 
#'  legend=list("Maximum entropy weighted", "Not weighted"))
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with 4 parameters and Gompertz model of growth weigted to maximized entropy
NULL
