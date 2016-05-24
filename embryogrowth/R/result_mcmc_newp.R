#' Result of the mcmc using the nest database with anchored parameters
#' @title Result of the mcmc using the nest database with anchored parameters
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name result_mcmc_newp
#' @description Fit using the nest database with anchored parameters
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage result_mcmc_newp
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' newp <- GenerateAnchor(nests=formated, number.anchors=7)
#' pfixed <- c(rK=2.093313)
#' resultNest_newp <- searchR(parameters=newp, fixed.parameters=pfixed, 
#'   temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#'   test=c(Mean=39.33, SD=1.92))
#' data(resultNest_newp)
#' pMCMC <- TRN_MHmcmc_p(resultNest_newp, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_newp <- GRTRN_MHmcmc(result=resultNest_newp,  
#'   parametersMCMC=pMCMC, n.iter=1000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_newp)
#' data(resultNest_4p)
#' newp <- GenerateAnchor(nests=resultNest_4p, number.anchors=7)
#' # Here the confidence interval is built based on anchored parameters
#' plotR_hist(resultNest_4p, parameters=newp, SE=result_mcmc_newp$SD, 
#'  ylim=c(0,0.4), ylimH=c(0,0.4))
#' # Here the confidence interval is built based on parametric SSM equation
#' data(result_mcmc_4p)
#' plotR_hist(resultNest_4p, SE=result_mcmc_4p$SD,
#'  ylim=c(0,0.4), ylimH=c(0,0.4))
#' plot(result_mcmc_newp, las=1, xlim=c(0,30), parameters="294", 
#' breaks=c(0, 1.00095, 2.0009, 3.00085, 4.0008, 5.00075, 6.0007, 7.00065, 8.0006, 9.00055, 
#' 10.0005, 11.00045, 12.0004, 13.00035, 14.0003, 15.00025, 16.0002, 17.00015, 18.0001, 
#' 19.00005, 20))
#' plot(result_mcmc_newp, las=1, xlim=c(0,30), parameters="296.333333333333")
#' plot(result_mcmc_newp, las=1, xlim=c(0,30), parameters=3)
#' # Confidence interval based on quantiles
#' plotR_hist(resultNest_4p, parameters=newp, SE=NULL, 
#'  ylim=c(0,0.4), ylimH=c(0,0.4))
#' CI <- apply(result_mcmc_newp$resultMCMC[[1]], 2, quantile, probs=c(0.025, 0.5, 0.975))
#' plot_add(as.numeric(colnames(CI))-273.15, CI[1,], lty=2, type="l")
#' plot_add(as.numeric(colnames(CI))-273.15, CI[2,], lty=2, type="l")
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with anchored parameters
NULL
