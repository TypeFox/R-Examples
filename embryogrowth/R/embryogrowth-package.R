#' Fit a parametric function that describes dependency of embryo growth to temperature
#'
#' \tabular{ll}{
#'  Package: \tab embryogrowth\cr
#'  Type: \tab Package\cr
#'  Version: \tab 6.2 - build 475\cr
#'  Date: \tab 2016-01-08\cr
#'  License: \tab GPL (>= 2)\cr
#'  LazyLoad: \tab yes\cr
#'  }
#' @title The package embryogrowth
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType package
#' @name embryogrowth-package
#' @description Tools to analyze the embryo growth and the sexualisation thermal reaction norms.\cr
#' The lastest version of this package can always been installed using:\cr
#' install.packages("http://www.ese.u-psud.fr/epc/conservation/CRAN/embryogrowth.tar.gz", repos=NULL, type="source")
#' @references Girondot, M. & Kaska, Y. 2014. A model to predict the thermal 
#'          reaction norm for the embryo growth rate from field data. Journal of
#'          Thermal Biology. 45, 96-102.
#' @seealso Delmas, V., Prevot-Julliard, A.-C., Pieau, C. & Girondot, M. 2008. 
#'          A mechanistic model of temperature-dependent sex determination 
#'          in a Chelonian, the European pond turtle. Functional 
#'          Ecology, 22, 84-93.
#' @seealso Girondot, M., Ben Hassine, S., Sellos, C., Godfrey, M. & Guillon, 
#'          J.-M. 2010. Modeling thermal influence on animal growth and sex 
#'          determination in Reptiles: being closer of the target gives new 
#'          views. Sexual Development, 4, 29-38.
#' @seealso Girondot, M. 1999. Statistical description of temperature-dependent 
#'          sex determination using maximum likelihood. Evolutionary Ecology 
#'          Research, 1, 479-486.
#' @seealso Girondot, M., & Kaska, Y. (2014). Nest temperatures in a loggerhead-
#'          nesting beach in Turkey is more determined by sea surface temperature 
#'          than air temperature. Journal of Thermal Biology, 47, 13-18.
#' @keywords Temperature Embryo Ecology Growth Gompertz Sex-determination
#' @examples
#' \dontrun{
#' library("embryogrowth")
#' packageVersion("embryogrowth")
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(115.758929130522, 428.649022170996, 503.687251738993, 
#' 12.2621455821612, 306.308841227278, 116.35048615105), .Names = c("DHA", 
#' "DHH", "DHL", "DT", "T12L", "Rho25"))
#' # or
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
#' out <- as.mcmc(result_mcmc_4p)
#' # This out obtained after as.mcmc can be used with coda package
#' # plot() can use the direct output of GRTRN_MHmcmc() function.
#' plot(result_mcmc_4p, parameters=1, xlim=c(0,550))
#' plot(result_mcmc_4p, parameters=3, xlim=c(290,320))
#' # summary() permits to get rapidly the standard errors for parameters
#' summary(result_mcmc_4p)
#' se <- result_mcmc_4p$SD
#' }

NULL
