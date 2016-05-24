#' STRN estimates the parameters that best describe the sexualisation thermal reaction norm within the TSP
#' @title Estimate the parameters that best describe the sexualisation thermal reaction norm within the TSP
#' @author Marc Girondot
#' @return The list with object return by optim() 
#' @param Initial_STRN Values for initial model of Sexualisation Thermal Reaction Norm
#' @param EmbryoGrowthTRN The Embryo growth Thermal Reaction Normal obtained with searchR()
#' @param tsd The model used to predict sex ratio, obtained from tsd()
#' @param Sexed The number of sexed embryos
#' @param Males The number of males embryos
#' @param Females The number of females embryos
#' @param Temperatures The temperature from out of info.nests to be used
#' @param ... Parameters used for optim()
#' @description Estimate the parameters that best describe the sexualisation thermal reaction norm within the TSP.\cr
#' The Temperatures parameter is a character string which can be:\cr
#' \itemize{
#'   \item \code{TimeWeighted.temperature.mean}
#'   \item \code{TSP.TimeWeighted.temperature.mean}
#'   \item \code{TSP.MassWeighted.temperature.mean}
#'   \item \code{TSP.STRNWeighted.temperature.mean}
#'   \item \code{TSP.MassWeighted.STRNWeighted.temperature.mean}
#'   \item \code{MiddleThird.TimeWeighted.temperature.mean}
#'   }
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' MedIncubation_Cc <- subset(STSRE_TSD, Species=="Caretta caretta" & 
#' RMU=="Mediterranean" & Sexed!=0)
#' Med_Cc <- with(MedIncubation_Cc, tsd(males=Males, females=Females, 
#'  temperatures=Incubation.temperature, par=c(P=29, S=-0.01)))
#' plot(Med_Cc, xlim=c(25, 35))
#' # Initial_STRN <- rep(1, 7)
#' # names(Initial_STRN) <- as.character(seq(from=20, to=35, length=7))
#' Initial_STRN <- structure(c(1, 143.248982215757, -25.7029976477549, -0.00489843027318209,
#' -8.94560833594928, 135.781961273868, 71.2176230826628), 
#' .Names = c("20", "22.5", "25", "27.5", "30", "32.5", "35"))
#' males <- males <- c(7, 0, 0, 0, 0, 5, 6, 3, 5, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' sexed <- rep(10, length(males))
#' fitSTRN <- STRN(Initial_STRN, EmbryoGrowthTRN=resultNest_4p, tsd=Med_Cc, 
#' Sexed=sexed, Males=males, 
#' Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean")
#' CTE <- info.nests(NestsResult=resultNest_4p, 
#' SexualisationTRN=fitSTRN$par, out="summary")
#' plot_add(x=CTE$TSP.MassWeighted.STRNWeighted.temperature.mean, y=males/sexed, 
#' col="red", pch=19)
#' legend("topright", legend=c("CTE with Sexualisation TRN"), 
#' pch=19, col=c(red"))
#' plotR(parameters=fitSTRN$par, main="Sexualisation TRN")
#' # Initial_STRN <- resultNest_4p$par
#' Initial_STRN <- structure(c(4230.10750319997, 510.543319171189, 1015.78663983953,
#' 118.189709917707), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' males <- c(7, 0, 0, 0, 0, 5, 6, 3, 5, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' sexed <- rep(10, length(males))
#' fitSTRN <- STRN(Initial_STRN, EmbryoGrowthTRN=resultNest_4p, tsd=Med_Cc, 
#' Sexed=sexed, Males=males, 
#' Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean")
#' CTE <- info.nests(NestsResult=resultNest_4p, 
#' SexualisationTRN=fitSTRN$par, out="summary")
#' plot(Med_Cc, xlim=c(25, 35))
#' plot_add(x=CTE$TSP.MassWeighted.STRNWeighted.temperature.mean, y=males/sexed, 
#' col="red", pch=19)
#' legend("topright", legend=c("CTE with Sexualisation TRN"), 
#' pch=19, col=c("red"))
#' plotR(parameters=fitSTRN$par, main="Sexualisation TRN")
#' }
#' @export

STRN <- function(Initial_STRN=NULL, 
EmbryoGrowthTRN=stop("Embryo Growth Thermal Reaction Norm must be provided"), 
tsd=stop("A result from the function tsd() must be provided"),
Sexed=NULL, Males=NULL, Females=NULL, Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean", ...)

{

  if (is.null(Initial_STRN)) {pSTRN=EmbryoGrowthTRN$par} else {pSTRN=Initial_STRN}

if (is.null(Sexed)) {Sexed <- Males+Females}
if (is.null(Males)) {Males <- Sexed-Females}
if (is.null(Females)) {Females <- Sexed-Males}

  if (identical(numeric(0), Sexed+Males+Females)) {
  stop("Error in Males, Females or Sexed data")
}

  repeat {
    
    L <- list(...)
    L <- modifyList(list(hessian=FALSE, method="BFGS", control=list(trace=1), 
                         par=pSTRN, fn=.STRN_fit, EmbryoGrowthTRN=EmbryoGrowthTRN, 
                         tsd=tsd, Sexed=Sexed, Males=Males, Temperatures=Temperatures), L)
    
    resultfitSTRN <- do.call(optim, L)

    if (resultfitSTRN$convergence==0) break
    
    pSTRN<-resultfitSTRN$par
    print("Convergence is not achieved. Optimization continues !")
    print(dput(pSTRN))
  }
  
  return(invisible(resultfitSTRN))
}
