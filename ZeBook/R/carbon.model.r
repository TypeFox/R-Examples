################################################################################
# "Working with dynamic models for agriculture"
# R script for practical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2012-04-23
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
# the Soil Carbon model function - Update Z
################################################################################
#' @title The CarbonSoil model - calculate change in soil carbon for one year
#' @description Simple dynamic model of soil carbon content, with a time step of one year. The equations that describe the dynamics of this system are adapted from the Henin-Dupuis model described in Jones et al. (2004). The soil carbon content is represented by a single state variable: the mass of carbon per unit land area in the top 20 cm of soil in a given year (Z, kg.ha-1). It is assumed that soil C is known in some year, taken as the initial year. 
#' The yearly change in soil C is the difference between input from crop biomass and loss.
#' @param Zy : Soil carbon content for year 
#' @param R : the fraction of soil carbon content lost per year
#' @param b : the fraction of yearly crop biomass production left in the soil
#' @param Uy : the amount of C in crop biomass production in the given year
#' @return Soil carbon content for year+1. 
#' @export
# Arguments: Zy=Z(year), 2 parameters, Uy=U(year).
carbonsoil.update<-function(Zy,R,b,Uy)
{
	# Calculate the rates of state variables (dZ)
	dZ = -R*Zy+b*Uy
	# Update the state variables Z
	Zy1=Zy+dZ
  # Return Zy1=Z(year+1)
	return(Zy1)
}

################################################################################
#' @title The CarbonSoil model - calculate daily values over designated time period
#' @description Simple dynamic model of soil carbon content, with a time step of one year. The equations that describe the dynamics of this system are adapted from the Henin-Dupuis model described in Jones et al. (2004). The soil carbon content is represented by a single state variable: the mass of carbon per unit land area in the top 20 cm of soil in a given year (Z, kg.ha-1). It is assumed that soil C is known in some year, taken as the initial year. 
#' The yearly change in soil C is the difference between input from crop biomass and loss.
#' @param R : the fraction of soil carbon content lost per year
#' @param b : the fraction of yearly crop biomass production left in the soil
#' @param U : the amount of C in crop biomass production (constant or time series)
#' @param Z1 : initial soil carbon content
#' @param duration :  duration of simulation (year))
#' @return Soil carbon years. 
#' @export
carbonsoil.model<-function(R,b,U,Z1,duration)
{
	# U : if U is a constant (not a dataframe), convert into dataframe.
	if (!is.data.frame(U)) {U=data.frame("year"=1:duration, "U"=rep(U,duration))}
	# Initialize variables
	# 1 state variable, as 1 vectors initialized to NA
	# Z : mass of carbon per area (top 20 cm of soil) (kg ha?1)
	Z = rep(NA,duration)
		# Initialize state variable when starting simulation on year 0
	Z[1] = Z1
		# Integration loop
	for (year in 1:(duration-1))
		{
		# using the update function.
		Z[year+1] = carbonsoil.update(Z[year],R,b,U$U[year])
		}
		# End loop
  return(data.frame(year=1:duration,Z=Z[1:duration]))
}

# End of file
