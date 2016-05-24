
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2016-01-13 15:46:23 +0100 (Wed, 13 Jan 2016) $
# $Rev: 344 $
# ----------------------------------------------------------------

##########################################################################
##-----------------------GetParameterShortName--------------------------##
##########################################################################

GetParameterShortName <- function(longname,
                                  model.par.names) {
  ## Get the short name of specified meteorological parameter
  ## from user.input or get default name to read in NetCDF files properly.
  ##
  ## Args:
  ##   longname: CF-convention name of meteorological parameter to get the
  ##             short name from.
  ##   model.par.names: Named vector of shortnames. Names are
  ##                    CF-convention parameters.
  ##
  ## Returns:
  ##    Short parameter name for CF conventions parameter name.
  ##    This short name should be exactly the variable name
  ##    in the NetCDF file. "Name" attribute remains long name.
  ##
  ## History:
  ##   2010-10-22 | Original code (thm)
  ##   2010-11-03 | user.input adaptation (thm)
  ##   2011-07-21 | redocumentation (geh)

  ## If there is no parameter name definition in InitModelDictionary.R
  ## for a given model, get the default short name. Else extract specified
  ## short name from InitModelDictionary.R.
  if (is.null(model.par.names)) {
    shortname <- GetParameterDefaultShortName(longname)
  } else {
    shortname <- model.par.names[longname]
    names(shortname) <- longname
  }

  ## error if not short name is specified
  if (is.null(shortname))
    stop("NO MATCHING SHORTNAME COULD BE FOUND")

  return(shortname)
}
##########################################################################


##########################################################################
##-------------------GetParameterDefaultShortName-----------------------##
##########################################################################

GetParameterDefaultShortName <- function(longname) {
  ## Converts given parameterstring to short parameter string
  ## eg 'air_temperature' -> 'tas'.
  ## When dealing with new parameters, this is the place to add the default
  ## behaviour.
  ##
  ## Args:
  ##   longname: Vector of CF-convention parameter name
  ##
  ## Returns:
  ##   Default parameter shortname for NetCDF input.

  ## initialize shortname with "NA"s in order to check whether all
  ## parameters have been asigned a short name
  shortname <- rep(NA, length(longname))
  names(shortname) <- longname

  ## get shortnames
  ## get shortname for air_temperature
  is.parameter <- (longname == "air_temperature")
  shortname[is.parameter] <- "tas"

  ## get shortname for precipitation_amount
  is.parameter <- (longname == "precipitation_amount")
  shortname[is.parameter] <- "pr"

  ## get shortname for relative_humidity
  is.parameter <- (longname == "relative_humidity")
  shortname[is.parameter] <- "hurs"

  ## get shortname for global_radiation
  is.parameter <- (longname == "global_radiation")
  shortname[is.parameter] <- "rsds"

  ## get shortname for wind_speed
  is.parameter <- (longname == "wind_speed")
  shortname[is.parameter] <- "wss"

  ## get shortname for eastward_wind
  is.parameter <- (longname == "eastward_wind")
  shortname[is.parameter] <- "ua"

  ## get shortname for northward_wind
  is.parameter <- (longname == "northward_wind")
  shortname[is.parameter] <- "va"

  ## get shortname for "air_pressure_at_sea_level"
  is.parameter <- (longname == "air_pressure_at_sea_level")
  shortname[is.parameter] <- "psl"

  ## get shortname for "specific humidity"
  is.parameter <- (longname == "specific_humidity")
  shortname[is.parameter] <- "hus"

  ## get shortname for sensible heat flux
  is.parameter <- (longname == "surface_upward_sensible_heat_flux")
  shortname[is.parameter] <- "hfss"

  ## get shortname for minimum air temperature
    is.parameter <- (longname == "air_temperature_minimum")
  shortname[is.parameter] <- "tasmin"

  ## get shortname for maximum air temperature
    is.parameter <- (longname == "air_temperature_maximum")
  shortname[is.parameter] <- "tasmax"

  ## get shortname for surface temperature
    is.parameter <- (longname == "surface_temperature")
  shortname[is.parameter] <- "ts"

  ## some indices...

  ## positive_degree_days
  is.parameter <- (longname == "positive_degree_days")
  shortname[is.parameter] <- "pdd"

  ## snow_accumulation
  is.parameter <- (longname == "snow_accumulation")
  shortname[is.parameter] <- "snow_acc"


  ## if there is at least one "NA"-shortname left (i.e. no longame match for
  ## at least one parameter
  if (any(is.na(shortname)))
    stop("    ERROR: AT LEAST ONE SHORTNAME IS UNDEFINED BECAUSE LONGNAME ",
         "COULD NOT BE FOUND IN \"GetParameterDefaultShortName\". PLEASE ",
         "EXTEND THIS FUNCTION WITH THE NEW PARAMETER.")

  return(shortname)
}
##########################################################################
