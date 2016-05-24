modifyNcdfSetMissing <- function(
##title<< Set missing value attribute to a netCDF file
  con ##<<  file connection to modify
  , var ##<<  variable name (or index) of the variable to modify
  , value = -9999 ##<< value of the missing value attribute
)
  ##description<< This function sets the missing_value and the _Fill_value of a
  ##              variable in a netCDF file to a given value (-9999 by default).
{
  if (var == 'all') {
    vars <- infoNcdfVars(con)$name
  } else {
    vars = var
  }
  for (varT in vars)
    modifyNcdfDefAtts(con, varT, atts = list(missing_value = -9999,
                                   `_FillValue` =  -9999))
}
