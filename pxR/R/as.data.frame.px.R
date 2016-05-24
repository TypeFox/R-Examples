#################################################################
# 
# File:         as.data.frame.px.R 
# Purpose:      extracts a df from a px object
#
# Created:      20110801
# Authors:      cjgb, opl, fvf
#
# Modifications: 
#    20120323, cjgb: added error check: variables and codes should have the same length
#    20120402, cjgb: warnings can be either errors or warnings depending on paranoia level
#    20120402, cjgb: adapted to the new px object format (where DATA is already a df)
#
#################################################################

as.data.frame.px <- function( x, ..., use.codes = FALSE, warnings.as.errors = TRUE, direction = "long"){

  # stores the user settings and resets them on exit
#   initial.warning.option <- options("warn")$warn
#   on.exit( options(warn=initial.warning.option) )   
# 
#   if ( warnings.as.errors )
#     options(warn=2)
  
  dat <- x$DATA$value  # already a data frame
  
  ## maybe we need to change values to codes
  if (is.logical(use.codes) && use.codes)
    use.codes <- names(x$CODES)
  
  if (! is.logical(use.codes))
    for( var.name in intersect( use.codes, intersect(colnames(dat), names(x$CODES) ) ) )
        dat[[var.name]] <- mapvalues(dat[[var.name]], 
                                     from = x$VALUES[[var.name]], 
                                     to   = x$CODES[[var.name]])
     
  ## do we need to reshape?
  if (direction == "wide")
    dcast(dat, list(x$HEADING$value, x$STUB$value))
  else
    dat
}

