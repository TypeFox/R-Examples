## sanitychecks.R

###############################################################
## sanity check definitions

row_sense_is_feasible <- function( x )
  all( x %in% available_row_sense() )

