wDC <- function(DCf = 100, rw = 200, a = 0.75, b = 0.75) {
  #############################################################################
  # Description:
  #   Computes the over wintering Drought Code (DC) value.
  #   All variables names are laid out in the same manner as Lawson & Armitage
  #   (2008).
  
  #   Lawson B.D. and Armitage O.B. 2008. Weather Guide for the Canadian Forest 
  #   Fire Danger Rating System. Natural Resources Canada, Canadian Forest 
  #   Service, Northern Forestry Centre, Edmonton, Alberta. 84 p. 
  #   http://cfs.nrcan.gc.ca/pubwarehouse/pdfs/29152.pdf
  #
  # Args:
  #     DCf: Final Fall DC Value
  #      rw: winter precipitation (mm)
  #       a: user-selected value accounting for carry-over fraction
  #       b: user-selected value accounting for wetting efficiency fraction
  # Returns:
  #     DCs: Overwintered Drought Code (Spring startup DC value)
  #
  #############################################################################
  #Eq. 3 - Final fall moisture equivalent of the DC
  Qf <- 800 * exp(-DCf / 400)
  #Eq. 2 - Starting spring moisture equivalent of the DC
  Qs <- a * Qf + b * (3.94 * rw)
  #Eq. 4 - Spring start-up value for the DC
  DCs <- 400 * log(800 / Qs)
  #Constrain DC
  DCs <- ifelse(DCs < 15, 15, DCs)
  return(DCs)
}

