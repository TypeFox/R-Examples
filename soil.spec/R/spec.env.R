#' Enviroment for standard parameters
#'
#' @author Andrew Sila \email{asila@cgiar.org} and Tomislav Hengl \email{tom.hengl@wur.nl}
#' Note: Generated onLoad

spec.opts <- new.env(hash=TRUE)

spec.env <- function(MIR = c(390, 7500), NIRS = c(3900, 12500), NIRP = c(4000,10000), VISNIR1 = c(420, 960), VISNIR2 = c(1020, 1770), VISNIR3 = c(1830, 2480), icraf.htsxt = c(3578, 7497.964, 599.76), icraf.alpha = c(2542, 3998.12872, 399.387991), icraf.mpa = c(2307, 12493.2, 3598.69), CO2.band = c(2350.8,2379.8), signif.digit=5, attributes = c("ORCCNS", "PHIHOX", "ALUM3S", "ECAM3S", "EXKM3S", "EMGM3S", "ENAM3S", "EXB", "NITCNS", "SNDLDF"), mdnames = c("MID", "Instrument_name", "Instrument_URL", "Laboratory_name", "Laboratory_contact", "Laboratory_URL", "Material_class", "Wavenumber_conversion", "Wavenlength_unit", "Location_error"), show.env = FALSE){
   pl.lst <- list(
     MIR = MIR,
     NIRS = NIRS,
     NIRP = NIRP,
     VISNIR1 = VISNIR1,
     VISNIR2 = VISNIR2,
     VISNIR3 = VISNIR3,
     icraf.htsxt = icraf.htsxt,
     icraf.alpha = icraf.alpha,
     icraf.mpa = icraf.mpa,
     CO2.band = CO2.band,
     signif.digit = signif.digit,
     attributes = attributes,
     mdnames = mdnames
   )
   x <- lapply(names(pl.lst), function(x){ assign(x, pl.lst[[x]], envir=spec.opts) })
   if(show.env){
     return(pl.lst)
   }
}

spec.env()