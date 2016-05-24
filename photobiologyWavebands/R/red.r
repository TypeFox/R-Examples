#' Definition of red waveband
#'
#' Red radiation according to "ISO" (610-760 nm) or as commonly defined in plant
#' photobiology, "Smith10" (655-665 nm), "Smith20" (650-670 nm), "Inada"
#' (600-700 nm), "Warrington" (625-675 nm), and "Sellaro" (620-680 nm). No
#' weighting applied.
#'
#' @param std a character string, defaults to "ISO".
#'
#' @return a waveband object wavelength defining a wavelength range.
#'
#' @export
#'
#' @references Aphalo, P. J., Albert, A., Björn, L. O., McLeod, A. R., Robson,
#' T. M., Rosenqvist, E. (Eds.). (2012). Beyond the Visible: A handbook of best
#' practice in plant UV photobiology (1st ed., p. xxx + 174). Helsinki:
#' University of Helsinki, Department of Biosciences, Division of Plant Biology.
#' ISBN 978-952-10-8363-1 (PDF), 978-952-10-8362-4 (paperback). Open access PDF
#' download available at http://hdl.handle.net/10138/37558
#'
#' ISO (2007) Space environment (natural and artificial) - Process for
#' determining solar irradiances. ISO Standard 21348. ISO, Geneva.
#'
#' Murakami, K., Aiga I. (1994) Red/Far-red photon flux ratio used as an index
#' number for morphological control of plant growth under artificial lighting
#' conditions. Proc. Int. Symp. Artificial Lighting, Acta Horticulturae, 418,
#' ISHS 1997.
#'
#' Sellaro, R., Crepy, M., Trupkin, S. A., Karayekov, E., Buchovsky, A. S.,
#' Rossi, C., & Casal, J. J. (2010). Cryptochrome as a sensor of the blue/green
#' ratio of natural radiation in Arabidopsis. Plant physiology, 154(1), 401-409.
#' doi:10.1104/pp.110.160820
#'
#' Smith, H. (1982) Light quality, photoperception and plant strategy. Annual
#' Review of Plant Physiology, 33:481-518.
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' Red()
#' Red("ISO")
#' Red("Smith")
#' Red("Sellaro")
#'
#' @family unweighted wavebands
#'
Red <- function(std="ISO"){
  label="Red"
  if (std=="Smith") {
    warning("The definition of 'Smith' defaults to 'Smith10', to restore old behaviour use 'Smith20'.")
    std <- "Smith10"
  }
  if (std=="ISO") {
    return(new_waveband(610, 760, wb.name=paste("Red", std, sep="."), wb.label=label))
  } else if (std=="Smith20"){
    return(new_waveband(650, 670, wb.name=paste("Red", std, sep="."), wb.label="R"))
  }else if (std=="Smith10"){
    return(new_waveband(655, 665, wb.name=paste("Red", std, sep="."), wb.label="R"))
  }else if (std=="Inada"){
    return(new_waveband(600, 700, wb.name=paste("Red", std, sep="."), wb.label="R"))
  }else if (std=="Warrington"){
    return(new_waveband(625, 675, wb.name=paste("Red", std, sep="."), wb.label="R"))
  } else if (std=="Sellaro"){
    return(new_waveband(620, 680, wb.name=paste("Red", std, sep="."), wb.label=label))
  } else {
    warning("'std' argument value not implemented.")
    return(NA)
  }
}
