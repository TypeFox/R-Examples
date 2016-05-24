#' Redox Titration Curve
#' 
#' This function calculates and plots the titration curve for a 
#' reducing agent analyte using an oxidizing agent as the titrant. 
#' The calculation uses a single master equation that finds the volume 
#' of titrant needed to achieve a fixed potential, as outlined in 
#' R. de Levie's \emph{Principles of Quantitative Chemical Analysis} 
#' (McGraw-Hill, 1997). 
#' 
#' @param conc.analyte Molar concentration of the analyte; defaults 
#' to 0.010 M.
#' 
#' @param vol.analyte Initial volume, in mL, of the solution 
#' containing the analyte; defaults to 25.00 mL.
#' 
#' @param pot.analyte Standard state or formal potential for the 
#'analyte's half-reaction in V; defaults to 0.77 V.
#' 
#' @param elec.analyte The number, n, of electrons lost by the analyte 
#' in its oxidation half-reaction; defaults to 1.
#' 
#' @param conc.titrant Molar concentration of the titrant; defaults 
#' to 0.010 M.
#' 
#' @param pot.titrant Standard state or formal potential for the 
#' titrant's half-reaction in V; defaults to 1.7 V.
#' 
#' @param elec.titrant The number, n, of electrons gained by the 
#' analyte in its reduction half-reaction; defaults to 1.
#' 
#' @param plot Logical; if TRUE, plots the titration curve.
#' 
#' @param eqpt Logical; if TRUE, draws a vertical line at the titration
#' curve's equivalence point.
#' 
#' @param overlay Logical; if TRUE, adds the current titration curve
#' to the existing titration curve.
#' 
#' @param \dots Additional arguments to pass to \code{plot()} function.
#' 
#' @return A two-column data frame that contains the volume of titrant
#' in the first column and the solution's potential in the second column. 
#' Also produces a plot of the titration curve with options to display 
#' the equivalence point and to overlay titration curves.
#' 
#' @author David T. Harvey, DePauw University. \email{harvey@@depauw.edu}
#' 
#' @export
#' 
#' @importFrom graphics plot lines 
#' 
#' @examples
#' ### Simple titration curve with equivalence point
#' ex12 = redox_titration(eqpt = TRUE)
#' head(ex12)
#' 
#' ### Overlay titration curves using different potentials for tirant
#' redox_titration(pot.titrant = 1.7, eqpt = TRUE)
#' redox_titration(pot.titrant = 1.5, overlay = TRUE)
#' redox_titration(pot.titrant = 1.3, overlay = TRUE)

redox_titration = function(conc.analyte = 0.01, vol.analyte = 25, 
                         pot.analyte = 0.77, elec.analyte = 1, 
                         conc.titrant = 0.01, pot.titrant = 1.7, 
                         elec.titrant = 1, plot = TRUE, eqpt = FALSE,
                         overlay = FALSE, ...) {
  veq = (elec.analyte * conc.analyte * vol.analyte)/
    (elec.titrant * conc.titrant)
  k.analyte = 10^-(pot.analyte/0.05916)
  k.titrant = 10^-(pot.titrant/0.05916)
  potential = seq(-3, 3, 0.01)
  h = 10^-(potential/0.05916)
  alpha.analyte = (k.analyte^elec.analyte)/
    (h^elec.analyte + k.analyte^elec.analyte)
  alpha.titrant = (h^elec.titrant)/
    (h^elec.titrant + k.titrant^elec.titrant)
  volume = vol.analyte * 
    (elec.analyte * conc.analyte * alpha.analyte)/
    (elec.titrant * conc.titrant * alpha.titrant)
  df = data.frame(volume, potential)
  df = df[df$volume > 0.1 & df$volume < 2 * veq, ]
  rownames(df) = 1:nrow(df)
  if (plot == TRUE) {
  if (overlay == FALSE) {
    plot(df$volume, df$potential, type = "l", lwd = 2,
         xlim = c(0, 1.5 * veq), 
         ylim = c(pot.analyte - 0.5, 1.2 * pot.titrant), 
         xlab = "volume of titrant (ml)", ylab = "potential (V)",
         xaxs = "i", ...)
  } else {
    lines(df$volume, df$potential, type = "l", lwd = 2, ...)
  }
  if (eqpt == TRUE) {
    x = c(veq, veq)
    y = c(-3, 3)
    lines(x, y, type = "l", lty = 2, col = "red")
  }
  }
  invisible(df)
}
