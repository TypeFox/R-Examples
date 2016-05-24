#' Precipitation Titration Curve: Mixture of Analytes
#' 
#' This function calculates and plots the precipitation titration curve 
#' for a mixture of two analytes using a titrant that form precipitates 
#' with 1:1 stoichiometries. The calculation uses a single master equation
#' that finds the volume of titrant needed to achieve a fixed 
#' concentration of titrant, expressed as pTitrant, as outlined in 
#' R. de Levie's \emph{Principles of Quantitative Chemical Analysis} 
#' (McGraw-Hill, 1997).
#' 
#' @param conc.analyte1 Molar concentration of the first analyte; 
#' defaults to 0.050 M.
#' 
#' @param conc.analyte2 Molar concentration of the second analyte; 
#' defaults to 0.050 M.
#' 
#' @param conc.titrant Molar concentration of the titrant; 
#' defaults to 0.050 M. 
#' 
#' @param vol.analyte The initial olume, in mL, of the solution 
#' containing the analyte; defaults to 25.00 mL. 
#' 
#' @param pksp1 The pKsp value for the first analyte's precipitate; 
#' defaults to 16.08, which is the pKsp for AgI.
#' 
#' @param pksp2 The pKsp value for the second analyte's precipitate; 
#' defaults to 11.97, which is the pKsp for AgSCN.
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
#' in the first column and the solution's pTitrant in the second column. 
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
#' ### Simple titration curve with equivalence points
#' ex15 = ppt_mixture(eqpt = TRUE)
#' head(ex15)
#' 
#' ### Overlay mixture titration curves using different pKsp values 
#' ppt_mixture(pksp1 = 16, pksp2 = 12, eqpt = TRUE)
#' ppt_mixture(pksp1 = 14, pksp2 = 10, overlay = TRUE)

ppt_mixture = function(conc.analyte1 = 0.05, conc.analyte2 = 0.05, 
                       vol.analyte = 25, conc.titrant = 0.05, 
                       pksp1 = 16.08, pksp2 = 11.97, plot = TRUE,
                       eqpt = FALSE, overlay = FALSE, ...) {
  veq1 = conc.analyte1 * vol.analyte/conc.titrant
  veq2 = conc.analyte2 * vol.analyte/conc.titrant
  veq = veq1 + veq2
  ksp1 = 10^-pksp1
  ksp2 = 10^-pksp2
  y.lim = max(c(pksp1, pksp2))
  p.titrant = seq(1, y.lim, 0.01)
  titrant = 10^-p.titrant
  analyte1 = ksp1/titrant
  analyte2 = ksp2/titrant
  volume1 = vol.analyte * (conc.analyte1 - analyte1 + titrant)/
    (conc.titrant + analyte1 - titrant)
  volume2 = vol.analyte * 
    (conc.analyte1 + conc.analyte2 - analyte1 - analyte2 + titrant)/
    (conc.titrant + analyte1 + analyte2 - titrant)
  df = data.frame(volume1, volume2)
  volume = apply(df, 1, max)
  df = data.frame(volume, p.titrant)
  df = df[df$volume > 0 & df$volume < 2 * veq, ]
  rownames(df) = 1:nrow(df)
  if (plot == TRUE) {
  if (overlay == FALSE) {
    plot(df$volume, df$p.titrant, type = "l", lwd = 2, 
         xlim = c(0, 1.5 * veq), ylim = c(0, y.lim), 
         xlab = "volume of titrant (ml)", ylab = "pTitrant",
         xaxs = "i", yaxs = "i", ...)
  } else {
    lines(df$volume, df$p.titrant, type = "l", lwd = 2, ...)
  }
  if (eqpt == TRUE) {
    x1 = c(veq1, veq1)
    x2 = c(veq, veq)
    y = c(-1, y.lim + 1)
    lines(x1, y, type = "l", lty = 2, col = "red")
    lines(x2, y, type = "l", lty = 2, col = "red")
  }
  }
  invisible(df)
}
