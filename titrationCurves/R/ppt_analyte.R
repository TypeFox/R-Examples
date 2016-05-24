#' Precipitation Titration Curve: Monitoring pAnalyte
#' 
#' This function calculates and plots the precipitation titration curve 
#' for an analyte and a titrant that form a precipitate with a 1:1
#' stoichiometry. The calculation uses a single master equation
#' that finds the volume of titrant needed to achieve a fixed
#' concentration of the analyte, expressed as pAnalyte, as outlined in 
#' R. de Levie's \emph{Principles of Quantitative Chemical Analysis} 
#' (McGraw-Hill, 1997).
#' 
#' @param conc.analyte Molar concentration of the analyte; defaults to
#' 0.025 M.
#' 
#' @param conc.titrant Molar concentration of the titrant; defaults to
#' 0.050 M. 
#' 
#' @param vol.analyte The initial volume, in mL, of the solution 
#' containing the analyte; defaults to 50.00 mL. 
#' 
#' @param pksp The pKsp value for the precipitate; defaults to 16.08,
#' which is the pKsp for AgI.
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
#' in the first column and the solution's pAnalyte in the second column. 
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
#' ex13 = ppt_analyte(eqpt = TRUE)
#' head(ex13)
#' 
#' ### Overlay titration curves using different pKsp values 
#' ppt_analyte(pksp = 16, eqpt = TRUE)
#' ppt_analyte(pksp = 14, overlay = TRUE)
#' ppt_analyte(pksp = 12, overlay = TRUE)

ppt_analyte = function(conc.analyte = 0.025, conc.titrant = 0.05, 
                    vol.analyte = 50, pksp = 16.08, plot = TRUE, 
                    eqpt = FALSE, overlay = FALSE, ...) {
  veq = conc.analyte * vol.analyte/conc.titrant
  ksp = 10^-pksp
  p.analyte = seq(1, pksp, 0.01)
  analyte = (10^-p.analyte)
  titrant = ksp/analyte
  volume = vol.analyte * (conc.analyte - analyte + titrant)/
    (conc.titrant + analyte - titrant)
  df = data.frame(volume, p.analyte)
  df = df[df$volume > 0 & df$volume < 2 * veq, ]
  rownames(df) = 1:nrow(df)
  if (plot == TRUE) {
  if (overlay == FALSE) {
    plot(df$volume, df$p.analyte, type = "l", lwd = 2, 
         xlim = c(0, 1.5 * veq), ylim = c(0, pksp), 
         xlab = "volume of titrant (ml)", ylab = "pAnalyte",
         xaxs = "i", yaxs = "i", ...)
  } else {
    lines(df$volume, df$p.analyte, type = "l", lwd = 2, 
          xlim = c(0, 1.5 * veq), ylim = c(0, pksp), 
          xlab = "volume of titrant (ml)", ylab = "pAnalyte", ...)
  }
  if (eqpt == TRUE) {
    x = c(veq, veq)
    y = c(-1, pksp + 1)
    lines(x, y, type = "l", lty = 2, col = "red")
  }
  }
  invisible(df)
}