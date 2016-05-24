#' Titration Curve for a Strong Base
#' 
#' This function calculates and plots the titration curve for a 
#' monoprotic strong base analyte using a monoprotic strong acid as 
#' the titrant. The calculation uses a single master equation
#' that finds the volume of titrant needed to achieve a fixed pH, 
#' as outlined in R. de Levie's \emph{Principles of Quantitative 
#' Chemical Analysis} (McGraw-Hill, 1997).
#'  
#' @param conc.base Molar concentration of the strong base analyte;
#' defaults to 0.10 M.
#' 
#' @param  conc.acid Molar concentration of the strong acid titrant;
#' defaults to 0.10 M.
#' 
#' @param pkw The pKw (or pKs) value for the solvent; defaults to water
#' as a solvent with a pKw of 14.
#' 
#' @param vol.base Initial volume, in mL, of the solution that 
#' contains the strong base analyte; defaults to 50.00 mL.
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
#' in the first column and the solution's pH in the second column. Also
#' produces a plot of the titration curve with options to display the
#' equivalence point and to overlay titration curves.
#' 
#' @author David T. Harvey, DePauw University. \email{harvey@@depauw.edu}
#' 
#' @export
#' 
#' @importFrom graphics plot lines 
#' 
#' @examples
#' ### Simple titration curve with equivalence point
#' ex2 = sb_sa(eqpt = TRUE)
#' head(ex2)
#' 
#' ### Overlay titration curves
#' sb_sa(conc.acid = 0.10)
#' sb_sa(conc.acid = 0.15, overlay = TRUE)
#' sb_sa(conc.acid = 0.20, overlay = TRUE)

sb_sa = function(conc.base = 0.1, conc.acid = 0.1, pkw = 14, 
                 vol.base = 50, plot = TRUE,
                 overlay = FALSE, eqpt = FALSE, ...) {
  veq = conc.base * vol.base / conc.acid
  kw = 10 ^ -pkw
  ph = seq(pkw, 1,-0.01)
  h = 10 ^ -ph
  oh = kw / h
  delta = h - oh
  volume = vol.base * (conc.base + delta) / (conc.acid - delta)
  df = data.frame(volume, ph)
  df = df[df$volume > 0 & df$volume <2 * veq, ]
  rownames(df) = 1:nrow(df)
  if (plot == TRUE) {
  if (overlay == FALSE) {
    plot(df$volume, df$ph, type = "l", lwd = 2, xlim =
         c(0, 1.5 * veq), ylim = c(0, pkw), 
       xlab = "volume of strong acid (mL)", ylab = "pH", 
       xaxs = "i", yaxs = "i", ...)
} else {
  lines(df$volume, df$ph, type = "l", lwd = 2, ...)
}
if (eqpt == TRUE) {
  x = c(veq, veq)
  y = c(0, pkw + 1)
  lines(x, y, type = "l", lty = 2, col = "red")
}
  }
invisible(df)
}