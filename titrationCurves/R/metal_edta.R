#' Complexation Titration Curve
#' 
#' This function calculates and plots the titration curve for a 
#' metal ion analyte using EDTA as the titrant. The calculation uses 
#' a single master equation that finds the volume of titrant needed to 
#' achieve a fixed concentration of the metal ion, pM, as outlined in 
#' R. de Levie's \emph{Principles of Quantitative Chemical Analysis} 
#' (McGraw-Hill, 1997).
#' 
#' @param conc.metal Molar concentration of the metal ion analyte;
#' defaults to 0.10 M.
#' 
#' @param conc.edta Molar concentration of the EDTA titrant;
#' defaults to 0.10 M.
#' 
#' @param vol.metal Initial volume, in mL, of the solution that contains 
#' the metal ion analyte; defaults to 50.00 mL.
#' 
#' @param ph The pH of the solution, which is used to calculate the
#' fraction of EDTA present in its fully deprotonated form; defaults to
#' a pH of 10.
#' 
#' @param logkf The log of the formation constant, Kf, for the metal-EDTA
#' complex; defaults to 8.79, which is the value for the complex of
#' Mg2+ and EDTA.
#' 
#' @param alpha.metal The fraction of the metal ion analyte that is not
#' complexed by an auxilary complexing agent; defaults to 1, the value
#' when there is no secondary complexing agent present.
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
#' in the first column and the solution's pMetal in the second column. 
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
#' ex11 = metal_edta(eqpt = TRUE)
#' head(ex11)
#' 
#' ### Overlay titration curves using different pH values
#' metal_edta(ph = 12, eqpt = TRUE)
#' metal_edta(ph = 10, overlay = TRUE)
#' metal_edta(ph = 8, overlay = TRUE)

metal_edta = function(conc.metal = 0.1, conc.edta = 0.1, 
                      vol.metal = 50, ph = 10, logkf = 8.79, 
                      alpha.metal = 1, plot = TRUE, eqpt = TRUE, 
                      overlay = FALSE, ...) {
  ka1 = 1
  ka2 = 0.032
  ka3 = 0.01
  ka4 = 0.0022
  ka5 = 6.9e-07
  ka6 = 5.8e-11
  veq = conc.metal * vol.metal/conc.edta
  h = 10^-ph
  alpha.edta = (ka1 * ka2 * ka3 * ka4 * ka5 * ka6)/
    (h^6 + 
       ka1 * h^5 + 
       ka1 * ka2 * h^4 + 
       ka1 * ka2 * ka3 * h^3 + 
       ka1 * ka2 * ka3 * ka4 * h^2 + 
       ka1 * ka2 * ka3 * ka4 * ka5 * h + 
       ka1 * ka2 * ka3 * ka4 * ka5 * ka6
     )
  kf = 10^logkf
  kf.cond = alpha.edta * alpha.metal * kf
  p.metal = seq(0, logkf + 2, 0.01)
  metal = (10^-p.metal)/alpha.metal
  alpha.metal.edta = (kf.cond * metal)/(1 + kf.cond * metal)
  volume = vol.metal * (conc.metal - metal)/
    (alpha.metal.edta * conc.edta + metal)
  df = data.frame(volume, p.metal)
  df = df[df$volume > 0 & df$volume < 2 * veq, ]
  rownames(df) = 1:nrow(df)
  if (plot == TRUE) {
  if (overlay == FALSE) {
    plot(df$volume, df$p.metal, type = "l", lwd = 2, 
         xlim = c(0, 1.5 * veq), ylim = c(0, logkf + 2), 
         xlab = "volume of titrant (ml)", ylab = "pMetal",
         xaxs = "i", yaxs = "i", ...)
  } else {
    lines(df$volume, df$p.metal, type = "l", lwd = 2, ...)
  }
  if (eqpt == TRUE) {
    x = c(veq, veq)
   y = c(-1, logkf + 3)
    lines(x, y, type = "l", lty = 2, col = "red")
  }
  }
  invisible(df)
}