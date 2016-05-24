#' Titration Curve for Weak Acid Mixture
#' 
#' This function calculates and plots the titration curve for a 
#' mixture of two monoprotic weak acid analyte using a monoprotic 
#' strong base as the titrant. The calculation uses a single master 
#' equation that finds the volume of titrant needed to achieve a fixed 
#' pH, as outlined in R. de Levie's \emph{Principles of Quantitative 
#' Chemical Analysis} (McGraw-Hill, 1997).
#' 
#' @param conc.acid1 Molar concentration of the first monoprotic 
#' weak acid analyte; defaults to 0.10 M.
#' 
#' @param conc.acid2 Molar concentration of the second monoprotic 
#' weak acid analyte; defaults to 0.10 M.
#' 
#' @param  conc.base Molar concentration of the strong base titrant;
#' defaults to 0.10 M.
#' 
#' @param pka1 The pKa value for the first monoprotic weak acid 
#' analyte; defaults to a pKa of 5.
#' 
#' @param pka2 The pKa value for the second monoprotic weak acid 
#' analyte; defaults to a pKa of 8.
#' 
#' @param pkw The pKw (or pKs) value for the solvent; defaults to water
#' as a solvent with a pKw of 14.
#' 
#' @param vol.acid Initial volume, in mL, of the solution that 
#' contains the weak acid analytes; defaults to 50.00 mL.
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
#' ### Simple titration curve with equivalence points
#' ex9 = wamix_sb(eqpt = TRUE)
#' head(ex9)
#' 
#' ### Overlay titration curves using different pKa values
#' wamix_sb(pka1 = 5, pka2 = 8, eqpt = TRUE)
#' wamix_sb(pka1 = 4, pka2 = 7, overlay = TRUE)
#' wamix_sb(pka1 = 6, pka2 = 9, overlay = TRUE)

wamix_sb = function(conc.acid1 = 0.1, conc.acid2 = 0.1, 
                    conc.base = 0.1, pka1 = 5, pka2 = 8, pkw = 14, 
                    vol.acid = 50, plot = TRUE, eqpt = FALSE, 
                    overlay = FALSE, ...) {
  veq1 = conc.acid1 * vol.acid/conc.base
  veq2 = conc.acid2 * vol.acid/conc.base
  ka1 = 10^-pka1
  ka2 = 10^-pka2
  kw = 10^-pkw
  ph = seq(1, pkw, 0.01)
  h = 10^-ph
  oh = kw/h
  delta = h - oh
  alpha1 = ka1/(ka1 + h)
  alpha2 = ka2/(ka2 + h)
  volume = vol.acid * 
    (conc.acid1 * alpha1 + conc.acid2 * alpha2 - delta)/
    (conc.base + delta)
  df = data.frame(volume, ph)
  df = df[df$volume > 0 & df$volume < 2 * (veq1 + veq2), ]
  rownames(df) = 1:nrow(df)
  if (plot == TRUE) {
  if (overlay == FALSE){
    plot(df$volume, df$ph, type = "l", lwd = 2, 
         xlim = c(0, 1.5 * (veq1 + veq2)), ylim = c(0, pkw), 
         xlab = "volume of strong base (mL)", ylab = "pH",
         xaxs = "i", yaxs = "i", ...)
  }else{
    lines(df$volume, df$ph, type = "l", lwd = 2, ...)
  }
  if (eqpt == TRUE) {
    x1 = c(veq1, veq1)
    x2 = c(veq1 + veq2, veq1 + veq2)
    y = c(0, pkw + 1)
    lines(x1, y, type = "l", lty = 2, col = "red")
    lines(x2, y, type = "l", lty = 2, col = "red")
  }
  }
  invisible(df)
}