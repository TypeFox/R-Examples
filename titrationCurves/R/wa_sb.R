#' Titration Curve for a Weak Acid
#' 
#' This function calculates and plots the titration curve for a 
#' monoprotic weak acid analyte using a monoprotic strong base as 
#' the titrant. The calculation uses a single master equation
#' that finds the volume of titrant needed to achieve a fixed pH, 
#' as outlined in R. de Levie's \emph{Principles of Quantitative 
#' Chemical Analysis} (McGraw-Hill, 1997).
#' 
#' @param conc.acid Molar concentration of the weak acid analyte;
#' defaults to 0.10 M.
#' 
#' @param  conc.base Molar concentration of the strong base titrant;
#' defaults to 0.10 M.
#' 
#' @param pka The pKa value for the weak acid analyte; defaults to a
#' pKa of 5.
#' 
#' @param pkw The pKw (or pKs) value for the solvent; defaults to water
#' as a solvent with a pKw of 14.
#' 
#' @param vol.acid Initial volume, in mL, of the solution that 
#' contains the weak acid analyte; defaults to 50.00 mL.
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
#' ex3 = wa_sb(eqpt = TRUE)
#' head(ex3)
#' 
#' ### Overlay titration curves using different pKa values
#' wa_sb(pka = 5, eqpt = TRUE)
#' wa_sb(pka = 7, overlay = TRUE)
#' wa_sb(pka = 9, overlay = TRUE)
#' 
#' ### Overlay titration curve for strong acid and weak acid
#' sa_sb(eqpt = TRUE)
#' wa_sb(overlay = TRUE)

wa_sb = function(conc.acid = 0.1, conc.base = 0.1, pka = 5, pkw = 14, 
                 vol.acid = 50, plot = TRUE,
                 eqpt = FALSE, overlay = FALSE, ...) {
  veq = conc.acid * vol.acid/conc.base
  ka = 10^-pka
  kw = 10^-pkw
  ph = seq(1, pkw, 0.01)
  h = 10^-ph
  oh = kw/h
  delta = h - oh
  alpha = ka/(ka + h)
  volume = vol.acid * (conc.acid * alpha - delta)/(conc.base + delta)
  df = data.frame(volume, ph)
  df = df[df$volume > 0 & df$volume < 2 * veq, ]
  rownames(df) = 1:nrow(df)
  if (plot == TRUE) {
  if (overlay == FALSE) {
    plot(df$volume, df$ph, type = "l", lwd = 2, xlim =
           c(0, 1.5 * veq), ylim = c(0, pkw),
            xlab = "volume of strong base (mL)", ylab = "pH",
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