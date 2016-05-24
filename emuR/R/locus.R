##' Calculate locus equations for two-dimensional data
##' 
##' The function plots a locus equation and returns associated statistical
##' information.
##' 
##' A locus equation is a straight line regression fitted with lm() in which
##' the F2- values typically at the vowel onset are regressed on those of the
##' target. The slope can be used to give an indication of target-on-onset
##' coarticulatory influences.
##' 
##' The best estimate of the locus frequency is where the locus equation
##' bisects the line target = onset.
##' 
##' @param target a numerical vector typically of F2 values at the vowel target
##' @param onset a numerical vector typically of the same length as target of
##' F2 values at the vowel onset
##' @param labels.vow an optionally character vector for plotting labels at the
##' points (target, onset) of the same length as target
##' @param yxline optionally plot the line target = onset.  Defaults to True.
##' @param plotgraph a logical vector for specifying whether the data should be
##' plotted. Defaults to True.
##' @param axes A logical vector indicating whether the axes should be plotted
##' @param ... graphical options \link{par}
##' @return A list containing regression diagnostics of the function lm() that
##' can be accessed with summary() and the estimated locus frequency in
##' \$locus. A plot of values in the onset x target plane with superimposed
##' locus equation and line onset=target.
##' @author Jonathan Harrington
##' @keywords math
##' @examples
##' 
##' 
##'  # calculate an F2-locus equation for initial [d] 
##' # preceding lax vowels produced by female speaker "68".
##' # the onset is taken at the vowel onset; the
##' # vowel target is taken at the vowel's temporal midpoint.
##' 
##' # identify initial "d" of speaker "68"
##' temp <- vowlax.left == "d" & vowlax.spkr == "68"
##' # get the F2 value at the vowel's temporal midpoint
##' targ <- dcut(vowlax.fdat[temp,2], .5, prop=TRUE)
##' # F2 value at the vowel's acoustic onset.
##' on <- dcut(vowlax.fdat[temp,2], 0, prop=TRUE)
##' 
##' # locus equation plot
##' result <- locus(targ, on, vowlax.l[temp])
##' # statistical diagnostics of the regression line (locus equation)
##' summary(result)
##' # intercept and slope
##' result$coeff
##' # best estimate of the locus frequency, i.e. the
##' # point of bisection of on = TRUEarg with the regression line
##' result$locus
##' 
##' 
##' @export locus
`locus` <- function (target, onset, labels.vow = NULL, 
                     yxline = TRUE, plotgraph = TRUE, 
                     axes=TRUE,  ...) 
{
  # target: vector of target freqs
  # onset: vector of onset freqs
  # labels.vow: optional vowel labels for plotting
  # xlim, ylim: optional range for x and y-axes
  # xlab, ylab: optional label for axes
  # plot: if T, produces a plot of  target x  onset
  # with a superimposed regression line through
  # the scatter with superimposed line target = onset
  # returns: statistical results of the fitting
  # the regression line and the locus frequency in $locus
  regr <- lm(onset ~ target)
  mat <- rbind(cbind(1, -1), cbind(1, -regr$coef[2]))
  vec <- c(0, regr$coef[1])
  regr$locus <- ((solve(mat) %*% vec)[1])
  if (plotgraph) {
    if (is.null(labels.vow)) 
      labels.vow <- rep("x", length(target))
    graphics::plot(target, onset,  type = "n", axes=FALSE,  ...)
    if(axes)
    {
      graphics::axis(side=1)
      graphics::axis(side=2)
    }
    if(is.character(labels.vow))
      graphics::text(target, onset, labels.vow, ...)
    else if(is.numeric(labels.vow))
      graphics::points(target, onset, pch=labels.vow, ...)
    graphics::abline(regr, ...)
    if (yxline) 
      graphics::abline(0, 1, lty = 2)
  }
  
  regr
}
