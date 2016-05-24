#' Joint plot for in and out-migration fields
#'
#' This migration field diagram makes easy to visualize both direction of migration. E.g. points above the diagonal "are outward redistributors, while those below that line are inward redistributors."
#' @param m migration matrix
#' @param method measurement of in and out-migration
#' @param title plot title
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @export
#' @references \itemize{
#'   \item Source code was adopted from Michael Ward and Kristian Skrede Gleditsch (2008) \emph{Spatial Regression Models}. Thousand Oaks, CA: Sage. \url{http://privatewww.essex.ac.uk/~ksg/code/srm_enhanced_code_v5.R} with the permission of the authors.
#'   \item Case study and use case: Andrei Rogers and Stuart Sweeney (1998) Measuring the Spatial Focus of Migration Patterns. \emph{The Professional Geographer} \bold{50}, 232--242
#' }
#' @importFrom calibrate textxy
#' @examples \dontrun{
#' data(migration.world)
#' par(mfrow = c(2, 1))
#' migration.field.diagram(migration.world)
#' migration.field.diagram(migration.world, method = 'acv')
#' }
migration.field.diagram <- function(m, method = c('gini', 'acv'), title = 'Migration field diagram', xlab = 'Out-migration', ylab = 'In-migration') {

    ## compute indices
    method <- match.arg(method)
    switch(method,
           'gini' = {
               .in  <- migration.gini.in(m)
               .out <- migration.gini.out(m)
           },
           'acv'  = {
               .in  <- migration.cv.in(m)
               .out <- migration.cv.out(m)
           })

    ## z-scale
    .in  <- as.numeric(scale(.in))
    .out <- as.numeric(scale(.out))

    ## plot main graph
    plot(.out, .in, xlim = c(-3,3), ylim = c(-3,3), pch = 20, las = 1, bty = 'n', xlab = xlab, ylab = ylab, main = title)

    ## 1+2 standard deviation boxes
    lines(c(-2, -2, 2, 2, -2), c(-2, 2, 2, -2, -2))
    lines(c(-1, -1, 1, 1, -1), c(-1, 1, 1, -1, -1))
    lines(c(-2, 2), c( 0, 0))
    lines(c( 0, 0), c(-2, 2))

    ## annotations
    text(-2,  3, "('Pure outward')")
    text( 2,  3, "('Intensive')")
    text(-2, -3, "('Extensive')")
    text( 2, -3, "('Pure inward')")
    polygon(x = c(-1,  0,  0, -1), y = c(-1, -1,  0,  0), col = 'Light Blue 3')
    polygon(x = c( 0,  1,  1,  0), y = c( 0,  0,  1,  1), col = 'Light Blue 3')
    polygon(x = c( 0, -1, -1,  0), y = c( 0,  0,  1,  1), col = 'Light Blue 3')
    polygon(x = c( 0,  1,  1,  0), y = c( 0,  0, -1, -1), col = 'Light Blue 3')

    ## model
    fit   <- lm(.in ~ .out)

    ## confidence interval
    xgrid <- seq(-3, 3, length.out = 20)
    pred  <- predict(fit, data.frame(.out = xgrid), interval = 'confidence')
    polygon(x = c(xgrid, rev(xgrid)), y=c(pred[,3], rev(pred[,2])), col = 'Light Blue 3', border = TRUE)

    ## points
    points(.out, .in, pch = 20)
    textxy(.out, .in, labs = rownames(m))

    ## regression parameters
    rmse  <- round(sqrt(mean(resid(fit) ^ 2)), 2)
    coefs <- coef(fit)
    b0    <- round(coefs[1], 2)
    b1    <- round(coefs[2], 2)
    r2    <- round(summary(fit)$r.squared, 2)
    mtext(bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ r^2 == .(r2) * "," ~~ RMSE == .(rmse)), side = 3, line = 0)

    ## densities
    sldensity <- density(.in)
    lines(sldensity$y + 2, sldensity$x, lty = 2, col = 'Light Blue 3')
    ddensity <- density(.out)
    lines(ddensity$x, ddensity$y + 2, lty = 2, col = 'Light Blue 3', xlim = c(-2, 2))

    ## regression line
    lines(xgrid, pred[,1], type = 'l', lty = 2, col = 'Blue 4', lwd = 2)

    ## rugs
    rug(jitter(.out, factor = 2), col = 'Light Blue 3')
    rug(.in, side = 2, col = 'Light Blue 3')

}
