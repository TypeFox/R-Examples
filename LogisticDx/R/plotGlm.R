#' @name plot.glm
#' @rdname plotGlm
#' @aliases plot.glm
#' @method plot glm
#' @export
#'
#' @include dx.R
#' @include genBinom.R
#'
#' @title Plot diagnostics for a binomial \code{glm} model
#'
#' @description Standard diagnostic plots.
#'
#' @param x A regression model with class \code{glm} and
#' \code{x$family$family == "binomial"}.
#' @param y Not used. Present for compatibility with
#' generic \code{plot()} function.
#' @param ... Additional arguments, which can be
#' passed to the plotting functions. See:
#' \cr
#' ?graphics::plot.default
#' \cr
#' ?graphics::symbols
#' \cr
#' ?graphics::par
#' @param toPdf \itemize{
#'  \item If \code{toPdf=TRUE} the output will be directed
#'       to a \code{.pdf} file.
#'  \item If \code{toPdf=FALSE} a new device is opened for each plot.
#' }
#' @param file Filename if writing to \code{.pdf} as above,
#'  e.g. \code{"plots.pdf"}.
#' @param palette Palette of colors to use as the
#'  'fill'/ 'background' colors for the plots.
#' \cr
#' The options are taken from
#' \href{http://colorbrewer2.org/}{color_brewer}.
#' \cr
#' @param usePalette Use the colorscheme in palette above.
#' \itemize{
#'  \item If \code{usePalette=TRUE} (the default), this colorscheme
#'         will be passed to the argument \code{bg} below:
#'         \itemize{
#'          \item \code{graphics::plot.default(bg= )}
#'          \item \code{graphics::symbols(bg= )}
#'          }
#'  \item If \code{usePalette=FALSE}, then the color specified
#'        in \code{bg} below will be used instead.
#' }
#' @param bg The 'fill' or background color(s) to use, if
#' \code{usePalette=FALSE}.
#' \cr
#' This can be a \code{vector} of colors.
#' @param col The 'edge' or 'foreground' color used
#' to outline points in the plot.
#' \cr
#' The default, \code{"white"} is used to make overlapping points
#' easier to see.
#' \cr
#' This is passed as an argument to
#' \itemize{
#'  \item \code{graphics::plot.default(col= )}
#'  \item \code{graphics::symbols(fg= )}
#' }
#' @param alpha Transparency for colors above.
#' \cr
#' Should be in the range \code{0} (transparent) to \code{1} (opaque). See:
#' \cr
#' ?grDevices::adjustcolor
#' @param cex \bold{C}haracter \bold{ex}pansion.
#'  \cr
#'  A multiplier used for size of the plotting symbols/ characters. See:
#'  \cr
#'  ?graphics::par
#' @param pch \bold{P}lotting \bold{ch}aracter.
#'  \cr
#'  The symbol/ character to for the plot.
#'  \cr
#'  The default, \code{pch=21} shows filled circles at each point. See:
#' \cr
#' ?graphics::points
#' @param cex.main \bold{C}haracter \bold{ex}pansion for
#' the plot title and the labels for the axes.
#' @param inches Width of circles for the bubble plot. See
#' \cr
#' ?graphics::symbols
#' @param identify If \code{TRUE} will give option to identify
#' individual points on a number of the plots produced.
#' \cr
#' The number which appears next to the point corresponds
#' to the relevant row as given by \code{\link{dx}}.
#' \cr
#' This may be useful for identifying outliers. See:
#' \cr
#' ?graphics::identify
#' @param devNew If \code{devNew==TRUE} (the default),
#'  \code{dev.new} will be called before each plot.
#'  \cr
#'  This is useful in interactive mode.
#'  \cr
#'  \code{devNew==FALSE} is used for vignette building by \code{package:knitr}.
#'
#' @return There is one point per observation.
#'
#' The following show \bold{probability} \eqn{P_i}{P[i]} on the \eqn{x}-axis:
#'
#' \item{\eqn{P_i \times h_i}{P[i] vs. h[i]}}{
#'  Probability vs. leverage.
#' }
#'
#' \item{\eqn{P_i \times \Delta P \chi^2_i}{P[i] vs. dChisq[i]}}{
#'  Probability vs. the change in the standardized
#'  Pearsons chi-squared
#'  with observation \eqn{i} excluded.
#' }
#' \item{\eqn{P_i \times \Delta D_i}{P[i] vs. dDev[i]}}{
#'  Probability vs. the change in the standardized deviance
#'  with observation \eqn{i} excluded.
#' }
#' \item{\eqn{P_i \times \Delta \hat{\beta}_i}{P[i] vs. dBhat[i]}}{
#'  Probability vs. the change in the standardized
#'  maximum likelihood estimators of the model coefficients
#'  with observation \eqn{i} excluded.
#' }
#' \item{\eqn{P_i \times \Delta P \chi^2_i}{P[i] vs. dChisq[i]}}{
#'  Bubbleplot of
#'  probability vs. the change in the standardized
#'  Pearsons chi-squared
#'  with observation \eqn{i} excluded.
#'  \cr
#'  The area \eqn{A_i}{A[i]} of each circle is
#'  proportional to \eqn{\Delta \hat{\beta}_i}{dBhat[i]}:
#'   \deqn{A_i = \pi r_i^2 \quad r_i = \sqrt{\frac{\Delta \hat{\beta}_i}{P_i}}}{
#'         A[i] = pi r[i]^2, r[i] = (dBhat[i] / P[i])^0.5}
#'  For details see:
#'  \cr
#'  ?graphics::symbols
#' }
#'
#' The following show \bold{leverage} \eqn{h_i}{h[i]} on the \eqn{x}-axis:
#'
#' \item{\eqn{h_i \times \Delta P \chi^2_i}{h[i] vs. dChisq[i]}}{
#'  Leverage vs. the change in the standardized
#'  Pearsons chi-squared
#'  with observation \eqn{i} excluded.
#' }
#' \item{\eqn{h_i \times \Delta D_i}{h[i] vs. dDev[i]}}{
#'  Leverage vs. the change in the standardized deviance
#'  with observation \eqn{i} excluded.
#' }
#' \item{\eqn{h_i \times \Delta \hat{\beta}_i}{h[i] vs. dBhat[i]}}{
#'  Leverage vs. the change in the standardized
#'  maximum likelihood estimators of the model coefficients
#'  with observation \eqn{i} excluded.
#' }
#'
#' The correlation of
#' \eqn{\Delta \chi^2_i, \Delta D_i \mathrm{and} \hat{\beta}_i}{
#'      dChisq, dDev and dBhat}.
#' is shown in a \code{pairs} plot. See:
#' \cr
#' ?graphics::pairs
#'
#' The \bold{Value} of \code{\link{dx}} is also returned, invisibly.
#'
#' @note A choice of colors can be found with e.g.
#' \cr
#' grDevices::colours()[grep("blue", grDevices::colours())]
#'
#' @keywords hplot
#' 
#' @examples
#' ## H&L 2nd ed. Table 4.9. Figures 5.5-5.8. Pages 177-180.
#' data(uis)
#' uis <- within(uis, {
#'     NDRGFP1 <- 10 / (NDRGTX + 1)
#'     NDRGFP2 <- NDRGFP1 * log((NDRGFP1 + 1) / 10)
#' })
#' summary(g1 <- glm(DFREE ~ AGE + NDRGFP1 + NDRGFP2 + IVHX +
#'                   RACE + TREAT + SITE +
#'                   AGE:NDRGFP1 + RACE:SITE,
#'                   family=binomial, data=uis))
#' plot(g1)
#' ## H&L. Similar to Figure 5.3.
#' set.seed(133)
#' (g1 <- glm(sample(c(0, 1), size=100,
#'                   replace=TRUE, prob=c(0.5, 0.5))
#'            ~ 0 + I(0.08 * rnorm(n=100, mean=0, sd=sqrt(9))),
#'            family=binomial))$coef # approx. 0.8
#' plot(g1)
plot.glm <- function(x,
                     y=NULL,
                     ...,
                     toPdf=FALSE,
                     file="dxPlots.pdf",
                     palette=c("Dark2", "Set2", "Accent", "Blues"),
                     usePalette=TRUE,
                     bg=NULL,
                     col="white",
                     alpha=0.4,
                     cex=2,
                     pch=21,
                     cex.main=1.5,
                     inches=0.25,
                     identify=FALSE,
                     devNew=TRUE){
    stopifnot(inherits(x, "glm"))
    stopifnot(x$family$family=="binomial")
    if (toPdf && identify) stop (
        "Cannot use 'identify' (interactive) when writing to a file")
    ## for R CMD check
    n <- P <- h <- dChisq <- dBhat <- dDev <- NULL
    bg <- if(usePalette){
        palette <- match.arg(palette)
        suppressWarnings(RColorBrewer::brewer.pal(n=Inf,
                                                  name=palette))
    } else {
        bg
    }
    bg <- grDevices::adjustcolor(bg, alpha=alpha)
    if(toPdf) grDevices::pdf(file)
    ## get diagnostics for model
    dx1 <- dx(x)
    opar <- graphics::par
    on.exit(graphics::par(opar))
    on.exit(if(toPdf) grDevices::dev.off())
    ## bottom, left, top, right
    b1 <- 4
    l1 <- 7
    t1 <- 6
    r1 <- 0.5
    graphics::par(mar=c(b1, l1, t1, r1))
###
###-------------------------------------------
###
    ## probabiility by leverage
    m1 <- expression(
        atop(
            paste("Probability ",
                  italic(P[i]) %*% "leverage ",
                  italic(h[i])),
            atop(
                paste(
                    paste(0.1 < italic(P[i])) > 0.9,
                    "   " %->% "   ",
                    italic(h[i])
                    %prop%
                    paste(x[i] - mu[x])),
                atop(
                    paste(
                        italic(h[i]) %~~%
                        " distance of covariate pattern ",
                        x[i], " from mean ", mu[x]),
                    italic(h[i]) == "diagonal of hat matrix"))))
    graphics::plot.default(dx1[, P], dx1[, h],
                           main=m1,
                           cex.main=cex.main,
                           col=col,
                           bg=bg,
                           pch=pch,
                           cex=cex,
                           xlab="",
                           ylab="",
                           ...)
    graphics::mtext(expression(italic(h[i])),
                    side=2, line=3, las=1, cex=cex.main)
    graphics::mtext(expression(paste("probability ", italic(P[i]))),
                    side=1, line=3, cex=cex.main)
    if (identify) graphics::identify(dx1[, P], dx1[, h])
###
### reset plotting parameters each time
    dn1 <- function(b1=b1, l1=l1, t1=t1, r1=r1){
        grDevices::dev.new()
        graphics::par(mar=c(b1, l1, t1, r1))
    }
###
###-------------------------------------------
### Probability by ...
    ## dXsq
    if(!toPdf && devNew) dn1(b1, l1, t1, r1)
    m1 <- expression(
        atop(
            paste("Probability ",
                  italic(P[i]) %*%
                  paste("scaled change in Pearson chi-sq ",
                        paste(s, Delta, P, chi[i] ^2))),
            list(paste(italic(Pr[i]) ==
                       frac(paste(italic(y[i] - mu[y])),
                            italic(sigma[y]))),
                 paste(paste(s, Delta, P, chi[i] ^2),
                       " = ",
                       paste(frac(italic(Pr[i]),
                                  sqrt(paste(1 - italic(h[i])))))))))
    graphics::plot.default(dx1[, P], dx1[, dChisq],
                           main=m1,
                           cex.main=cex.main,
                           xlab="",
                           ylab="",
                           col=col,
                           bg=bg,
                           pch=pch,
                           cex=cex,
                           ...)
    graphics::mtext(expression(paste(s, Delta, P, chi[i] ^2)),
                    side=2, line=3, las=1, cex=cex.main)
    graphics::mtext(expression(paste("probability ", italic(P[i]))),
                    side=1, line=3, cex=cex.main)
    if (identify) graphics::identify(dx1[, P], dx1[, dChisq])
    ## dDev
    if(!toPdf && devNew) dn1(b1, l1, t1, r1)
    m1 <- expression(
        atop(
            paste("Probability ",
                  italic(P[i]) %*%
                  paste("scaled change in deviance ",
                        paste(Delta, D[i]))),
            list(paste(italic(dr[i]) ==
                       paste(sign, "(", italic(y[i]) - italic(hat(y)[i]), ")",
                             sqrt(italic(d[i])))),
                 paste(paste(s, Delta, D[i]),
                       " = ",
                       paste(frac(italic(dr[i]),
                                  sqrt(paste(1 - italic(h[i])))))))))
    graphics::plot.default(dx1[, P], dx1[, dDev],
                           main=m1,
                           cex.main=cex.main,
                           xlab="",
                           ylab="",
                           col=col,
                           bg=bg,
                           pch=pch,
                           cex=cex,
                           ...)
    graphics::mtext(expression(paste(s, Delta, D[i])),
                    side=2, line=3, las=1, cex=cex.main)
    graphics::mtext(expression(paste("probability ", italic(P[i]))),
                    side=1, line=3, cex=cex.main)
    if (identify) graphics::identify(dx1[, P], dx1[, dDev])
    ## dBhat
    if (!toPdf && devNew) dn1(b1, l1, t1, r1)
    m1 <- expression(
        atop(
            paste("Probability " ,
                  italic(P[i]) %*%
                  paste("scaled change in coefficients ",
                        paste(s, Delta, hat(beta)[i]))),
            paste(s, Delta, hat(beta[i]),
                  " = ",
                  frac(italic(sPr[i] ^2 * h[i]),
                       1 - italic(h[i])))))
    graphics::plot.default(dx1[, P], dx1[, dBhat],
                           main=m1,
                           cex.main=cex.main,
                           xlab="",
                           ylab="",
                           col=col,
                           bg=bg,
                           pch=pch,
                           cex=cex,
                           ...)
    graphics::mtext(expression(paste(s, Delta, hat(beta)[i])),
                    side=2, line=3, las=1, cex=cex.main)
    graphics::mtext(expression(paste("probability ", italic(P[i]))),
                    side=1, line=3, cex=cex.main)
    if (identify) graphics::identify(dx1[, P], dx1[, dBhat])
    ##
    ##-------------------------------------------
    ## bubble plot - prob by dChisq, with area = dBhat
    if(!toPdf && devNew) dn1(b1, l1, t1, r1)
    radius1 <- sqrt(dx1[, dBhat] / dx1[, P])
    m1 <- expression(
        atop(
            paste("Probability ",
                  italic(P[i]) %*%
                  paste("scaled change in Pearson chi-sq ",
                        paste(s, Delta, P, chi[i] ^2))),
            list(
                area %prop% paste(s, Delta, hat(italic(beta[i]))),
                radius == sqrt(frac(
                paste(s, Delta, hat(beta[i])), italic(P[i]))))))
    graphics::symbols(dx1[, P], dx1[, dChisq],
                      main=m1,
                      cex.main=cex.main,
                      circles=radius1,
                      inches=inches,
                      fg=col,
                      bg=bg,
                      xlab="",
                      ylab="",
                      ...)
    graphics::mtext(expression(paste(s, Delta, P, chi[i] ^2)),
                    side=2, line=3, las=1, cex=cex.main)
    graphics::mtext(expression(paste("probability ", italic(P[i]))),
                    side=1, line=3, cex=cex.main)
    if (identify) graphics::identify(dx1[, P], dx1[, dChisq])
###
###---------------------------
### leverage by...
    ## dChisq
    if(!toPdf && devNew) dn1(b1, l1, t1, r1)
    m1 <- expression(
        atop(
            paste("Leverage ",
                  italic(h[i]) %*%
                  paste("scaled change in Pearson chi-sq ",
                        paste(s, Delta, P, chi[i] ^2))),
            list(
                italic(h[i]) %~~% x[i] - mu[x],
                paste(paste(s, Delta, P, chi[i]^2) == italic(sPr[i])^2))))
    graphics::plot.default(dx1[, h], dx1[, dChisq],
                           main=m1,
                           cex.main=cex.main,
                           col=col,
                           bg=bg,
                           pch=pch,
                           cex=cex,
                           xlab="",
                           ylab="",
                           ...)
    graphics::mtext(expression(paste(s, Delta, P, chi[i] ^2)),
                    side=2, line=3, las=1, cex=cex.main)
    graphics::mtext(expression(paste("leverage ", italic(h[i]))),
                    side=1, line=3, cex=cex.main)
    if (identify) graphics::identify(dx1[, h], dx1[, dChisq])
    ## dDev
    if(!toPdf && devNew) dn1(b1, l1, t1, r1)
    m1 <- expression(
        atop(
            paste("Leverage ",
                  italic(h[i]) %*%
                  paste("scaled change in deviance ",
                        paste(s, Delta, D[i]))),
            list(paste(italic(dr[i]) ==
                       paste(sign, "(", italic(y[i]) - italic(hat(y)[i]), ")",
                             sqrt(italic(d[i])))),
                 paste(paste(s, Delta, D[i]),
                       " = ",
                       paste(frac(italic(dr[i]),
                                  sqrt(paste(1 - italic(h[i])))))))))
    graphics::plot.default(dx1[, h], dx1[, dDev],
                           main=m1,
                           cex.main=cex.main,
                           col=col,
                           bg=bg,
                           pch=pch,
                           cex=cex,
                           xlab="",
                           ylab="",
                           ...)
    graphics::mtext(expression(paste(s, Delta, D[i])),
                    side=2, line=3, las=1, cex=cex.main)
    graphics::mtext(expression(paste("leverage ", italic(h[i]))),
                    side=1, line=3, cex=cex.main)
    if (identify) graphics::identify(dx1[, h], dx1[, dDev])
    ## dBhat
    if(!toPdf && devNew) dn1(b1, l1, t1, r1)
    m1 <- expression(
        atop(
            paste("Leverage ",
                  italic(h[i]) %*%
                  paste("scaled change in coefficients ",
                        paste(s, Delta, hat(beta)[i]))),
            list(
                italic(h[i]) %~~% x[i] - bar(x),
                paste(s, Delta, hat(beta[i]),
                      " = ",
                      frac(italic(sPr[i]^2 * h[i]),
                           1 - italic(h[i]))))))
    graphics::plot.default(dx1[, h], dx1[, dBhat],
                           main=m1,
                           cex.main=cex.main,
                           col=col,
                           bg=bg,
                           pch=pch,
                           cex=cex,
                           xlab="",
                           ylab="",
                           ...)
    graphics::mtext(expression(paste(s, Delta, hat(beta)[i])),
                    side=2, line=3, las=1, cex=cex.main)
    graphics::mtext(expression(paste("leverage ", italic(h[i]))),
                    side=1, line=3, cex=cex.main)
    if (identify) graphics::identify(dx1[, h], dx1[, dBhat])
### pairs plot
    m1 <- expression(
        paste("Correlation between ",
              paste(s, Delta, P, chi[i] ^2),
              ", ",
              paste(s, Delta, D[i]),
              " and ",
              paste(s, Delta, hat(beta)[i])))
    graphics::pairs(dx1[ ,list(dChisq, dDev, dBhat)],
                    main=m1,
                    cex.main=cex.main,
                    labels=c(
                        expression(paste(s, Delta, P, chi[i] ^2)),
                        expression(paste(s, Delta, D[i])),
                        expression(paste(s, Delta, hat(beta)[i]))),
                    bg=bg,
                    col=col,
                    pch=pch,
                    cex=cex)
###
    return(invisible(dx1))
}
