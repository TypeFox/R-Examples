#' Probability plot
#'
#' Requires \code{plotSEMM_setup} be run first. Generates a plot which expresses
#' the mixing probabilities for each latent class conditioned on the latent predictor.
#'
#' @aliases plotSEMM_probability
#' @param SEMLIdatapks object returned from \code{\link{plotSEMM_setup}}
#' @param EtaName Label of the latent predictor.  If no value is provided, defaults to Eta1.
#' @param lnty Determines the line types used for the class lines.  If no value is provided,
#'   defaults to 3.  See \code{\link{par}} for information about line type.
#' @param lncol Determines the line colors used for the class lines.  If no value is
#'   provided, defaults to 1.  See \code{\link{par}} for information about line type.
#' @param title Titles the graph.
#' @param leg Logical variable.  If TRUE, a legend accompanies the graph.  If FALSE,
#'   no legend appears.  Defaults to TRUE.
#' @param cex par(cex) value. Default is 1.5
#' @param ... addition inputs, mostly from plotSEMM_GUI()
#' @author Bethany Kok and Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords hplot color
#' @export plotSEMM_probability
#' @seealso \code{\link{plotSEMM_setup}}, \code{\link{plotSEMM_contour}}
#' @examples
#' \dontrun{
#' # 2 class empirical example on positive emotions and heuristic processing in
#' # Pek, Sterba, Kok & Bauer (2009)
#' pi <- c(0.602, 0.398)
#'
#' alpha1 <- c(3.529, 2.317)
#'
#' alpha2 <- c(0.02, 0.336)
#'
#' beta21 <- c(0.152, 0.053)
#'
#' psi11 <- c(0.265, 0.265)
#'
#' psi22 <- c(0.023, 0.023)
#'
#'
#' plotobj <- plotSEMM_setup(pi, alpha1, alpha2, beta21, psi11, psi22)
#'
#' plotSEMM_probability(plotobj)
#'
#' plotSEMM_probability(plotobj , EtaName = "Latent Predictor", lnty = 2, title = "Probability")
#' }
plotSEMM_probability <- function(SEMLIdatapks, EtaName = "Eta1", lnty = 3, lncol = 1,
                                 title = "", leg = TRUE, cex = 1.5, ...) {
    dots <- list(...)
    input <- dots$input
    if(!is.null(input$xlab)) EtaName <- input$xlab
    legend_location <- if(!is.null(input$legend_location)) input$legend_location else 'topright'
    if(legend_location == 'default') legend_location <- 'topright'
    if(legend_location == 'none') leg <- FALSE

    # plot 2 probabilities and exogenous variable
    def.par <- par(no.readonly = TRUE)  # save default, for resetting...
    nf <- layout(matrix(c(1, 2), 2, 1, byrow = TRUE), c(3, 1), c(1, 3), TRUE)
    layout.show(nf)  #make sure it's the format we want

    maintitle = deparse(substitute(title))
    if (substring(maintitle, 1, 1) == "\"") {
        maintitle = substring(maintitle, 2, nchar(maintitle) - 1)
    }

    # plot1 Exogenous
    par(mar = c(0, 4, 1, 1))
    plot(SEMLIdatapks$Eta1, SEMLIdatapks$agg_denEta1, type = "l", xlab = "", ylab = "", main = maintitle, axes = FALSE, cex.lab=cex, cex.axis=cex)
    for (i in 1:SEMLIdatapks$classes[1]) {
        lines(SEMLIdatapks$Eta1, SEMLIdatapks$class_denEta1[, i], lwd = 1, lty = (i + lnty), col = (i + lncol))
    }

    # legend
    if (leg == TRUE) {
        classes <- SEMLIdatapks$classes[1]

        text <- vector(mode = "character", length = classes + 1)
        text[1] <- "Estimate"

        lwd1 <- vector(mode = "numeric", length = (classes + 1))
        lwd1[1] <- 2

        for (i in 2:(classes + 1)) {
            text[i] <- paste("Class", (i - 1), sep = " ")
            lwd1[i] <- 1
        }

        lty1 <- vector(mode = "numeric", length = (classes + 1))

        col1 <- vector(mode = "numeric", length = (classes + 1))

        col1[1] <- 1
        lty1[1] <- 1
        for (i in 2:(classes + 1)) {
            col1[i] <- (i + lncol - 1)
            lty1[i] <- (i + lnty - 1)
        }


        legend(x = legend_location, legend = text, horiz = FALSE, lwd = lwd1, lty = lty1, col = col1, , bty = "n", cex=cex)
    }

    # plot2 Conditional Probabilities
    par(mar = c(5, 5, 0, 0))
    xlabel = deparse(substitute(EtaName))
    if (substring(xlabel, 1, 1) == "\"") {
        xlabel = substring(xlabel, 2, nchar(xlabel) - 1)
    }

    ylabel = paste(paste("Probability(Class | ", xlabel, sep = ""), ")", sep = "")

    plot(SEMLIdatapks$Eta1, SEMLIdatapks$class_prob[, 1], type = "n", xlab = xlabel, ylab = ylabel, main = "", cex.lab=cex, cex.axis=cex)
    for (i in 1:SEMLIdatapks$classes[1]) {
        lines(SEMLIdatapks$Eta1, SEMLIdatapks$class_prob[, i], lwd = 1, lty = (i + lnty), col = (i + lncol))
    }


    par(def.par)  #reset to default
}

