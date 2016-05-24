#' Nonlinear regression function
#' 
#' Requires \code{plotSEMM_setup} be run first.  Generates (a) the potential nonlinear 
#' regression function; (b) bivariate distribution of the latent variables; 
#' (c) marginal distributions of the latent variables; (d) within class linear 
#' regression functions; and (e) within class marginal distributions for the latent variables. 
#' 
#' @aliases plotSEMM_contour
#' @param SEMLIdatapks object returned from \code{\link{plotSEMM_setup}}
#' @param EtaN2 Label for the X axis.  If no value is provided, defaults to "Eta2." 
#' @param EtaN1 Label for the Y axis.  If no value is provided, defaults to "Eta1." 
#' @param classinfo Logical variable. TRUE shows the lines for each class as well as the combined estimate.  
#'   FALSE shows only the combined estimate.  If no value is provided, defaults to TRUE. 
#' @param lnty Determines the line types used for the class lines.  If no value is provided, defaults to 3. 
#'   See \code{\link{par}} for information about line type.
#' @param lncol Determines the line colors used for the class lines.  If no value is provided, defaults to 1. 
#'   See \code{\link{par}} for information about line type. 
#' @param title Titles the graph. 
#' @param leg Logical variable.  If TRUE, a legend accompanies the graph.  If FALSE, no legend appears. 
#'   Defaults to TRUE.
#' @param cex par(cex) value. Default is 1.5
#' @param ... addition inputs, mostly from plotSEMM_GUI()
#' @author Bethany Kok and Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords hplot color
#' @export plotSEMM_contour
#' @examples 
#' \dontrun{
#' ## code for latent variables with two classes
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
#' 
#' plotSEMM_contour(plotobj)
#' 
#' plotSEMM_contour(plotobj, EtaN1 = "Latent Predictor", 
#'    EtaN2 = "Latent Outcome", classinfo = FALSE, lncol = 5) 
#' }
plotSEMM_contour <- function(SEMLIdatapks, EtaN2 = "Eta2", EtaN1 = "Eta1", 
                             classinfo = TRUE, lnty = 3, lncol = 1, title = "", 
                             leg = TRUE, cex = 1.5, ...) {
    dots <- list(...)
    input <- dots$input
    if(!is.null(input$xlab)) EtaN1 <- input$xlab
    if(!is.null(input$ylab)) EtaN2 <- input$ylab
    legend_location <- if(!is.null(input$legend_location)) input$legend_location else 'topright'
    if(legend_location == 'default') legend_location <- 'topright'
    if(legend_location == 'none') leg <- FALSE
    if(!is.null(input$class_info)) classinfo <- input$class_info
    if(!classinfo) leg <- FALSE
    
    # This is to make all 3 graphs on 1 panel for contour
    def.par <- par(no.readonly = TRUE)  # save default, for resetting... 
    nf <- layout(matrix(c(1, 0, 2, 3), 2, 2, byrow = TRUE), c(3, 1), c(1, 3), TRUE)
    
    xlabel = deparse(substitute(EtaN1))
    if (substring(xlabel, 1, 1) == "\"") {
        xlabel = substring(xlabel, 2, nchar(xlabel) - 1)
    }
    
    ylabel = deparse(substitute(EtaN2))
    if (substring(ylabel, 1, 1) == "\"") {
        ylabel = substring(ylabel, 2, nchar(ylabel) - 1)
    }
    
    maintitle = deparse(substitute(title))
    if (substring(maintitle, 1, 1) == "\"") {
        maintitle = substring(maintitle, 2, nchar(maintitle) - 1)
    }
    
    
    # plot1 Exogenous
    par(mar = c(0, 5, 1, 1))
    plot(SEMLIdatapks$Eta1, SEMLIdatapks$agg_denEta1, type = "l", xlab = "", ylab = "", main = maintitle, axes = FALSE, cex.lab=cex, cex.axis=cex)
    if (classinfo == TRUE) {
        for (i in 1:SEMLIdatapks$classes[1]) {
            lines(SEMLIdatapks$Eta1, SEMLIdatapks$class_denEta1[, i], lwd = 1, lty = (i + lnty), col = (i + lncol))
        }
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
        
        legend(x = legend_location, legend = text, horiz = FALSE, lwd = lwd1, 
                   lty = lty1, col = col1, , bty = "n", cex=cex)
    }
    
    # plot3 contour
    par(mar = c(4, 5, 0, 0))
    plot(SEMLIdatapks$Eta1, SEMLIdatapks$Eta2, type = "n", xlab = xlabel, ylab = ylabel, main = "", cex.lab=cex, cex.axis=cex)
    contour(SEMLIdatapks$Eta1, SEMLIdatapks$Eta2, SEMLIdatapks$contour, drawlabels = TRUE, add = TRUE, nlevels = 20)
    lines(SEMLIdatapks$Eta1, SEMLIdatapks$agg_pred, lwd = 2, lty = 1)
    if (classinfo == TRUE) {
        for (i in 1:SEMLIdatapks$classes[1]) {
            lines(SEMLIdatapks$Eta1, SEMLIdatapks$class_pred[, i], lwd = 1, lty = (i + lnty), col = (i + lncol))
        }
    }
    
    # plot Endogenous
    par(mar = c(4, 0, 1, 1))
    plot(SEMLIdatapks$agg_denEta2, SEMLIdatapks$Eta2, type = "l", xlab = "", ylab = "", main = "", axes = FALSE, cex.lab=cex, cex.axis=cex)
    if (classinfo == TRUE) {
        for (i in 1:SEMLIdatapks$classes[1]) {
            lines(SEMLIdatapks$class_denEta2[, i], SEMLIdatapks$Eta2, lwd = 1, lty = (i + lnty), col = (i + lncol))
        }
    }
    
    
    
    par(def.par)  #reset to default
}
 
