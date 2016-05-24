#' Plot score for prcomp result
#' 
#' Plot score for \code{\link[stats]{prcomp}} (PCA) result
#' 
#' @param prcompResult output object from \code{\link[stats]{prcomp}} function
#' @param xPC an integer indicating PC component on x-axis
#' @param yPC an integer indicating PC component on y-axis
#' @param cex (optional) size of points on graphs 
#' @param cex.legend (optional) size of fonts in legend 
#' @param label (optional) a character vector or expression specifying the text to be written.
#' @param pos (optional, applicable when label is given) a position specifier for the text. If specified this overrides 
#' any adj value given. Values of 1, 2, 3 and 4, respectively indicate positions below, 
#' to the left of, above and to the right of the specified coordinates.
#' @param group a vector of numeric, character or factor class separating 
#' the samples into groups. Correspond to point color. 
#' @param group2 The second group, can be a vector of numeric, character or factor class separating 
#' the samples into groups. Correspond to point shape. 
#' @param legendlocation (optional)location of legend on graph. 
#' Look up \code{\link[graphics]{legend}} for more details.
#' @param legendoutside (optional) set to TRUE if you want to put legend on 
#' the outside of the plot. The legend location is defaulted to topright. 
#' @param rightwhitespace (optional) set width for white space for legend. 
#' Only applicable if legendoutside = TRUE
#' @param col point color palette
#' @param pch point type palette
#' @param ... additional arguments for \code{\link[graphics]{par}}
#' 
#' @return A figure is returned on the graphic device
#' 
#' @seealso
#' \code{\link{plotScorem}}
#' 
#' @examples
#' data(applejuice)
#' applejuice_uf <- unfold(applejuice) # unfold list into matrix
#' result <- prcomp(applejuice_uf) 
#' plotScore(result) # plot PC1 vs PC2 score
#' plotScore(result, pch = 3, col = "blue") # change shape and color
#' 
#' # get country of apple production
#' country <- sapply(strsplit(names(applejuice), split = "-"), "[", 1) 
#' plotScore(result, label = country) # add label
#' 
#' # or plot by group
#' plotScore(result, xPC = 1, yPC = 3, group = country) 
#' 
#' # custom point types and color
#' plotScore(result, xPC = 1, yPC = 3, group = country, pch = c(1,2), col = c("green", "black"))
#' 
#' # move legend outside
#' plotScore(result, xPC = 1, yPC = 3, group = country, legendoutside = TRUE)
#' 
#' # two groups
#' cultivar <- sapply(strsplit(names(applejuice), split = "-"), "[", 2) 
#' plotScore(result, group = country, group2 = cultivar)
#' 
#' # make the points more transparent
#' \dontrun{
#' require(scales)
#' plotScore(result, group = country, group2 = country, col = alpha(generateColor(2), 0.7))
#' }
#' 
#' @importFrom graphics abline legend par plot points text
#' @export
#' 
plotScore <-
    function(prcompResult, xPC = 1, yPC = 2, group = NULL, group2 = NULL,
             cex = 1.5, cex.legend = 1, label = NULL, pos = 4, col = NULL, pch = NULL,
             legendlocation = "bottomright",
             legendoutside = FALSE,
             rightwhitespace = 0,
             ...){
        
        # check if group information is provided
        has_group <- !is.null(group)
        has_group2 <- !is.null(group2)
        
        # get information from prcompResult
        score <- prcompResult$x
        numSample <- dim(score)[1]
        
        # check validity of group 
        if (has_group) {
            
            # check length
            if (length(group) %% numSample != 0) {
                stop("The dimension of group and sample do not match. 
           Please check your group variable.")
            }
            
            # turn into factor if it isn't already
            if (!is.factor(group) & is.null(attributes(group)$levels)) {
                group <- as.factor(group)
            }
            numLevels_color <- nlevels(group)
        } 
        
        # do the same for group2
        if (has_group2) {
            
            # check length
            if (length(group2) %% numSample != 0) {
                stop("The dimension of group2 and sample do not match. 
           Please check your group2 variable.")
            }
            
            # turn into factor if it isn't already
            if (!is.factor(group2) & is.null(attributes(group2)$levels)) {
                group2 <- as.factor(group2)
            }
            numLevels_point <- nlevels(group2)
        } else if (has_group){
            # if there is only group supplied, let group2 be group
            group2 <- group
            numLevels_point <- nlevels(group2)
        }
        
        if (is.null(group)) group <- 1
        if (is.null(group2)) group2 <- 1
        if (!(exists("numLevels_color"))) numLevels_color <- 1
        if (!(exists("numLevels_point"))) numLevels_point <- 1
        
        # color and point type base on group
        col.palette <- generateColor(numLevels_color, if (!is.null(col))col)
        col <- col.palette[group]
        
        pch.palette <- generatePoint(numLevels_point, if (!is.null(pch))pch)
        pch <- pch.palette[group2]
        
        
        # prepare plotting information
        xLabel <- prcompname(prcompResult, xPC)
        yLabel <- prcompname(prcompResult, yPC)
        
        # set plotting area when there is legend outside
        if (isTRUE(legendoutside)) {
            # modify space if legendoutside == TRUE
            par(mar = c(5.1, 4.1, 4.1, 2.1+6+rightwhitespace))
        }
        
        # plot
        plot(score[, xPC], score[, yPC], xlab = xLabel, ylab = yLabel,
             cex = cex, col = col, pch = pch, ...)
        
        abline(v = 0, h = 0, lty = 2, col = "grey39")
        
        # put in label if label is provided
        if (!is.null(label)){
            
            # check input
            if (length(label) != dim(score)[1]) {
                stop("The dimension of label and sample do not match. \n
             Please check your label variable.")
            }
            text(score[, xPC], score[, yPC], labels = label,
                 pos = pos)
        }
        
        # legend for group
        if (has_group | has_group2){
            # turn unclassed factor back 
            group <- levels(group) 
            group2 <- levels(group2)
            
            if (has_group & has_group2) {
                legend_labels <- c(group, NA, group2)
                col.palette <- c(col.palette, NA, rep("black", numLevels_point))
                pch.palette <- c(rep(45, numLevels_color), NA, pch.palette)
            } 
            if (has_group & !has_group2) legend_labels <- group
            if (!has_group & has_group2)  legend_labels <- group2
            
            # plot legend inside or outside
            if (legendoutside){
                # plot legend outside
                # legendlocation will be overwritten if provided
                par(xpd = T)
                legend(par()$usr[2], par()$usr[4], legend = legend_labels, pch = pch.palette,
                           pt.cex = cex, col = col.palette, xpd = TRUE, cex = cex.legend)
                
                # reset mar 
                par(mar = c(5.1, 4.1, 4.1, 2.1))
                
            } else {
                # plot legend inside
                legend(legendlocation, legend = legend_labels, pch = pch.palette,
                       pt.cex = cex, col = col.palette, xpd = TRUE, cex = cex.legend)
            }
        }
        par(xpd = FALSE)
    }
