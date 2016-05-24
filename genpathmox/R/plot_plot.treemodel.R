
#'@title Comparative plot between nodes from the Pathmox Segmentation Trees: PLS-PM
#'
#'@description
#'Plot method for objects of class \code{"treemodel"}. Barplots of path
#'coefficients of terminal nodes with respect to those of the global (root)
#'model
#'
#'@details
#'This function aims to visualize the comparison between path coefficients of
#'the terminal nodes against the path coefficients of the global model in the
#'root node. \cr When \code{comp.by="nodes"} a graphic window is displayed for
#'each endogenous latent variable of the PLS model, and barplots of nodes are
#'shown. \cr When \code{comp.by="latents"} a graphic window is displayed for
#'each endogenous relationship of the PLS model, and barplots of independent
#'latent variables are shown.
#'
#'@param x An object of class \code{"treemodel"} returned by
#'\code{\link{pls.treemodel}}.
#'@param comp.by One of "nodes" or "latents". This argument indicates the type
#'of barplots comparison.
#'@param nodes.names Optional vector of names for the terminal nodes (must be a
#'vector of length equal to the number of terminal nodes).
#'@param ordered A logical value indicating whether the barplots are shown in
#'increasing ordered.
#'@param decreasing A logical value indicating if the sort order should be
#'increasing or decreasing.
#'@param color Optional vector of colors for the bars. When \code{color=NULL}
#'rainbow colors are used.
#'@param show.box A logical value indicating whether a box is drawn around each
#'barplot.
#'@param border The color to be used for the border of the bars. Use
#'\code{border=NA} to omit borders.
#'@param cex.names Expansion factor for axis names (bar labels).
#'@param cex.axis Expansion factor for numeric axis labels.
#'@param short.labs Logical value indicating if the labels of the barplots
#'should be abbreviated (\code{TRUE} by default).
#'@param short.min Integer number indicating the minimum length of the
#'abbreviations for the labels. Only used when \code{short.labs=TRUE}.
#'@param cex.main Allows to fix the size of the main. Equal to 1 to default 
#' @param \dots Further arguments are ignored.
#'
#' @author Giuseppe Lamberti
#' 
#' @references Sanchez, G. (2009) \emph{PATHMOX Approach: Segmentation Trees in
#' Partial Least Squares Path Modeling.} PhD Dissertation. 
#' 
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#' @method plot treemodel
#' @S3method plot treemodel
#' 
#'@examples
#'  \dontrun{
#'  ## example of PLS-PM in alumni satisfaction
#'  
#'  data(fibtele)
#'  
#'  # select manifest variables
#'  data.fib <-fibtele[,12:35]
#'  
#'  # define inner model matrix
#'  Image 			= rep(0,5)
#'	 Qual.spec	= rep(0,5)
#'	 Qual.gen		= rep(0,5)
#'	 Value			= c(1,1,1,0,0)
#'	 Satis			= c(1,1,1,1,0)
#'  inner.fib <- rbind(Image,Qual.spec, Qual.gen, Value, Satis)
#'  colnames(inner.fib) <- rownames(inner.fib)
#'  
#'  # blocks of indicators (outer model)
#'  outer.fib <- list(1:8,9:11,12:16,17:20,21:24)
#'  modes.fib  = rep("A", 5)
#'  
#'  # apply plspm
#'  pls.fib <- plspm(data.fib, inner.fib, outer.fib, modes.fib)
#'                  
#'
#'  # re-ordering those segmentation variables with ordinal scale 
#'   seg.fib= fibtele[,2:11]
#'  
#'	 seg.fib$Age = factor(seg.fib$Age, ordered=T)
#'	 seg.fib$Salary = factor(seg.fib$Salary, 
#'		levels=c("<18k","25k","35k","45k",">45k"), ordered=T)
#'	 seg.fib$Accgrade = factor(seg.fib$Accgrade, 
#'		levels=c("accnote<7","7-8accnote","accnote>8"), ordered=T)
#'	 seg.fib$Grade = factor(seg.fib$Grade, 
#'		levels=c("<6.5note","6.5-7note","7-7.5note",">7.5note"), ordered=T)
#'
#'  # Pathmox Analysis
#'  fib.pathmox=pls.pathmox(pls.fib,seg.fib,signif=0.05,
#'				deep=2,size=0.2,n.node=20)
#' 
#'  fib.comp=pls.treemodel(pls.fib,fib.pathmox)
#'  plot(fib.comp)
#'  
#'  }
#'

plot.treemodel	<- function (x, comp.by = "nodes", nodes.names = NULL, ordered = TRUE, 
    decreasing = FALSE, color = NULL, show.box = TRUE, border = NA, 
    cex.names = 0.75, cex.axis = 0.75, short.labs = TRUE, short.min = NULL,cex.main=1, 
    ...) 
{
    if (length(union(comp.by, c("nodes", "latents"))) > 2) 
        stop("Invalid argument 'comp.by'")
    if (!is.null(nodes.names)) {
        if (length(nodes.names) != (ncol(x$paths) - 1)) {
            cat("NOTICE: Invalid argument 'nodes.names'. Default names were used", 
                "\n")
        }
        else {
            colnames(x$paths) <- c("Root_Node", nodes.names)
        }
    }
    if (short.labs) 
        if (mode(short.min) != "numeric" || length(short.min) != 
            1 || (short.min%%1) != 0) 
            short.min <- 7
    IDM <- x$IDM
    lvs <- nrow(IDM)
    idm <- as.vector(IDM)
    idm[idm == 1] <- 1:sum(IDM)
    SDM <- matrix(idm, lvs, lvs)
    B <- x$paths[, -1] - x$paths[, 1]
    rs <- c(1, 1, 1, 2, 2, 2, 2, 2, 3, 2, 3, 3, 4, 4)
    cs <- c(1, 2, 3, 2, 3, 3, 4, 4, 3, 5, 4, 4, 4, 4)
    index.mat <- cbind(1:14, rs, cs)
    endo <- rowSums(IDM)
    endo[endo != 0] <- 1
    who.endo <- which(endo != 0)
    if (comp.by == "nodes") {
        nn <- ncol(B)
        if (is.null(color)) 
            colors = rainbow(lvs, s = 0.5, v = 0.7)
        else colors = rep(color, length.out = lvs)
        for (i in who.endo) {
            n.ind <- which(IDM[i, ] == 1)
            indep <- SDM[i, which(SDM[i, ] != 0)]
            arg.names <- colnames(IDM)[n.ind]
            if (short.labs) 
                arg.names <- abbreviate(arg.names, minlength = short.min)
            cols <- colors[n.ind]
            dev.new()
            par(mfrow = index.mat[nn + 1, 2:3])
            par(mar = c(3, 3, 3, 3))
            barplot(x$paths[indep, 1], names.arg = arg.names, 
                main = paste("Global:", rownames(IDM)[i]), col = cols, 
                border = border,cex.axis=cex.axis,cex.names=cex.names)
            abline(h = 0)
            if (show.box) 
                box("figure")
            ylim <- 1.15 * c(min(B[indep, ]), max(B[indep, ]))
            for (n in 1:nn) {
                if (ordered) {
                  barplot(sort(B[indep, n], decreasing), names.arg = arg.names[order(B[indep, 
                    n], decreasing = decreasing)], main = colnames(B)[n], 
                    border = border, col = cols[order(B[indep, 
                      n], decreasing = decreasing)], cex.names = cex.names, 
                    cex.main = cex.main, cex.axis = cex.axis, ylim = ylim, 
                    ...)
                }
                else {
                  barplot(B[indep, n], names.arg = arg.names, 
                    main = colnames(B)[n], border = border, col = cols, 
                    cex.names = cex.names, cex.main = cex.main, cex.axis = cex.axis, 
                    ylim = ylim, ...)
                }
                abline(h = 0)
                if (show.box) 
                  box("figure")
            }
        }
    }
    else {
        ncols <- ncol(B)
        if (is.null(color)) 
            colors = rainbow(ncols, s = 0.5, v = 0.7)
        else colors = rep(color, length.out = ncols)
        arg.names <- colnames(B)
        if (short.labs) 
            arg.names <- abbreviate(arg.names, minlength = short.min)
        for (i in who.endo) {
            indep <- SDM[i, which(SDM[i, ] != 0)]
            dev.new()
            par(mfrow = index.mat[length(indep), 2:3])
            par(mar = c(3, 3, 3, 3))
            ylim <- 1.15 * c(min(B[indep, ]), max(B[indep, ]))
            for (k in indep) {
                if (ordered) {
                  barplot(sort(B[k, ], decreasing), names.arg = arg.names[order(B[k, 
                    ], decreasing = decreasing)], main = rownames(B)[k], 
                    border = NA, cex.main = cex.main, cex.names = cex.names, 
                    cex.axis = cex.axis, ylim = ylim, col = colors[order(B[k, 
                      ], decreasing = decreasing)], ...)
                }
                else {
                  barplot(B[k, ], names.arg = arg.names, main = rownames(B)[k], 
                    border = NA, cex.main = cex.main, cex.names = cex.names, 
                    cex.axis = cex.axis, ylim = ylim, col = colors, 
                    ...)
                }
                abline(h = 0)
                if (show.box) 
                  box("figure")
            }
        }
    }
}
