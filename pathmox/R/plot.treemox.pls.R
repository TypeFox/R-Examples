#'@title Comparative plot between nodes from a PATHMOX or TECHMOX tree
#'
#'@description
#'Plot method for objects of class \code{"treemox.pls"}. Barplots of path
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
#'@param x An object of class \code{"treemox.pls"} returned by
#'\code{\link{treemox.pls}}.
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
#'@param \dots Arguments to be passed to/from other methods.
#'@seealso \code{\link{treemox.pls}}, \code{\link{pathmox}},
#'\code{\link{techmox}}
#' @method plot treemox.pls
#' @S3method plot treemox.pls
#'@examples
#'
#'  \dontrun{
#'  ## example of PLS-PM in customer satisfaction analysis
#'  ## model with seven LVs and reflective indicators
#'  data(csimobile)
#'  
#'  # select manifest variables
#'  data_mobile = csimobile[,8:33]
#'  
#'  # define path matrix (inner model)
#'  IMAG = c(0, 0, 0, 0, 0, 0, 0)
#'  EXPE = c(1, 0, 0, 0, 0, 0, 0)
#'  QUAL = c(0, 1, 0, 0, 0, 0, 0)
#'  VAL = c(0, 1, 1, 0, 0, 0, 0)
#'  SAT = c(1, 1, 1, 1, 0, 0, 0)
#'  COM = c(0, 0, 0, 0, 1, 0, 0)
#'  LOY = c(1, 0, 0, 0, 1, 1, 0)
#'  mob_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, COM, LOY)
#'  
#'  # blocks of indicators (outer model)
#'  mob_blocks = list(1:5, 6:9, 10:15, 16:18, 19:21, 22:24, 25:26)
#'  mob_modes = rep("A", 7)
#'  
#'  # apply plspm
#'  mob_pls = plspm(data_mobile, mob_path, mob_blocks, modes = mob_modes, 
#'                  scheme = "factor", scaled = FALSE)
#'
#'  # re-ordering those segmentation variables with ordinal scale (Age and Education)
#'  csimobile$Education = factor(csimobile$Education, 
#'      levels=c("basic","highschool","university"),
#'      ordered=TRUE)
#'  
#'  # select the segmentation variables
#'  seg_vars = csimobile[,1:7]
#'  
#'  # Pathmox Analysis
#'  mob_pathmox = pathmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)
#'  
#'  # applying function treemox.pls
#'  mob_nodes <- treemox.pls(mob_pls, mob_pathmox)
#'
#'  # default plot
#'  # comparative barplots of endogenous latent variables between nodes
#'  plot(mob_nodes, comp.by="nodes")  
#'
#'  # comparative barplots of nodes between latent regressors
#'  plot(mob_nodes, comp.by="latents", decreasing=TRUE)
#'  }
#'
plot.treemox.pls <-
function(x, comp.by="nodes", nodes.names=NULL, ordered=TRUE, decreasing=FALSE, 
       color=NULL, show.box=TRUE, border=NA, cex.names=.75, cex.axis=.75, short.labs=TRUE, short.min=NULL, ...)
{
    # ARGUMENTS
    # x: an object of class treemox.pls
    # comp.by: one of "nodes" or "latents"
    # nodes.names: optional vector of names for the nodes (length of terminal nodes)
    # ordered: logical value indicating whether bars are ordered (increasing order)
    # decreasing: logical value indicating if the sort order should be increasing or decreasing
    # color: optional vector of colors for the bars
    # show.box: logical value indicating whether a box is drawn around each barplot
    # border: color of the borders
    # cex.names: expansion factor of labels (x-axis)
    # cex.axis: expansion factor of y-axis
    # short.labs: logical value indicating if labels should be abbreviated
    # short.min: integer number indicating the minimum length of the abbreviations

    if (length(union(comp.by,c("nodes","latents")))>2)
        stop("Invalid argument 'comp.by'")
    if (!is.null(nodes.names)) {
        if (length(nodes.names)!=(ncol(x$paths)-1)) {
            cat("NOTICE: Invalid argument 'nodes.names'. Default names were used", "\n")   
        } else {
            colnames(x$paths) <- c("Root_Node",nodes.names)
        }        
    }
    if (short.labs)
        if (mode(short.min)!="numeric" || length(short.min)!=1 || (short.min%%1)!=0)
            short.min <- 7

    IDM <- x$IDM
    lvs <- nrow(IDM)
    idm <- as.vector(IDM)
    idm[idm==1] <- 1:sum(IDM)
    SDM <- matrix(idm,lvs,lvs)
    B <- x$paths[,-1] - x$paths[,1]
    rs <- c(1,1,1,2,2,2,2,2,3,2,3,3,4,4)
    cs <- c(1,2,3,2,3,3,4,4,3,5,4,4,4,4)
    index.mat <- cbind(1:14,rs,cs)
    endo <- rowSums(IDM)
    endo[endo!=0]<-1
    who.endo <- which(endo!=0)

    if (comp.by=="nodes")
    {    
        nn <- ncol(B)
        if (is.null(color))  colors=rainbow(lvs,s=.5,v=.7)  else  colors=rep(color,length.out=lvs)
        for (i in who.endo)
        {
            n.ind <- which(IDM[i,]==1)
            indep <- SDM[i,which(SDM[i,]!=0)]
            arg.names <- colnames(IDM)[n.ind]
            if (short.labs) 
                arg.names <- abbreviate(arg.names, minlength=short.min)
            cols <- colors[n.ind]            
            dev.new()
            par(mfrow=index.mat[nn+1,2:3])
            par(mar=c(3,3,3,3)) 
            barplot(x$paths[indep,1], names.arg=arg.names, main=paste("Global:",rownames(IDM)[i]), col=cols,
                    border=border)
            abline(h=0)
            if (show.box) box("figure")
            ylim <- 1.15 * c(min(B[indep,]),max(B[indep,]))  
            for (n in 1:nn)
            {   
                if (ordered) {
                    barplot(sort(B[indep,n], decreasing), names.arg=arg.names[order(B[indep,n], decreasing=decreasing)], 
                            main=colnames(B)[n], border=border, col=cols[order(B[indep,n], decreasing=decreasing)], 
                            cex.names=cex.names, cex.main=1, cex.axis=cex.axis, ylim=ylim, ...)
                } else {
                     barplot(B[indep,n], names.arg=arg.names, main=colnames(B)[n], border=border, 
                            col=cols, cex.names=cex.names, cex.main=1, cex.axis=cex.axis, ylim=ylim, ...)
                }
                abline(h=0)
                if (show.box) box("figure")
            }   
        }
    } else {
        # comp.by="latents"
        ncols <- ncol(B)
        if (is.null(color))  colors=rainbow(ncols,s=.5,v=.7)  else  colors=rep(color,length.out=ncols)
        arg.names <- colnames(B)
        if (short.labs) 
            arg.names <- abbreviate(arg.names, minlength=short.min)
        for (i in who.endo)
        {
            indep <- SDM[i,which(SDM[i,]!=0)]    
            dev.new()
            par(mfrow=index.mat[length(indep),2:3])
            par(mar=c(3,3,3,3)) 
            ylim <- 1.15 * c(min(B[indep,]),max(B[indep,]))
            for (k in indep)
            {
                if (ordered) {
                    barplot(sort(B[k,], decreasing), names.arg=arg.names[order(B[k,], decreasing=decreasing)], 
                            main=rownames(B)[k], border=NA, cex.main=1,cex.names=cex.names, cex.axis=cex.axis, 
                            ylim=ylim, col=colors[order(B[k,], decreasing=decreasing)], ...)
                } else {
                    barplot(B[k,], names.arg=arg.names, main=rownames(B)[k], border=NA, cex.main=1,
                            cex.names=cex.names, cex.axis=cex.axis, ylim=ylim, col=colors, ...)
                }                  
                abline(h=0)
                if (show.box) box("figure")
            }   
        }
    } 
}

