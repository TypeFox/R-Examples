#' Function to visualise the same graph but with multiple graph node colorings according to input data matrix
#'
#' \code{visNetMul} is supposed to visualise the same graph but with multiple colorings according to input data matrix
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param data an input data matrix used to color-code vertices/nodes. One column corresponds to one graph node coloring. The input matrix must have row names, and these names should include all node names of input graph, i.e. V(g)$name, since there is a mapping operation. After mapping, the length of the patern vector should be the same as the number of nodes of input graph. The way of how to color-code is to map values in the pattern onto the whole colormap (see the next arguments: colormap, ncolors, zlim and colorbar)
#' @param height a numeric value specifying the height of device
#' @param margin margins as units of length 4 or 1
#' @param border.color the border color of each figure
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param colorbar logical to indicate whether to append a colorbar. If pattern is null, it always sets to false
#' @param colorbar.fraction the relative fraction of colorbar block against the figure block
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param glayout either a function or a numeric matrix configuring how the vertices will be placed on the plot. If layout is a function, this function will be called with the graph as the single parameter to determine the actual coordinates. This function can be one of "layout.auto", "layout.random", "layout.circle", "layout.sphere", "layout.fruchterman.reingold", "layout.kamada.kawai", "layout.spring", "layout.reingold.tilford", "layout.fruchterman.reingold.grid", "layout.lgl", "layout.graphopt", "layout.svd" and "layout.norm". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param mtext.side on which side of the mtext plot (1=bottom, 2=left, 3=top, 4=right)
#' @param mtext.adj the adjustment for mtext alignment (0 for left or bottom alignment, 1 for right or top alignment)
#' @param mtext.cex the font size of mtext labels
#' @param mtext.font the font weight of mtext labels
#' @param mtext.col the color of mtext labels
#' @param ... additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.
#' @return 
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{visNet}}, \code{\link{visNetAnimate}}
#' @include visNetMul.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/80)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#'
#' # 3) visualise the module with vertices being color-coded by scores
#' nnodes <- vcount(subg)
#' nsamples <- 10
#' data <- matrix(runif(nnodes*nsamples), nrow=nnodes, ncol=nsamples)
#' rownames(data) <- V(subg)$name
#' visNetMul(g=subg, colormap="bwr", data=data, glayout=layout.fruchterman.reingold)

visNetMul <- function (g, data, height=7, margin=rep(0.1,4), border.color="#EEEEEE", colormap=c("bwr","jet","gbr","wyr","br","yr","rainbow","wb"), ncolors=40, zlim=NULL, colorbar=T, colorbar.fraction=0.25, newpage=T, glayout=layout.fruchterman.reingold, mtext.side=3, mtext.adj=0,mtext.cex=1,mtext.font=2,mtext.col="black", ...)
{   

    ## check input graph
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    if(is.null(V(ig)$name)){
        V(ig)$name <- as.character(V(ig))
    }

    ## check input data
    if(is.matrix(data) | is.data.frame(data) | is.vector(data)){
        data <- as.matrix(data)
    }else if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }

    if(is.null(rownames(data))) {
        stop("The function must require the row names of the input data.\n")
    }else if(any(is.na(rownames(data)))){
        warning("Data with NA as row names will be removed")
        data <- data[!is.na(rownames(data)),]
    }
    cnames <- colnames(data)
    if(is.null(cnames)){
        cnames <- seq(1,ncol(data))
    }
    
    ## check mapping between input data and graph
    ind <- match(rownames(data), V(ig)$name)
    nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
    if(length(nodes_mapped)!=vcount(ig)){
        stop("The function must require that the row names of input data contain all those in the input graph.\n")
    }
    data <- as.matrix(data[nodes_mapped,])
    
    ## determine the color range
    vmin <- floor(stats::quantile(data, 0.05))
    vmax <- ceiling(stats::quantile(data, 0.95))
    if(vmin < 0 & vmax > 0){
        vsym <- abs(min(vmin, vmax))
        vmin <- -1*vsym
        vmax <- vsym
    }
    if(!is.null(zlim)){
        if(zlim[1] < floor(min(data)) | zlim[2] > ceiling(max(data))){
            #zlim <- c(vmin,vmax)
        }
    }else{
        zlim <- c(vmin,vmax)
    }
    
    ######################################################################################
    colNum <- ceiling(sqrt(length(cnames))) + 1
    rowNum <- ceiling((length(cnames)+1)/colNum) + 1
    tolNum <- colNum*rowNum
    if(tolNum < length(cnames)+colNum+rowNum-1){
        rowNum <- rowNum+1
    }
    
    ## a matrix object specifying the location of the next N figures on the output device. Each value in the matrix must be 0 or a positive integer
    layout_vec <- vector()
    k <- 0
    t <- 0
    for(j in 1:rowNum){
        for(i in 1:colNum){
            t <- t+1
            if(j == 1){
                layout_vec[t] <- 0
            }else{
                if(i == colNum){
                    if(j == rowNum){
                        layout_vec[t] <- length(cnames)+1
                    }else{
                        layout_vec[t] <- 0
                    }
                }else{
                    k <- k+1
                    if(k<=length(cnames)){
                        layout_vec[t] <- k
                    }else{
                        layout_vec[t] <- 0
                    }
                }
            }
        }
    }
    layout_matrix <- matrix(layout_vec, rowNum, colNum, byrow=T)
    
    ## relative heights and widths
    frac <- colorbar.fraction
    ## relative heights for each row
    row_first <- frac*10/rowNum
    row_rest <- (10-row_first)/(rowNum-1)
    layout_heights <- c(row_first, rep(row_rest, rowNum-1))
    ## relative widths for each column
    col_last <- frac*10/colNum
    col_rest <- (10-col_last)/(colNum-1)
    layout_widths <- c(rep(col_rest, colNum-1), col_last)
    
    ######################################################################################
    if (newpage){
        grDevices::dev.new(width=height*colNum/rowNum, height=height)
    }
    graphics::par(mfrow=c(rowNum,colNum), mar=margin)
    graphics::layout(layout_matrix, widths=layout_widths, heights=layout_heights)
    if(is.function(glayout)){
        glayout_fix <- glayout(ig)
    }else{
        glayout_fix <- glayout
    }
    for(k in 1:length(cnames)){
        visNet(ig, glayout=glayout_fix, pattern=data[,k], colormap=colormap, ncolors=ncolors, zlim=zlim, colorbar=F, newpage=F, ...)
        graphics::mtext(sprintf("%s",cnames[k]), line=-1.5, side=mtext.side, adj=mtext.adj, cex=mtext.cex, font=mtext.font, col=mtext.col)
        graphics::box("figure",col=border.color)
    }
    #graphics::box("outer", col="black", lwd=4)
    
    ######################################################################################
    ## colorbar
    if(colorbar){
        
        ## empty plot
        plot(c(0,1),c(0,1),xlab="", ylab="", axes=F, type="n")
            
        palette.name <- visColormap(colormap=colormap)
        colors <- palette.name(ncolors)
        lab.scale <- length(colors)/(zlim[2]-zlim[1])
        for (i in 1:length(colors)) {
            yValue <- (i-1)/ncolors
            hValue <- 1/ncolors
            xValue <- 0.25
            wValue <- 0.25
        
            ## for rect
            xleft <- xValue
            ybottom <- yValue
            xright <- xValue+wValue
            ytop <- yValue+hValue
            graphics::rect(xleft,ybottom,xright,ytop, col=colors[i], border="transparent")
        
            if(i == 1 | i == 1+length(colors)/2){
                tx <- (i-1)/lab.scale + zlim[1]
                graphics::text(x=xright+0.1, y=ybottom, labels=tx, cex=1)
            }else if(i==length(colors)){
                tx <- i/lab.scale + zlim[1]
                graphics::text(x=xright+0.2, y=ytop, labels=tx, cex=1)
            }
        }
    }
    
    invisible()
}