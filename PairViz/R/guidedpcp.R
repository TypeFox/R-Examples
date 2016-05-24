pcp <- 
function (data, order = NULL, panel.colors = NULL, col = 1, lty = 1, 
    horizontal = TRUE, mar = NULL, scale=TRUE,axis.width=0,axis.grid.col="grey70",connect=TRUE,...) 
{
    if (is.null(mar)) 
        if (horizontal == TRUE) 
            mar <- c(2, 2, 2, 2) + 0.1
        else mar <- c(2, 5, 2, 2) + 0.1
    par("mar"=mar)
    if (!is.null(order)) {
    	axis.labels <- colnames(data)[order]
        data <- data[, order]
        if (is.matrix(panel.colors)) 
            panel.colors <- panel.colors[order, order]
            }
    else axis.labels <- colnames(data)
    if (is.matrix(panel.colors)) 
        panel.colors <- panel.colors[col(panel.colors) == row(panel.colors) + 1]
    if (is.vector(panel.colors)) 
        if (ncol(data) - 1 != length(panel.colors)) 
            stop("dimensions do not match")
    if (scale==TRUE) x <- apply(data, 2, function(x) (x - min(x))/(max(x) - min(x)))
    else x <- data
    p <- ncol(x)
     indx <- 1:p
    if ((length(axis.width)==1) && (axis.width == 0))
      bx <- x
    else {
      bx <- x[,rep(1:p,times=rep(2,times=p))]
      if (length(axis.width)==1)
        indx <- as.vector(sapply(indx,function(x)c(x-axis.width/2,x+axis.width/2)))
      else indx <- as.vector(t(cbind(indx-axis.width/2,indx+axis.width/2)))
    	}
     linesr <- range(bx)
     pts <- rep(1,length.out=p)
     if (horizontal == TRUE) {
        matplot(indx, t(bx), xlab = "", ylab = "", axes = FALSE, 
            type = "n", xaxs="i",pch=pts,...)
        axis(1, at = 1:p, labels = axis.labels)
        if (!(is.null(panel.colors))) 
            for (i in 1:(p - 1)) rect(i, 0, i + 1, 1, lty = 0, 
                col = panel.colors[i])
         if (!is.null(axis.grid.col))
        for (i in 1:p) lines(c(i, i), linesr, col = axis.grid.col)
        if (connect)
        matpoints(indx, t(bx), type = "l", col = col, lty = lty, pch=pts,...)
        else if (axis.width == 0)
        matpoints(indx, t(bx),  col = col,  pch=pts,...)
         else for (i in seq(1,length(indx),2))
                  matpoints(indx[i:(i+1)], t(bx[,i:(i+1)]), type = "l", col = col, lty = lty, pch=pts,...)
         
            
            }
    else {
        matplot(t(bx), rev(indx), xlab = "", ylab = "", axes = FALSE, 
            type = "n",yaxs="i",pch=pts, ...)
        axis(2, at = p:1, labels = axis.labels, las = 2)
        if (!(is.null(panel.colors))) 
            for (i in 1:(p - 1)) rect(0, i, 1, i + 1, lty = 0, 
                col = panel.colors[p - i])
        if (!is.null(axis.grid.col))
          for (i in 1:p) lines(linesr, c(i, i), col = axis.grid.col)
        if (connect)
        matpoints(t(bx), rev(indx), type = "l", col = col, lty = lty, pch=pts, ...)
        else if (axis.width == 0)
        matpoints(t(bx), rev(indx),  col = col,  pch=pts,...)
        else for (i in seq(1,length(indx),2))
                matpoints(t(bx[,i:(i+1)]), rev(indx[i:(i+1)]), type = "l", col = col, lty = lty, pch=pts, ...)
            
            
    }
    invisible()
}


    
    
    
guided_pcp <- function(data, edgew=NULL, path = NULL, pathw=NULL,zoom=NULL,pcpfn=pcp,
     pcp.col = 1,lwd=0.5, panel.colors=NULL,  pcp.mar=c(1.5,2,2,2), pcp.scale=TRUE, bar.col=1:9,bar.axes=FALSE, bar.mar=NULL, 
     bar.ylim=NULL,reorder.weights=TRUE,
    layout.heights=NULL, layout.widths=c(10,1),
     main=NULL,legend=FALSE,cex.legend = 1,legend.mar=c(1,4,1,1),...){
    
  	if (is.null(path)) path <- 1:ncol(data)
	if (is.function(path)) 
	   o <- find_path(-edgew,path,...)
	else o <- path
    if (is.null(pathw) || (nrow(as.matrix(pathw)) != length(o) -1))
       if (!is.null(edgew)) pathw <- path_weights(edgew,o,...)
       
    pathw <- as.matrix(pathw)
    if (!is.null(zoom)) {
    	o <- o[zoom]
        pathw <- as.matrix(pathw[zoom,])
        pathw <- as.matrix(pathw[-nrow(pathw),])  
       	}
    if (reorder.weights==TRUE)
       w.ord <- order(apply(pathw,2,sum),decreasing=TRUE)
    else w.ord <- 1:ncol(pathw)
    
    if (legend ==TRUE){
       if (is.null(layout.heights)) 
        	layout.heights <- c(5,5,3)
        layout(matrix(c(1,1,2,4,3,3),nrow=3), heights=layout.heights, widths=layout.widths)
        }
     else {
       if (is.null(layout.heights)) layout.heights <- c(10,3)
     layout(matrix(c(1,2)), heights=layout.heights )
     }
     
    pcpfn(data,col=pcp.col,order=o,panel.colors=panel.colors,
    horizontal=TRUE,lwd=lwd,mar=pcp.mar,main=main,xaxs="i",scale=pcp.scale,...)
    if (is.null(bar.mar)) {
    	bar.mar <- pcp.mar
    	bar.mar[c(1,3)] <- bar.mar[c(1,3)]/3   	
    	}
    par("mar"=bar.mar)
     if (length(bar.col) == length(w.ord))
      bar.col <- bar.col[w.ord]
    p1 <- par("usr")
    p2 <- p1-1
    p2[3] <- min(apply(pathw,1,function(r) sum(pmin(r,0))))
    p2[4] <- max(apply(pathw,1,function(r) sum(pmax(r,0))))
    if (!is.null(bar.ylim)) p2[3:4] <- bar.ylim

    plot.new()
    par("usr"=p2)
    if (ncol(pathw)==1)
    barplot(pathw[,1],space=0,axisnames= FALSE,
        col=bar.col,axes=bar.axes,main="",xaxs="i",add=TRUE)
   else
    barplot(t(pathw[,w.ord]),space=0,axisnames= FALSE,
        col=bar.col,axes=bar.axes,main="",xaxs="i",add=TRUE)
     if (legend==TRUE) {
         par("mar"=legend.mar)
       barplot(rep(1,length(w.ord)),col=bar.col,horiz=TRUE,space=0, axes=FALSE,
       names.arg=colnames(pathw[,w.ord]),las=2,cex.names=cex.legend)
       }


     invisible()
    }

