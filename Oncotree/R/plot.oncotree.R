"plot.oncotree" <-
function(x, edge.weights=c("none","observed","estimated"), edge.digits=2, node.coords=NULL, plot=TRUE, cex=par("cex"), 
col.edge=par("col"), col.text=par("col"), col.weight=par("col"), ...) {
  edge.weights <- match.arg(edge.weights)
  if ((edge.weights=="estimated") & is.null(x$parent$est.weight))
  	stop("No estimated weights are available in the tree.")
  	
  #calculate node widths based in descendent's widths
    get.all.widths <- function(PA){
        get.node.width <- function(node,ww){
            children <- PA$child[PA$parent==node]
            ww2 <- ww
            if (length(children)==0) #node is a leaf
                {ww2[node] <- 1}
            else
                {for (ch in children)
                    {ww2 <- get.node.width(ch,ww2)
                    ww2[node] <- ww2[node] +ww2[ch]}
                }
            return(ww2)
            }
        w <- numeric()
        w[PA$child]<- 0 #initalize vector
        w <- get.node.width('Root', w)
        return(w)
				}
				
            
     get.node.coord <- function(PA, n, param, nodewidth, nodeorder){
        param2 <- param
        node <- nodeorder[1,n]
        nn <- as.numeric(nodeorder["nodenum", nodeorder["nodename",]==node])
        children <- PA$child[PA$parent==node]
        marker <- 0
        if (length(children)>0){ #node is not a leaf
	        for (ch in children){
	            ch.num <- as.numeric(nodeorder["nodenum", nodeorder["nodename",]==ch])
	            if (marker==0){
	              param2$mini[ch.num] <- param2$mini[nn]}
	            else {
	              param2$mini[ch.num] <- marker}
	            param2$maxi[ch.num] <- param2$mini[ch.num]+(param2$unitw[nn]*nodewidth[ch])
	            param2$unitw[ch.num] <- (param2$maxi[ch.num]-param2$mini[ch.num])/nodewidth[ch]
	            marker <- param2$maxi[ch.num]
            }
        }
        return(param2)
        }
          
    PA <- x$parent                  
    nmut <- x$nmut
	  nodeorder <- array(0,dim=c(3,nmut)) 
	  rownames(nodeorder) <- c("nodename", "nodenum", "Y")
    if (is.null(node.coords)) {
			plotinfo <- build.plot(PA)         
	    nodewidth <- get.all.widths(PA)
	    n <- 1
	#    nodewidth.ref <- as.matrix(nodewidth)
	    for (i in 1:nrow(plotinfo$levelgrp)) {
	        for (j in 1:ncol(plotinfo$levelgrp)) {
	#            name <- plotinfo$levelgrp[i,j]
	            if (plotinfo$levelgrp[i,j] != "0") {
	                nodeorder[1,n] <- plotinfo$levelgrp[i,j]
	                # node[3,] are the y coords
	                nodeorder[3,n] <- (nrow(plotinfo$levelgrp)-i+1)
	                nodeorder[2,n] <- n
	                n <- n+1}
	            }
	        }
        
	    mini <- rep(0,nmut)
	    maxi <- rep(0,nmut)
	    unitw <- rep(0,nmut)
	    mini[1] <- 0
	    maxi[1] <- 1
	    unitw[1] <- 1/nodewidth["Root"]
	    param <- list(mini=mini,maxi=maxi,unitw=unitw)
	    for (n in 1:nmut){
	        param <- get.node.coord(PA, n, param, nodewidth, nodeorder) }  
	 
	    X <- (param$mini + param$maxi)/2
	    nodeinfo <- rbind(nodeorder, X)
      xx <- as.numeric(nodeinfo[4,])
      yy <- as.numeric(nodeinfo[3,])
    }
    else {
      if (!all(dim(node.coords)==c(2,nmut)) )
			   stop("'node.coords' has incorrect dimensions")
      if (!all(colnames(node.coords) %in% PA$child))
         stop("the column names of 'node.coords' should be names of nodes of 'x'")
      nodeorder["nodename",] <- colnames(node.coords)
      nodeorder["nodenum",] <- 1:nmut
      xx <- node.coords[1,]
      yy <- node.coords[2,]    
      nodeorder["Y",] <- yy
		}
    
    if (plot){
	    pad.x <- xinch(2/72)  #2pt
	    pad.y <- yinch(2/72)  #2pt
	    plot(xx,yy, axes=FALSE, xlab="", ylab="", type='n', xlim=c(0,1), 
			     ylim=c(min(yy)-0.7, max(yy)+0.7),...)
			for (n in 1:nmut){
	        node <- nodeorder[1,n]
	        children <- PA$child[PA$parent==node]
	        if (length(children)==0) #node is a leaf
	        {}
	        else{
	            for (ch in children){
	            ch.num <- as.numeric(nodeorder["nodenum", nodeorder["nodename",]==ch])
	            lines(c(xx[n],xx[ch.num]), c(yy[n], yy[ch.num]), col=col.edge)
	            if (edge.weights!="none"){
		            ch.PA.num <- which(PA$child==ch)
		            ww <- if(edge.weights=="observed")
		            		   PA$obs.weight[ch.PA.num] 
		            	  else 
		            	     PA$est.weight[ch.PA.num]
		            mid.x <- (xx[n]+xx[ch.num])/2
		            mid.y <- (yy[n]+yy[ch.num])/2
		            w <- formatC(ww,edge.digits, flag="#")
		            rect(mid.x-0.7*strwidth(w, cex=cex)-pad.x,mid.y-0.7*strheight(w, cex=cex)-pad.y,
		                 mid.x+0.7*strwidth(w, cex=cex)+pad.x,mid.y+0.7*strheight(w, cex=cex)+pad.y,
		                 col="white", border=NA)
		            text(mid.x, mid.y ,w, cex=cex, col=col.weight )
	            }
	        }
	    }
	    }
	    vx <- nodeorder[1,]
	    rect(xx-0.7*strwidth(vx, cex=cex)-pad.x,yy-0.7*strheight(vx, cex=cex)-pad.y,
	         xx+0.7*strwidth(vx, cex=cex)+pad.x,yy+0.7*strheight(vx, cex=cex)+pad.y,
		                 col="white")
	    text(xx,yy, vx,  cex=cex, col=col.text)
    }
    res <- rbind(x=xx, y=yy)
    colnames(res) <- nodeorder[1,]
    invisible(res)
}

