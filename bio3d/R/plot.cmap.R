plot.cmap <- function(x, col=2, pch=16, 
                      main="Contact map", sub="",
                      xlim=NULL, ylim=NULL, 
                      xlab = "Residue index", ylab = xlab, 
                      axes=TRUE, ann=par("ann"),
                      sse=NULL, sse.type="classic", sse.min.length=5,
                      bot=TRUE, left=TRUE,
                      helix.col="gray20", sheet.col="gray80",
                      sse.border=FALSE, add=FALSE,
                      ...) {
  dims <- dim(x)

  if(is.null(xlim))
    xlim <- c(1, dims[1])
  if(is.null(ylim))
    ylim <- c(1, dims[2])

  if(!add) {
    plot.new()
  }
  else {
    axes <- FALSE
    xlab <- NA; ylab <- NA;
    main <- NA; sub <- NA;
    sse <- NULL;
  }
  plot.window(xlim=xlim, ylim=ylim, ...)
  
  inds <- which(x==1, arr.ind=TRUE)
  points(inds, pch=pch, col=col)
  
  if(!is.null(sse)) {
    ## Obtain SSE vector from PDB input 'sse'
    if(is.pdb(sse)) 
      sse$sse <- pdb2sse(sse)
    
    h <- bounds( which(sse$sse == "H") )
    e <- bounds( which(sse$sse == "E") )

    ## Remove short h and e elements that can crowd plots
    if(length(h) > 0) {
      inds <- which(h[,"length"]>=sse.min.length)
      h <- h[inds,,drop=FALSE]
    } else { h <- NULL }
    if(length(e) > 0) {
      inds <- which(e[,"length"]>=sse.min.length)
      e <- e[inds,,drop=FALSE]
    } else { e <- NULL }
    
    if(sse.type != "classic")
      warning("Only sse.type='classic' is currently available, 'fancy' coming soon")

    off <- c(0.006, 0.039)
    if(left) {
      ## Determine bottom and top of margin region 
      bo <- min(xlim) - (diff(xlim)*off[1]) 
      to <- min(xlim) - (diff(xlim)*off[2]) 

      
      if(length(h) > 0)
        rect(xleft=bo, xright=to, 
             ybottom=h[,"start"], ytop=h[,"end"],
             col=helix.col, border=sse.border)
      
      if(length(e) > 0)
        rect(xleft=bo, xright=to, 
             ybottom=e[,"start"], ytop=e[,"end"],           
             col=sheet.col, border=sse.border)
    }
    if(bot){
      	to <- min(ylim) - (diff(ylim)*off[1])
      	bo <- min(ylim) - (diff(ylim)*off[2])
        
      	if(length(h) > 0)
          rect(xleft=h[,"start"], xright=h[,"end"],
               ybottom=bo, ytop=to, col=helix.col, border=sse.border)
        
      	if(length(e) > 0)
          rect(xleft=e[,"start"], xright=e[,"end"],
               ybottom=bo, ytop=to, col=sheet.col, border=sse.border)
    }
  }

  if(axes) {
    box()
    at <- axTicks(1); at[1] = 1
    axis(1, at)
    axis(2, at)
  }
  
  if(ann) {
#    if(is.null(xlab))  xlab=xy$xlab
#    if(is.null(ylab))  ylab=xy$ylab
    title(main=main, sub=sub, 
          xlab=xlab, ylab=ylab, ...)
  }
  
}
