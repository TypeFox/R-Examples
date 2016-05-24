plotb3 <- function(x, resno=NULL, rm.gaps = FALSE, type="h",
                       main="", sub="",
                       xlim=NULL, ylim=NULL, ylim2zero=TRUE,
                       xlab = "Residue", ylab = NULL, 
                       axes=TRUE, ann=par("ann"),
                       col=par("col"),
                       sse=NULL, sse.type="classic", sse.min.length=5,
                       top=TRUE, bot=TRUE,
                       helix.col="gray20", sheet.col="gray80",
                       sse.border=FALSE, 
                       ...) {
   
  ## Check for gap positions 
  gaps.pos = gap.inspect(x)
  if(is.matrix(x)) x = x[1, ]   ## should support matrix in future

  if(!is.null(resno)) {
    if(is.pdb(resno)) {
      ## Take Calpha residue numbers from PDB input
      ca.inds <- atom.select(resno, "calpha", verbose = FALSE)
      resno <- resno$atom$resno[ca.inds$atom]
    }
    if(any(is.na(x))) {
       tmp.resno <- rep(NA, length(x))
       tmp.resno[!is.na(x)] = resno
       resno = tmp.resno
    }
    if(length(resno) != length(x)) {
      warning("Length of input 'resno' does not equal the length of input 'x'; Ignoring 'resno'")
      resno=NULL
    }
  }
  
  if(rm.gaps) {
     xy <- xy.coords(x[gaps.pos$f.inds])
     if(!is.null(resno)) resno <- resno[gaps.pos$f.inds]
  }
  else 
     xy <- xy.coords(x)

  if (is.null(xlim))
    xlim <- range(xy$x[is.finite(xy$x)])
  if (is.null(ylim))
    ylim <- range(xy$y[is.finite(xy$y)])
  if(ylim2zero) ylim[1]=0
  
  plot.new()
  plot.window(xlim, ylim, ...)
  points(xy$x, xy$y, col=col, type=type, ...)
  
  if(!is.null(sse)) {	

    ## Obtain SSE vector from PDB input 'sse'
    if(is.pdb(sse)) 
      sse$sse <- pdb2sse(sse)
 
    h <- bounds( which(sse$sse == "H") )
    e <- bounds( which(sse$sse == "E") )
    
    ## Remove short h and e elements that can crowd plots
    if(length(h) > 0) {
      h <- h[h[,"length"] >= sse.min.length,,drop=FALSE]
    } else { h <- NULL }
    if(length(e) > 0) {
      e <- e[e[,"length"] >= sse.min.length,,drop=FALSE]
    } else { e <- NULL }

    ## For gaps
    if(length(gaps.pos$t.inds) > 0) {

       # unwrap SSE after length filtering
       tmp.sse = rep(" ", length(x))
       tmp.inds = which(!is.na(x))
       if(length(h) > 0) tmp.sse[tmp.inds[unbound(h)]] = "H"
       if(length(e) > 0) tmp.sse[tmp.inds[unbound(e)]] = "E"
      
       # remove gaps if required
       if(rm.gaps) tmp.sse = tmp.sse[gaps.pos$f.inds] 

       # new SSE segments
       h <- bounds( which(tmp.sse == "H") )
       e <- bounds( which(tmp.sse == "E") )

    }

    if(sse.type != "classic")
      warning("Only sse.type='classic' is currently available, 'fancy' coming soon")
    
    if(top) {
      ## Determine bottom and top of margin region 
      bo <- max(ylim) + (diff(ylim)*0.001) # 0.1% 
      to <- max(ylim) + (diff(ylim)*0.04) # 4%
      
      if(length(h) > 0)
        rect(xleft=h[,"start"], xright=h[,"end"],
             ybottom=bo, ytop=to, col=helix.col, border=sse.border)
      
      if(length(e) > 0)
        rect(xleft=e[,"start"], xright=e[,"end"],
             ybottom=bo, ytop=to, col=sheet.col, border=sse.border)
    }
    if(bot){
      	to <- min(ylim) - (diff(ylim)*0.001)
      	bo <- min(ylim) - (diff(ylim)*0.04)
        
      	if(length(h) > 0)
          rect(xleft=h[,"start"], xright=h[,"end"],
               ybottom=bo, ytop=to, col=helix.col, border=sse.border)
        
      	if(length(e) > 0)
          rect(xleft=e[,"start"], xright=e[,"end"],
               ybottom=bo, ytop=to, col=sheet.col, border=sse.border)
    }
  }

  if(axes) {
    axis(2)
    box()
    at <- axTicks(1); at[1] = 1
    if(is.null(resno)) {
      axis(1, at)
    } else {
      labels <- resno[at]
      labels[is.na(labels)] <- "" # for gaps, no label
      axis(1, at=at, labels=labels)
    }
  }
  if(ann) {
    if(is.null(xlab))  xlab=xy$xlab
    if(is.null(ylab))  ylab=xy$ylab
    title(main=main, sub=sub, 
          xlab=xlab, ylab=ylab, ...)
  }
}

plot.bio3d <- function(...) { plotb3(...) }
