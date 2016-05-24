"hclustplot" <- function(hc, k=NULL, h=NULL, colors=NULL,
                           labels=NULL, fillbox=FALSE, 
                           heights = c(1, .3), mar = c(1, 1, 0, 1), ...) {

  if(!inherits(hc, "hclust"))
    stop("hc must be of type 'hclust'")
  if(is.null(k) & is.null(h) & is.null(colors))
    stop("provide either k or h to function 'cutree', or colors for manual coloring")
  
  mtext.names <- names(formals( mtext ))
  plot.dendrogram <- get("plot.dendrogram", envir = getNamespace("stats"))
  plot.names <- c(names(formals( plot.dendrogram )),
                  names(formals( plot.default )))
  
  dots <- list(...)
  mtext.args <- dots[names(dots) %in% mtext.names]
  plot.args <- dots[names(dots) %in% plot.names]
  par.args <- dots[!(names(dots) %in% unique(c(names(mtext.args), names(plot.args))))]

  mtext.args <- c(mtext.args, par.args)
  plot.args <- c(plot.args, par.args)

  ## set default and allowed mtext arguments
  if(!any(names(mtext.args)=="line")) {
    if( fillbox )
      mtext.args$line <- 0.5
    else
      mtext.args$line <- -0.25
  }
  if(!any(names(mtext.args)=="side"))
    mtext.args$side <- 1
  if(!any(names(mtext.args)=="las"))
    mtext.args$las <- 2
  if(any(names(mtext.args)=="col"))
    mtext.args$col <- NULL
  if(any(names(mtext.args)=="at"))
    mtext.args$at <- NULL

  ## set default and allowed plot.dendrogram arguments
  if(any(names(plot.args)=="axes")) {
    axes <- plot.args$axes
    plot.args$axes <- NULL
  }
  else {
    axes <- TRUE
  }
  if(any(names(plot.args)=="xaxs"))
    plot.args$xaxs <- NULL
  if(any(names(plot.args)=="leaflab"))
    plot.args$leaflab <- NULL
  if(any(names(plot.args)=="xlab"))
    plot.args$xlab <- NULL
  if(any(names(plot.args)=="horiz"))
    plot.args$horiz <- NULL
  
  ## print(mtext.args)
  ## print(par.args)
  ## print(plot.args)

  plot.labels <- TRUE
  if(is.logical(labels)) {
    if(labels)
      plot.labels <- TRUE
    else
      plot.labels <- FALSE
    labels <- NULL
  }
    
  if(is.null(labels)) {
    labels <- hc$labels
    if(is.null(labels))
      labels <- seq(1, length(hc$order))
  }
  else {
    if( length(hc$order) != length(labels) )
      stop("labels must be of same length as hc")
  }
  
  if(!is.null(colors)) {
    unq.cols <- unique(colors)
    grps <- unlist(lapply(colors, function(x) which(x==unq.cols)))
    labelColors <- unq.cols
  }
  else {
    grps <- cutree(hc, k=k, h=h)
    labelColors <- seq(1, length(unique(grps)))
  }
  	
  hcd         <- as.dendrogram(hc)
  cols <- labelColors[grps][hc$order]

  ## set margins
  mar.default <- c(1, 1, 0, 1)
  if(all(mar==mar.default)) {
    mar <- mar.default
    mar[1] <- mar[1] + ifelse(plot.labels, 3, 0)
    mar[1] <- mar[1] + ifelse(fillbox, 2, 0)
    mar[2] <- mar[2] + ifelse(!is.null(plot.args$ylab), 2, 0)
    mar[2] <- mar[2] + ifelse(axes, 2, 0)
    mar[3] <- ifelse(!is.null(plot.args$ylab), 4, 2)
  }

  ## colored filled boxes below the dendrogram
  if( fillbox ) { ##| plot.labels ) { 
    layout(as.matrix(c(2,1)), heights = heights)
  
    dev.hold()
    on.exit(dev.flush())
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    par(mar=c(mar[1], mar[2], 0, mar[4]))
    
    if( fillbox )
      image(cbind(1:length(grps)), col = cols, axes = FALSE)
    else
      frame()

    
    if(plot.labels) {
      do.call('mtext', c(list(text=labels[ hc$order ],
                              at=seq(0, 1, length.out=length(grps)),
                              col=cols),
                         mtext.args))
    }
  }

  else {
    layout(1)
  }

  ## dendrogram
  par(mar = c(ifelse(fillbox, 0, mar[1]),
        mar[2], mar[3], mar[4]))
  do.call('plot',  c(list(x=hcd, axes=FALSE, leaflab="none", xaxs="i"),
                     plot.args))
  if(axes)
    axis(2)

  ## labels when filled boxes are not drawn
  if(plot.labels & !fillbox ) {
    do.call('mtext', c(list(text=labels[ hc$order ],
                            at=seq(1, length(grps)),
                            col=cols),
                       mtext.args))
  }

  ##if(!is.null(sub)) {
  ##  mtext(sub, side=3, line=-0.5)
  ##}
  
}
