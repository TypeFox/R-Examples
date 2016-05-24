"plot.fasta" <- function(x, plot.labels=TRUE,
                         plot.bars=TRUE,
                         plot.lines=FALSE,
                         plot.axis=TRUE,
                         seq.index=NULL,
                         color.conserved=FALSE,
                         cutoff=0.5,
                         col=NULL,
                         bars.scale=2,
                         row.spacing=0.5, 
                         aln.col="grey50",
                         cex.text=1, add=FALSE, ... ) {

  if(!(inherits(x, "fasta") | inherits(x, "pdbs")))
    stop("input 'x' should be a list object as obtained from 'read.fasta'")

  row.height <- 1
  
  ## calculate sequence entropy
  ##ent <- entropy(x)
  ##hn <- ent$H.norm
  hn <- conserv(x$ali, method="entropy10")

  ## check for gaps
  gaps <- gap.inspect(x$ali)
  dims <- dim(x$ali)
  
  ## set x- and y-lim
  if(plot.labels)
    xlim <- c(1, dims[2])
  else
    xlim <- c(1, dims[2])

  if(plot.bars) {
    conserv.height <- row.height * bars.scale
    ylim <- c(0, (dims[1]*row.height) + conserv.height)
  }
  else {
    conserv.height <- NULL
    ylim <- c(0, dims[1]*row.height)
  }

  if(is.null(col)) {
    num.cols <- length(seq(cutoff, 1, by=0.1))-1
    col <- heat.colors(num.cols)
  }
  else
    num.cols <- length(col)

  ## start the plot
  if(!add) {
    plot.new()
    plot.window(xlim, ylim)
  }

  ##a: space between row bars, b: height of grey row bars
  b <- 1 / (row.height + row.spacing)
  a <- b * row.spacing
  
  ## plot the sequence 'boxes' so that residue i is plotted at
  ## in coordinates [i,i+1)

  all.start <- NULL; all.end <- NULL;
  for ( i in 1:nrow(x$ali) ) {
    nongap.inds <- which(gaps$bin[i,]==0)
    bs <- bounds(nongap.inds)
    
    for ( j in 1:nrow(bs) ) {
      xstart <- bs[j,"start"]
      xend <- bs[j,"end"]+1
           
      ystart <- ((i-1) * row.height) + 0.5
      yend <- ystart + (row.height-a)

      rect(xstart, ystart, xend, yend,
           col=aln.col, border=NA)
    }

    all.start <- c(all.start, ystart)
    all.end <- c(all.end, yend)
    
    ## label on the side of the sequence
    if(plot.labels)
      mtext(x$id[i], side=4, line=-2, at=ystart+(b*.5),
            las=1, cex=cex.text)
  }


  aa2col <- function(aa, default="grey50") {
    if(aa %in% c("A", "V", "L", "I", "M", "F", "W"))
      return("dodgerblue3")
    if(aa %in% c("S", "T", "N", "Q"))
      return("forestgreen")
    if(aa %in% c("D", "E"))
      return("purple")
    if(aa %in% c("R", "K"))
      return("red")
    if(aa %in% c("H", "Y"))
      return("cyan")
    if(aa %in% c("G"))
      return("orange")
    if(aa %in% c("P"))
      return("yellow")
    if(aa %in% c("C"))
      return("darksalmon")
    else {
      warning(paste("color for amino acid", aa, "not found"))
      return(default)
    }
  }

  if(color.conserved) {
    ## which columns to color
    inds <- intersect(which(hn>=cutoff), which(gaps$col<(nrow(gaps$bin)/2)))
  
    for ( j in 1:length(inds) ) {
      col.ind <- inds[j]
      
      ## only non-gap positions
      nongap.inds <- which( gaps$bin[, col.ind]==0 )
      
      if(length(nongap.inds)>0) {
        for (k in 1:length(nongap.inds) ) {
          row.ind <- nongap.inds[k]

          aa <- x$ali[row.ind, col.ind]
          #if(aa=="-")
          #  print(aa)
          
          aacol <- aa2col(aa, default=aln.col)

          xstart <- inds[j]
          xend <- xstart+1
          ystart <- ((row.ind-1) * row.height) + 0.5
          yend <- ystart + (row.height-a)
          
          rect(xstart, ystart, xend, yend, 
               col=aacol, border=NA)
        }
      }
    }
  }
  

  if(plot.bars || plot.lines)
    i=i+1
  
  ## plot conservation bar
  if(plot.bars) {
    ints <- rev(seq(0, 1, by=0.1))
    ints = ints[ints >= cutoff]
    
    for ( j in 1:(length(ints)-1) ) {
      cons.inds <- intersect(
        intersect(which(hn <= ints[j]),
                  which(hn > ints[j+1])),      
        gaps$f.inds)

      if(length(cons.inds)>0) {
        for (k in 1:length(cons.inds) ) {
          xstart <- cons.inds[k]
          xend <- xstart+1
          
          ystart <- ((i-1) * row.height) + 0.5
          yend <- ystart + (conserv.height*hn[cons.inds[k]])
          
          ##print( paste(ystart, yend, hn[cons.inds[k]] ))
          
          rect(xstart, ystart, xend, yend,
               col=col[j], border=NA)
        }
      }
      
    }
  }

  if(plot.lines) {
    ystart <- ((i-1) * row.height) + 0.5
    lines(1:dims[2]+0.5, ystart+(conserv.height*hn), type="o", cex=0.4)
  }

  if(plot.labels & (plot.lines | plot.bars)) {
    mtext("Conservation", side=4, line=-2., at=ystart+(b*.5),
          las=1, cex=cex.text)
  }
 

  ticks <- bounds(gaps$t.inds)
  starts <- c(1, ticks[,"start"])
  ends <- c(ticks[,"end"], dims[2])
  at <- c(starts, ends+1)
  
  if(plot.axis) {
    if(!is.null(seq.index)) {
      ig <- is.gap(x$ali[seq.index,])
      tick.labs <- rep(NA, length=dims[2])
      tick.labs[which(!ig)] = seq(1, length(which(!ig)))
      labs <- c(tick.labs[ starts ], tick.labs[ ends ]+1)
      labs[ 1 ] = 1
    }
    else {
      tick.labs <- seq(1, dims[2])
      labs <- c(tick.labs[ starts ], tick.labs[ ends ]+1)
    }
    
    if(!is.null(seq.index))
      axis(1, at=at, labels=labs, tick=FALSE, line=-2.)

    axis(1, line=0.5, at=at)
    mtext("Alignment index", 1, line=2, cex=1.0)
  }

  out <- list(at=sort(unique(at)), start=all.start, end=all.end)
  
  invisible(out)
}
