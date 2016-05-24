"plot.enma" <-
  function(x, 
           pdbs=NULL, conservation=NULL, variance=FALSE,
           spread = FALSE, offset = 1, 
           col=NULL, signif=FALSE,
           pcut=0.005, qcut=0.04,
           xlab="Alignment Position",
           ylab=c("Fluctuations", "Fluct.variance", "Seq.conservation"),
           xlim=NULL, ylim=NULL,
           mar = c(4, 5, 2, 2),
           ...) {

    if(!(inherits(x, "enma") | inherits(x, "matrix")))
      stop("provide a enma object as obtained from 'nma.pdbs'")
    
    if(spread & (variance | !is.null(conservation))) {
      warning(paste("incompatible arguments:", "\n",
                    " when 'spread=TRUE' conservation and variance will not be plotted"))
      variance <- FALSE
      conservation <- FALSE
    }
    
    ## configure what to plot
    if(inherits(x, "enma"))
      yval <- x$fluctuations
    else
      yval <- x

    ## indices to plot
    if(!is.null(col) && any(is.na(col)))
      row.inds <- which(!is.na(col))
    else
      row.inds <- 1:nrow(yval)

    ## indices for none-"all NA" columns
    gaps.tmp <- gap.inspect(yval[row.inds,,drop=FALSE])
    col.inds <- which(gaps.tmp$col < length(row.inds))
    
    ## colors
    if(is.null(col))
      col <- seq(1, nrow(yval))

    ## full dimensions of yval
    dims.full <- dim(yval)
    
    ## check for gaps
    gaps <- gap.inspect(yval)
    if(any(gaps$col>0))
      rm.gaps <- FALSE
    else
      rm.gaps <- TRUE

    
    ## check if pdbs match enma object
    gaps.pdbs <- NULL
    if(!is.null(pdbs)) {
      if(!inherits(pdbs, "pdbs")) {
        warning("argument 'pdbs' is not a 'pdbs' object (as obtained from pdbaln())")
        pdbs <- NULL
      }
      else {
        gaps.pdbs <- gap.inspect(pdbs$ali)
        
        if(rm.gaps)
          dims.pdbs <- dim(pdbs$ali[, gaps.pdbs$f.inds, drop=FALSE])
        else
          dims.pdbs <- dim(pdbs$ali)
        
        if(!identical(dims.full, dims.pdbs)) {
          warning("dimenension mismatch between modes and pdbs object")
          pdbs <- NULL
        }
      }
    }

    ## reduce all objects to match what we plot
    yval <- yval[row.inds, col.inds, drop=FALSE]
    if(!is.null(pdbs)) {
      if(rm.gaps)
        pdbs <- trim.pdbs(pdbs, row.inds=row.inds, col.inds=gaps.pdbs$f.inds)
      else
        pdbs <- trim.pdbs(pdbs, row.inds=row.inds, col.inds=col.inds)
    }
    col <- col[!is.na(col)]
    
    ## Sequence conservation
    h <- NULL
    cons.options <- c("similarity", "identity", "entropy22", "entropy10")
    if(is.null(conservation))
      conservation=FALSE
    
    if(!is.logical(conservation)) {
      
      if(length(conservation)>1) {
        h <- conservation
        conservation=TRUE
        if(length(h)!=dims.full[2L]) {
          warning("dimension mismatch of supplied 'conservation' vector")
          h <- NULL
          conservation=FALSE
        }
      }
      
      if(!is.logical(conservation) && length(conservation)==1) {
        if(all(conservation %in% cons.options)) {
          conserv.method <- conservation
        }
        else {
          warning("unknown option for 'conservation'")
          conserv.method <- "similarity"
        }
        conservation=TRUE
      }
    }
    else {
      conserv.method <- cons.options[1]
    }
    
    if(is.null(pdbs) & conservation) {
      conservation=FALSE
      warning("forcing 'conservation=FALSE': sequence conservation plot requires the corresponding 'pdbs' object")
    }
    
    if(conservation && is.null(h)) {
      h <- conserv(pdbs, method=conserv.method)
    }
        
    ## x- and ylim
    if(is.null(ylim)) {
      if(spread)
        ylim=c(0, length(unique(col))*offset)
      else
        ylim <- c(0,max(yval, na.rm=TRUE))
    }

    if(is.null(xlim)) {
      xlim <- c(0, ncol(yval))
    }

    ## SSE information
    dots <- list(...)
    sse.aln <- NULL
    if(!is.null(pdbs) & !spread) {
      if( "sse" %in% names(dots) )
        warning("SSE information from 'pdbs' will not be generated when 'sse' is provided")
      else
        sse.aln <- .pdbs2sse(pdbs, ind=1, rm.gaps=rm.gaps)
    }
    
    if( "sse" %in% names(dots) ) {
      sse.aln <- dots$sse
      dots$sse <- NULL
    }

    ## Perform test of significance
    ns <- levels(as.factor(col))
    if((length(ns) !=2) & signif) {
      warning("Number of states is not equal to 2. Ignoring significance test")
      signif <- FALSE
    }

    sig <- NULL
    if(signif) {
      inds1 <- which(col==ns[1])
      inds2 <- which(col==ns[2])

      if(length(inds1)>1 & length(inds2)>1) {
        p <- NULL; q <- NULL
        for(i in 1:ncol(yval)) {
          p <- c(p, t.test(yval[inds1,i],
                           yval[inds2,i],
                           alternative="two.sided")$p.value)
          m <- mean(yval[inds1,i])
          n <- mean(yval[inds2,i])
          q <- c(q, abs(m-n))
        }
        sig <- which(p<pcut & q>qcut)
      }

      ## Plot significance as shaded blocks
      if(is.null(sig))
        warning("Too few data points. Ignoring significance test")

      if(length(sig)==0) {
        ##warning("No significant differences found")
        sig <- NULL
      }
    }

    ## Configure plot
    nrows <- 1
    if(conservation)
      nrows=nrows+1
    if(variance)
      nrows=nrows+1
    
    if(nrows>length(ylab) & !is.null(ylab))
      warning("insufficient y labels")

    op <- par(no.readonly=TRUE)
    on.exit(par(op))

    if(nrows>1)
      par(mfrow=c(nrows,1), mar=mar)
    else
      par(mar=mar)

    plot.new()
    plot.window(xlim=xlim, ylim=ylim, ...)
  
    ## If significance test was performed successfully
    if(!is.null(sig)) {
      ##maxy <- max(yval, na.rm=TRUE)
      bds <- bounds(sig)
      ii <- 1:nrow(bds)
      rect(bds[ii,1], rep(0, length(ii)), bds[ii,2],
           rep(ylim[2], length(ii)),
           col=rep("lightblue", length(ii)), border=NA)
    }

    ## Plot fluctuations / deformations
    par(new=TRUE)
    if(!spread) {
      do.call('plotb3', c(list(x=yval[1,], xlab=xlab, ylab=ylab[1],
                               ylim=ylim, xlim=xlim, type='h', col=1, sse=sse.aln),
                          dots))
      
      ## Plot all lines (col==NA will not be plotted)
      for(i in 1:nrow(yval)) {
        lines( yval[i,], col=col[i], lwd=2, ... )
      }
    }
    else {
      do.call('.plot.enma.spread', c(list(x=yval,
                                          pdbs=pdbs, col=col,
                                          offset=offset,
                                          xlab=xlab, ylab=ylab[1],
                                          ylim=ylim, xlim=xlim),
                                     dots))
    }

    ## Fluctuation / deformations variance
    if (variance) {
      fluct.sd <- apply(yval, 2, var, na.rm=T)
      do.call('plotb3', c(list(x=fluct.sd,
                               xlab=xlab, 
                               ylab=ylab[2],
                               ##ylim=ylim,
                               xlim=xlim,
                               col=1), dots))
    }
    
    ## Plot sequence conservation / entropy
    if (conservation) {
      do.call('plotb3', c(list(x=h,
                               ylab=ylab[3],
                               xlab=xlab,
                               col=1), dots))
    }
    
    out <- list(signif=sig, sse=sse.aln)
    invisible(out)
  }


".plot.enma.spread" <- function(x, pdbs=NULL, col=NULL,
                                 xlab="Alignment Position",
                                 ylab="Fluctuations",
                                 xlim=NULL, ylim=NULL, offset=1, ...) {
  
  if(!inherits(x, "enma") & !inherits(x, "matrix"))
    stop("provide a enma object as obtained from 'nma.pdbs'")

  if(inherits(x, "enma"))
    fluct <- x$fluctuations
  else
    fluct <- x
  
  if(is.null(col))
    stop("group argument missing")

  if(length(col) != nrow(fluct))
    stop("dimension mismatch: col argument should be of same length as x")

  if(length(unique(col)) < 2)
    stop("provide > 2 unique groups")

  row.inds <- which(!is.na(col))
  newfluct <- fluct[row.inds,, drop=FALSE ]
  gaps <- gap.inspect(newfluct)
  
  col.inds <- which(gaps$col < length(row.inds))
  newfluct <- newfluct[, col.inds, drop=FALSE]
  
  ## check for gaps
  gaps <- gap.inspect(newfluct)
  if(any(gaps$col>0))
    rm.gaps <- FALSE
  else
    rm.gaps <- TRUE
  
  sse <- NULL
  dots <- list(...)
  if( "sse" %in% names(dots) ) {
    sse <- dots$sse
    dots$sse <- NULL
  }
  else {
    if(!is.null(pdbs)) {
      if(rm.gaps) {
        gs <- gap.inspect(pdbs$ali)
        pdbs <- trim.pdbs(pdbs, col.inds=gs$f.inds)
      }
      
      pdbs <- trim.pdbs(pdbs, row.inds=row.inds, col.inds=col.inds)
      sse <- .pdbs2sse(pdbs, ind=1, rm.gaps=rm.gaps)
    }
  }
  
  dims <- dim(newfluct)
  if(is.null(xlim))
    xlim <- c(0, dims[2])

  if(is.null(ylim))
    ylim <- c(0, (length(unique(col))*offset)-offset)

  plotb3(newfluct[1, ], col=1, type='l',
         ylab=ylab, xlab=xlab, 
         ylim=ylim, xlim=xlim, sse=sse, ...)
  
  col <- col[!is.na(col)]
  for(i in 1:length(unique(col))) {
    tmpinds <- which(col==i)
    off <- ((i-1)* offset )
    for(j in 1:length(tmpinds))
      lines(newfluct[tmpinds[j], ] + off, col=i)
  }
  
}


".pdbs2sse" <- function(pdbs, ind=1, rm.gaps=FALSE) {
  ind <- ind[1]
  if(file.exists(pdbs$id[ind]))
    id <- pdbs$id[ind]
  else if(file.exists(rownames(pdbs$ali)[ind]))
    id <- rownames(pdbs$ali)[ind]
  
  sse.aln <- NULL
  pdb.ref <- try(read.pdb(id), silent=TRUE)

  if(inherits(pdb.ref, "try-error"))
    pdb.ref <- try(read.pdb(substr(basename(id), 1, 4)), silent=TRUE)

  gaps.res <- gap.inspect(pdbs$ali)
  
  sse.ref <- NULL
  if(!inherits(pdb.ref, "try-error"))
    sse.ref <- try(dssp(pdb.ref), silent=TRUE)

  if(!inherits(sse.ref, "try-error") & !inherits(pdb.ref, "try-error")) {
    if(rm.gaps) {
      resid <- paste0(pdbs$resno[ind, gaps.res$f.inds], pdbs$chain[ind, gaps.res$f.inds])
    }
    else {
      resid <- paste0(pdbs$resno[ind, ], pdbs$chain[ind, ])
    }
    
    ## Helices
    resid.helix <- unbound(sse.ref$helix$start, sse.ref$helix$end)
    resid.helix <- paste0(resid.helix, rep(sse.ref$helix$chain, sse.ref$helix$length))
    inds        <- which(resid %in% resid.helix)
    
    ## inds points now to the position in the alignment where the helices are
    new.sse <- bounds( seq(1, length(resid))[inds] )
    if(length(new.sse) > 0) {
      sse.aln$helix$start  <- new.sse[,"start"]
      sse.aln$helix$end    <- new.sse[,"end"]
      sse.aln$helix$length <- new.sse[,"length"]
    }
    
    ## Sheets
    resid.sheet <- unbound(sse.ref$sheet$start, sse.ref$sheet$end)
    resid.sheet <- paste0(resid.sheet, rep(sse.ref$sheet$chain, sse.ref$sheet$length))
    inds        <- which(resid %in% resid.sheet)
    
    new.sse <- bounds( seq(1, length(resid))[inds] )
    if(length(new.sse) > 0) {
      sse.aln$sheet$start  <- new.sse[,"start"]
      sse.aln$sheet$end    <- new.sse[,"end"]
      sse.aln$sheet$length <- new.sse[,"length"]
    }

    ## SSE vector 
    sse <- rep(" ", length(resid))
    for(i in 1:length(sse.aln$helix$start))
      sse[sse.aln$helix$start[i]:sse.aln$helix$end[i]] <- "H"

    for(i in 1:length(sse.aln$sheet$start))
      sse[sse.aln$sheet$start[i]:sse.aln$sheet$end[i]] <- "E"

    sse.aln$sse <- sse
  }
  else {
    msg <- NULL
    if(inherits(pdb.ref, "try-error"))
      msg = c(msg, paste("File not found:", pdbs$id[1]))
    if(inherits(sse.ref, "try-error"))
      msg = c(msg, "Launching external program 'DSSP' failed")
    
    warning(paste("SSE cannot be drawn", msg, sep="\n  "))
  }

  return(sse.aln)
}
