group.heatmap.simple <- function(data, is.OTU=TRUE, meta=NULL, rank="g", row.factor=NULL,
                                 top=NULL, count=FALSE, drop.unclassified=FALSE,
                                 dendro="none", file=NULL, ext=NULL, 
                                 width=9, height=8, leg.x=-0.08, leg.y=0) {

  save <- !is.null(file)
  if ( is.OTU && !is.null(rank) ) {
    valid.OTU(data)
    .valid.meta(otu1=data, meta=meta)
    .valid.rank(rank)
    data.tax <- tax.abund(data, rank=rank, drop.unclassified=drop.unclassified,
                         count=count, top=top)  
  } else if ( !is.OTU || is.null(rank) ) {
    rank <- NULL    
    if ( count ) {
      stand.method=NULL
    } else {
      stand.method="total"
    }    
    data.tax <- data.revamp(data=list(df=data), is.OTU=is.OTU, 
                          stand.method=stand.method, 
                          ranks=rank, top=top)[[1]]
  }
  
  data.tax <- data.tax[match(rownames(meta), rownames(data.tax)),]
  if ( !identical(rownames(meta), rownames(data.tax)) ) {
    stop("samples are not the same in metadata and community data")
  }
 
  # we create a NULL col.factor for compatability with the internal
  # .valid.factor and .factor.settings functions
  # Note: a column factor is not valid here, as the factors are for samples
  # which are on the rows (the aggregated taxonomy is on the columns)
  col.factor <- NULL
  given.rfactor <- !is.null(row.factor)
  given.cfactor <- !is.null(col.factor)
  
  if (is.null(row.factor)) {
    rfactor <- NULL
  }
  
  if (given.rfactor) {
    rfactor <- .valid.factors(meta, row.factor, max.factors = 1)
    
    # if we have a row factor, but the length doesn't match the number of samples,
    # stop
    if (nrow(rfactor) != nrow(data.tax)) {
      stop(paste("The given data has", nrow(data.tax) - 1, 
                 "samples, but the row.factor has", nrow(rfactor), 
                 "points. These two values must match to continue."))
    }
  }
  
  if (!is.null(meta) && !all(rownames(data.tax) %in% rownames(meta))) {
    stop("the rownames for 'data' and 'meta' do not match.")
  }
  
  # general approach for this function: we build up a list of the arguments
  # we wish to pass to heatmap.2 and legend based on arguments supplied by 
  # the user and some default settings. Once we have created a list with all
  # the arguments, we pass that to do.call (with heatmap.2 or legend) to 
  # create the plot.
  
  if (save) {
    # does data validation for file
    .get.dev(file, ext, height=height, width=width)
  }
  
  # generate a palette with a large gradient
  cols <- RColorBrewer::brewer.pal(3, "YlOrRd") 
  gradient <- colorRampPalette(cols)(100)
  
  longest.row.label <- max(sapply(rownames(data.tax), FUN=nchar))
  longest.col.label <- max(sapply(colnames(data.tax), FUN=nchar))
  longest.factor.label <- max(sapply(levels(factor(meta[[row.factor]])), FUN=nchar))
  
  hmap.args <- list(x=as.matrix(data.tax), dendrogram=dendro, trace = "none",
                    key = TRUE, keysize = 1, density.info = c("none"), 
                    col=gradient, xpd=NA, 
                    # we set the margins based on the length of the labels
                    ### (this probably needs more tweaking!)
                    margins=c(0.6 * longest.col.label, 0.7 * longest.row.label ))
  
  hmap.args$dendrogram <- dendro
  # order based on metadata, if appropriate
  if (given.rfactor) {
    # ensure that the order of the samples matches the order of the metadata
    # exactly, after grouping all identical samples together
    data.tax <- data.tax[order(rfactor[ ,1]), ]
    hmap.args$x <- as.matrix(data.tax)
    
    if (dendro == "row" || dendro == "both") {
      # warning("if metadata is provided, clustering will not occur and no dendrogram will be shown for the rows.")
      leg.args <- list(x="right", inset=c(leg.x, leg.y), xpd=NA, cex=0.7)
      if ( dendro == "row" ) {
        hmap.args <- c(hmap.args, Rowv=TRUE, Colv=FALSE)
      } else {
        hmap.args <- c(hmap.args, Rowv=TRUE, Colv=TRUE)
      }
      #hmap.args$dendrogram <- "none"
      #leg.args <-list(x=locator(), xpd=NA)
      hmap.args$margins <- c(0.6 * longest.col.label, 0.7 * longest.row.label +longest.factor.label)
    } else {
      hmap.args <- c(hmap.args, Colv=TRUE, Rowv=FALSE)
      #leg.args <-list(x=locator(), xpd=NA)
      leg.args <- list(x="left", inset=c(leg.x, leg.y), xpd=NA, cex=0.7)
      hmap.args$margins <- c(0.6 * longest.col.label, 0.7 * longest.row.label )
    }
  }
  
  #leg.args <-list(x=locator(), xpd=NA)
  # leg.args <- list(x="right", inset=c(leg.x, leg.y), xpd=NA, cex=0.7)
  
  args <- .factor.settings(rfactor, NULL, hmap.args, leg.args)
  hmap.args <- args[[1]]
  leg.args <- args[[2]]
  
  do.call(heatmap.2, hmap.args)
  
  if (given.rfactor) {
    do.call(legend, leg.args)
  }
  
  if (save) {
    dev.off()
  }
  
  invisible()
}

dissim.heatmap <- function(data, is.OTU=TRUE, meta=NULL, row.factor=NULL, col.factor=NULL, 
                           stand.method="chi.square", dissim.method="euclidean",
                           file=NULL, ext=NULL, height=8, width=9, leg.x=-0.05, leg.y=0) {
  
  save <- !is.null(file)
  
  given.rfactor <- !is.null(row.factor)
  given.cfactor <- !is.null(col.factor)
  
  # general approach for this function: we build up a list of the arguments
  # we wish to pass to heatmap.2 and legend based on arguments supplied by 
  # the user and some default settings. Once we have created a list with all
  # the arguments, we pass that to do.call (with heatmap.2 or legend) to 
  # create the plot.
  
  # validate row and column factors, if supplied by user
  # (this is very ugly and should be cleaned up)
  if (given.rfactor) { 
    rfactor <- .valid.factors(meta, row.factor, max.factors = 1)
  } else {
    rfactor <- NULL
  }
  if (given.cfactor) {
    cfactor <- .valid.factors(meta, col.factor, max.factors = 1)
  } else {
    cfactor <- NULL
  }
  
  # abuse logical arithmetic
  num.factors <- given.rfactor + given.cfactor
  
  # calculate dissimilarity distances
  if ( is.OTU ) {
    otu.t <- transpose.OTU(data)
  } else {
    otu.t <- data
  }

  otu.t <- otu.t[match(rownames(meta), rownames(otu.t)),]
  if ( !identical(rownames(meta), rownames(otu.t)) ) {
    stop("samples are not the same in metadata and community data")
  }
 
  dissim <- decostand(otu.t, method=stand.method)
  distances <- vegdist(dissim, method=dissim.method)
  
  distances <- as.matrix(distances)
  
  # order the data based on how we will order the metadata
  
  if (given.rfactor) {
    sample.order <- order(rfactor[ ,1])
    
  } else if (given.cfactor) {
    sample.order <- order(cfactor[ ,1])
  }
  
  if (given.rfactor || given.cfactor) {
    distances <- distances[sample.order, sample.order]
  }
  
  # get colour palette
  mat.cols <- RColorBrewer::brewer.pal(11, "PiYG")
  
  # "basic" plot parameters; they mostly determine the non-meta data as well 
  # as the properties of the plotting environment
  hmap.args <- list(x=distances, Rowv=FALSE, Colv=FALSE, dendrogram="none",
                    trace="none", key=TRUE, keysize=1, density.info="none",
                    col=mat.cols, margins=c(width, height))
  
  if (save) {
    # does data validation for file
    .get.dev(file, ext, width=width, height=height)
  }
  
  leg.args <- list(x="left", inset=c(leg.x, leg.y), cex=0.7, xpd=NA)

  args <- .factor.settings(rfactor, cfactor, hmap.args, leg.args)
  
  hmap.args <- args[[1]]
  leg.args <- args[[2]]
  
  hmap.args$main <- paste0(.capitalize(stand.method), "/", .capitalize(dissim.method),
                           " Distance Matrix")
  
  do.call(heatmap.2, hmap.args)
  
  if (given.rfactor || given.cfactor) {
    do.call(legend, leg.args)
  }
  
  if (save) {
    dev.off()
  }
  
  invisible()
}


group.heatmap <- function(data, is.OTU=TRUE, meta, rank, factors, 
                          top=25, remove.unclassified=TRUE, 
                          stand.method=NULL,
                          dist.method="bray",
                          hclust.method="average",
                          dendro.row.status="yes",
                          dendro.col.status="hidden",
                          row.labels=TRUE, row.cex=1,
                          cut=NULL, file=NULL, ext=NULL, 
                          width=9, height=9) {
  
  # I have seen this discussion: http://yihui.name/en/2014/07/library-vs-require/
  # but I think returning an explanatory error message is worthwhile
  
 #if ( require("Heatplus") ) {
 #  Heatplus::annHeatmap2
 #  Heatplus::niceBreaks
  if(!requireNamespace('Heatplus')) {
#   } else {
    stop("package 'Heatplus' is required to use this function: follow the instructions at http://www.bioconductor.org/packages/release/bioc/html/Heatplus.html to install it.")
  }

  if(!requireNamespace('RColorBrewer')) {
    stop("package 'RColorBrewer' is required to use this function")
  }

  save <- !is.null(file)
  if ( is.OTU ) {
    valid.OTU(data)
    .valid.meta(otu1=data, meta=meta)
    if ( !is.null(rank) ) {
      .valid.rank(rank)
       rank <- rank
       # obtain tax.abund matrix in counts, remain all groups, including unclassified ones.
       tax <- tax.abund(data, rank=rank, drop.unclassified=FALSE, count=TRUE )
    } else {
      rank <- NULL
      tax <- data.revamp(data=list(df=data), is.OTU=is.OTU, 
                          stand.method=NULL,ranks=rank, top=NULL)[[1]]
    }
  } else {
    rank <- NULL    
    tax <- data
  }

  tax <- tax[match(rownames(meta), rownames(tax)),]
  if ( !identical(rownames(meta), rownames(tax)) ) {
    stop("samples are not the same in metadata and community data")
  }
 
  # order the taxa
  tax <- tax[, order(colSums(tax), decreasing=TRUE)]

  if(is.null(stand.method)) {
    stand <-tax
  } else {
    stand <- decostand(tax, stand.method)
  }
  
  # drop.unclassified using blacklist
  if(remove.unclassified && !is.null(rank) ) {
    remove.pat <- gsub(.get.rank.pat(rank), "", paste0(.blacklist(.get.rank.pat(rank)), "|no_taxonomy"))
    stand <- stand[ , !grepl(remove.pat, names(stand), ignore.case=TRUE), drop=FALSE]
  } else {
    stand<-stand
  }

  if ( ncol(stand) <= top ) {
    top <- ncol(stand)
  } else {
    top <- top
  }

  stand <- stand[, 1:top]
  
  row.clus <- hclust(vegdist(stand, method=dist.method), hclust.method)
  col.clus <- hclust(vegdist(t(stand), method=dist.method), hclust.method)
  
  # generate a palette with a large gradient
  base.cols <- RColorBrewer::brewer.pal(3, "YlOrRd")
  # this next part is a hack: the niceBreaks call is taken from the annHeatmap2
  # source; we use that call to calculate the minimum number of colours needed
  # which means that the majority of the palette is used (looks nice, distinguishes best)
  # cols.needed <- length(Heatplus::niceBreaks(range(as.matrix(stand), na.rm = TRUE), 256)) - 1
  cols.needed <- length(Heatplus::niceBreaks(range(as.matrix(stand), na.rm = TRUE), 256)) - 1
  cols <- colorRampPalette(base.cols)(cols.needed)
  
  #pastels <- function(n) {brewer.pal(n, "Pastel1")}
  col.pal <- c(RColorBrewer::brewer.pal(9, "Pastel1"), RColorBrewer::brewer.pal(8, "Accent"), .ram.pal(20))
  pastels <- function(n){col.pal[1:n]}
  
  meta.factors <- .valid.factors(meta, factors)
  
  # find the longest group label, we will allocate margin space based on this
  longest.label <- max(sapply(colnames(stand), FUN=nchar))

  if(!is.null(stand.method)) {
    scale=c("none")
  } else {
    scale=c("row")
  }

  ann.args <- list(x=as.matrix(stand), col=cols, scale=scale,
                   ann = list(Row = list(data = meta.factors)),
                   legend = 3,
       labels = list(Col=list(nrow=longest.label * 0.4), 
                                 Row=list(nrow=6)),
       dendrogram=list(Row=list(dendro=as.dendrogram(row.clus), status=dendro.row.status),
           Col=list(dendro=as.dendrogram(col.clus), status=dendro.col.status))
       )
  
  if ( row.labels ) {
    ann.args$labels = list(Col=list(nrow=longest.label * 0.4), 
                                 Row=list(nrow=6, cex=row.cex))
  } else {
    ann.args$labels = list(Col=list(nrow=longest.label * 0.4), 
                                 Row=list(labels=NULL,nrow=6, cex=row.cex))
  }
  
  # if the user wants to cut the dendrogram, add the correct arguments
  if (!is.null(cut)) {
    ann.args$cluster <- list(Row = list(cuth = cut, col = pastels))
  }
  
  ann.plot <- do.call(Heatplus::annHeatmap2, ann.args)
  #ann.plot <- do.call(annHeatmap2, ann.args)
  
  if (save) {
    .get.dev(file, ext, height=height, width=width)
  }
  
  plot(ann.plot)
  # ann.plot
  
  if (save) {
    dev.off()
  }
  
  invisible(ann.plot)
}

