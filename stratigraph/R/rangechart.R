"rangechart" <-
function(x = NULL,
         counts = NULL,
         depths = NULL,
         sample.labels = NULL,
         taxa = NULL,
         short.names = NULL,
         higher.grp = NULL,
         tax.cat = NULL,
         reorder = NULL,
         plot.points = FALSE,
         plot.depths.increasing.down = TRUE,
         llwd = 2,
         cex.xaxis = 0.5,
         cex.yaxis = 1,
         ...){

if(is.strat.column(x)){
  counts <- x$counts
  if(is.null(depths)) depths <- x$depths
  if(is.null(taxa)) taxa <- x$taxa
  if(is.null(short.names)) short.names <- x$short.names
  if(is.null(higher.grp)) higher.grp <- x$higher.grp
  if(is.null(tax.cat)) tax.cat <- x$tax.cat
}

# Read in the count data:
if(is.data.frame(counts) || is.matrix(counts)){
  counts <- apply(counts, 2, 'as.numeric')
}else if(is.character(counts) && length(counts) == 1){
  counts <- read.csv(counts, header = TRUE, skip = 0, colClasses = '')
  counts <- apply(counts, 2, 'as.numeric')
}else{
  stop('argument to counts not understood')
}

# Check that depths is the right length
if(is.null(depths)){
  warning('no depths provided; plotting samples at regular intervals')
  depths <- 1:nrow(counts)
}
depths <- as.numeric(depths)
if(length(depths) != nrow(counts)){
  stop(paste(length(depths), ' depths, and ', nrow(counts), ' rows in the count matrix.', sep = ''))
}

# Replace NA in counts with zeros
if(sum(is.na(counts)) > 0){
  warning(paste(sum(is.na(counts)), 'missing values in count matrix replaced with zeros'))
  counts[is.na(counts)] <- 0
}

# Remove empty columns
emptycols <- !(colSums(counts, na.rm = TRUE) > 0)
if(any(emptycols)){
  counts <- counts[,!emptycols]
  tax.cat <- tax.cat[!emptycols]
  taxa <- taxa[!emptycols]
  warning(paste(sum(emptycols),
                'columns with zero counts at all levels removed'))
}

# Deal with reorder
if(!is.null(reorder)){
  if(length(reorder) == 1 && is.character(reorder)){
    funny <- function(x) return((1:length(x))[x > 0])
    if(pmatch(reorder, 'fad.by.category', nomatch = FALSE)){
      fads <- depths[as.numeric(lapply(apply(counts, 2, funny), max))]
      reorder.vect <- sort(fads, decreasing = TRUE,
                             index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
      reorder.vect <- sort(as.character(tax.cat),
                             index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
    }else if(pmatch(reorder, 'lad.by.category', nomatch = FALSE)){
      lads <- depths[as.numeric(lapply(apply(counts, 2, funny), min))]
      reorder.vect <- sort(lads, decreasing = TRUE,
                      index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
      reorder.vect <- sort(as.character(tax.cat),
                      index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
    }else if(pmatch(reorder, 'by.count', nomatch = FALSE)){
      reorder.vect <- sort(colSums(counts), decreasing = TRUE,
                           index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
    }
  }else if(length(reorder) == ncol(counts)){
    reorder.vect <- as.numeric(reorder)
    counts <- counts[,reorder.vect]
  }else{
    stop('argument to reorder not understood')
  }
}

# Deal with tax.cat
if(!is.null(tax.cat)){
  if(length(tax.cat) != ncol(counts)){
    warning('taxon category labels seem to be the wrong length')
    tax.cat <- NULL
  }
}

if(is.null(tax.cat)){
  tax.cat <- rep('', ncol(counts))
}

ad <- a.datums(strat.column(counts = counts, depths = depths))

plot(1:ncol(counts), ylim = c(max(depths), min(depths)),
     type = 'n', xaxt = 'n', yaxt = 'n',
     xlab = '', ylab = 'Depth')
segments(1:ncol(counts), ad[,1],
         1:ncol(counts), ad[,2],
         lwd = llwd, col = tax.cat)
segments(1:ncol(counts), ad[,1],
         1:ncol(counts), rep(par()$usr[3], ncol(counts)),
         col = grey(0.5), lty = 3, lwd = 0.5)
axis(1, at = 1:ncol(counts), labels = colnames(counts),
     cex.axis = cex.xaxis, las = 3)
axis(2, las = 1, cex.axis = cex.yaxis)
axis(2, at = depths, labels = FALSE, tck = -0.01)
for(i in 1:ncol(counts)){
  plocs <- depths[(counts > 0)[,i]]
  points(rep(i, length(plocs)), plocs, col = tax.cat[i], ...)
}
} # End of function