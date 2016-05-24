################################
# dna_seg class and methods
################################
dna_seg <- function(x, ...){
  # support for list ?
  if (is.data.frame(x)){
    return(as.dna_seg(x, ...))
  } else if (is.list(x)) {
    dna_seg(as.data.frame(x, stringsAsFactors=FALSE))
  } else {
    stop(paste("Cannot coerce class", class(x), "to dna_seg"))
  }
}
# convert to dna_seg format. 
as.dna_seg <- function(df, col="blue", lty=1, lwd=1, pch=8, cex=1,
                       gene_type="arrows"){
  # check for class dna_seg, list, df
  if (is.dna_seg(dna_seg)) return(df)
  #if (is.list(df) && !is.data.frame(df)) df <- dna_seg(df)
  if (is.data.frame(df)) {
    # check that it has rows
    if (nrow(df) < 1){
      stop("Number of rows is 0, check data input")
    }
    # check for columns
    names <- c("name", "start", "end", "strand")
    if (!identical(names(df)[1:4], names))
      stop("Col names should start with name, start, end, strand")
    if (!all(sapply(df, function(x) all(!is.null(x)))))
      stop("NULL values not allowed in data.frame")
    if (!is.numeric(df$start) | !is.numeric(df$end)){
      stop("Start and end must be numeric")
    }
    if (is.factor(df$name)) df$name <- as.character(df$name)
    if (is.factor(df$strand)) df$strand <- as.character(df$strand)
    if (is.factor(df$col)) df$col <- as.character(df$col)
    if (is.factor(df$gene_type)) df$gene_type <- as.character(df$gene_type)
    # care for strand
    if (is.character(df$strand)) {
      df$strand[df$strand=="+"] <- 1
      df$strand[df$strand=="-"] <- -1
      df$strand <- as.numeric(df$strand)
    }
    if (!all(df$strand %in% c(-1, 1)))
      stop("Strand vector must be composed of 1, -1, - and +, uniquely")
    # col
    if (is.null(df$col)) df$col <- col
    # lwd & lty
    if (is.null(df$lty)) df$lty <- lty
    if (is.null(df$lwd)) df$lwd <- lwd
    if (is.null(df$pch)) df$pch <- pch
    if (is.null(df$cex)) df$cex <- cex
    # gene_type: not given in gene
    if (is.null(df$gene_type)){
      if (gene_type == "auto") gene_type <- auto_gene_type(nrow(df))
      df$gene_type <- gene_type
    }
    if (is.null(df$gene_type)) df$gene_type <- gene_type
    # check for correct argument types
    if (!is.character(df$name)) stop("Non-character name")
    if (!(is.numeric(df$start) && is.numeric(df$end)))
      stop("Non-numeric start or end")
    if (!all(is.numeric(c(df$lwd, df$lty, df$pch, df$cex))))
      stop("lwd, lty, pch and cex must be numeric")
    if (!is.character(df$gene_type))
      stop(paste("gene_type must be a character vector, made of:",
                 paste(gene_types(), collapse=", "), "or a function name"))
  }
  else {
    stop("Unable to handle this format")
  }
  class(df) <- c("dna_seg", "data.frame")
  return(df)
}
is.dna_seg <- function(dna_seg){
  inherits(dna_seg, "dna_seg")
}
# calculates dna_seg range
range.dna_seg <- function(x, ...){
  dna_seg <- x
  range(dna_seg$start, dna_seg$end, na.rm=FALSE)
}
#range.data.frame <- range.dna_seg
# trim dna_seg given x limits
trim.dna_seg <- function(x, xlim=NULL, ...){
  dna_seg <- x
  if (!is.null(xlim)){
    if (!is.numeric(xlim)) stop("xlim must be numeric")
    if (length(xlim) != 2) stop("xlim must be length 2")
    dna_seg <- dna_seg[dna_seg$start >= xlim[1] & dna_seg$end <= xlim[2],]
  }
  dna_seg
}
# reverse dna_seg
reverse.dna_seg <- function(x, ...){
  dna_seg <- x
  start <- -dna_seg$end
  dna_seg$end <- -dna_seg$start
  dna_seg$start <- start
  dna_seg$strand <- -dna_seg$strand
  dna_seg
}
# merges two dna_segs. So far returns the minimal common set of columns, maybe
# changed in future
c.dna_seg <- function(...){
  # a little helper function to grab colnames
  #fold <- function(f, x, L) (for(e in L) x <- f(x, e))
  fold <- function(x, fun) {
    if (length(x) == 1) return(fun(x))
    accumulator <- fun(x[[1]], x[[2]])
    if (length(x) == 2) return(accumulator)
    for(i in 3:length(x)) {
      accumulator <- fun(accumulator, x[[i]])
    }
    accumulator
  }
  # parse args
  x <- list(...)
  n <- length(x)
  if (!all(sapply(x, is.dna_seg))) 
    stop("All elements must be of class dna_seg")
  if (n == 1) return(x)
  cols <- lapply(x, names)
  com_cols <- fold(cols, intersect)
  dna_seg <- x[[1]][com_cols]
  for (i in 2:n){
    dna_seg <- rbind(dna_seg, x[[i]][com_cols])
  }
  dna_seg
}
