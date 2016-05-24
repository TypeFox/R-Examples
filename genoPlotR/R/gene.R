################################
# gene class and methods
################################
# deprecated
## gene <- function(name, start, end, strand, col=NULL, lwd=1, lty=1, pch=8, cex=1,
##                  gene_type="arrows"){
##   if (missing(name)) stop("Name must be provided")
##   if (missing(start)) stop ("Start must be provided")
##   if (missing(end)) stop ("End must be provided")
##   if (missing(strand)) stop ("Strand must be provided")
##   gene <- data.frame(name=name, start=start, end=end, strand=strand, col=col,
##                      lwd=lwd, lty=lty, pch=pch, cex=cex, gene_type=gene_type,
##                      stringsAsFactors=FALSE)
##   as.gene(gene)
## }
## as.gene <- function(gene, col="blue", lwd=1, lty=1, pch=8, cex=1,
##                     gene_type="arrows"){
##   # return self is already of the right class
##   if (is.gene(gene)) {
##     return(gene)
##   }
##   # check for presence of arguments
##   if (is.null(gene$name)) stop("No gene name")
##   if (is.null(gene$start)) stop("No start in gene")
##   if (is.null(gene$end)) stop("No end in gene")
##   if (is.null(gene$strand)) stop("No strand in gene")
##   if (is.null(gene$col)) gene$col <- col
##   if (is.null(gene$lwd)) gene$lwd <- lwd
##   if (is.null(gene$lty)) gene$lty <- lty
##   if (is.null(gene$pch)) gene$pch <- pch
##   if (is.null(gene$cex)) gene$cex <- cex
##   # gene_type: not given in gene
##   if (is.null(gene$gene_type)){
##     if (!is.character(gene_type) || !(gene_type %in% gene_types()))
##       stop(paste("gene_type muste be one of:",
##                  paste(gene_types(), collapse=", ")))
##     if (gene_type == "auto") gene_type <- "arrows"
##     gene$gene_type <- gene_type
##   }
##   # given in gene
##   else {
##     if (!is.character(gene$gene_type) || !all(gene_type %in% gene_types()))
##       stop(paste("gene_type muste be one of:",
##                  paste(gene_types(), collapse=", ")))   
##   }
##   # check for correct argument types
##   if (!is.character(gene$name)) stop("Non-character name")
##   if (!(is.numeric(gene$start) && is.numeric(gene$end)))
##     stop("Non-numeric start or end")
##   if (!all(is.numeric(c(gene$lwd, gene$lty, gene$pch, gene$cex))))
##     stop("lwd, lty, pch and cex must be numeric")
##   # care about strand
##   if (gene$strand == 1 | gene$strand == "+") {
##     gene$strand <- 1
##   } else if (gene$strand == -1 | gene$strand == "-") {
##     gene$strand <- -1
##   } else {
##     stop("Invalid gene strand: should be 1, -1, + or -")
##   }
##   # care about start/end
##   if (gene$end < gene$start) stop("End smaller than start. Use strand instead")
##   # no color check for the moment
##   # convert to df
##   gene <- as.data.frame(gene, stringsAsFactors=FALSE)
##   # class attribution, no check for other at the moment
##   class(gene) <- c("gene", "data.frame")
##   gene
## }
## is.gene <- function(gene){ 
##   inherits(gene, "gene")
## }
## # calculates x range
## rang.gene <- function(gene){
##   range(c(gene$start, gene$end))
## }
## rang.list <- rang.gene
