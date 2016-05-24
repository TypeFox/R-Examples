#' Rownames to symbol-probeID
#' 
#' This function takes an \code{exprs(eset)} matrix where the rownames are
#' probeset IDs and takes an annotated topTable output where you have an ID and
#' Symbol column and outputs a character vector with symbol_probeid for each
#' probeid in rownames(exprs(eset)). You can use this such that the output on a
#' heatmap contains the gene names concatenated to the probe ID in case you have
#' multiple symbols with the same probeID.
#' 
#' @author Stephen Turner
#' @keywords NA
#' 
#' @param exprset The output of \code{exprs(eset)}.
#' @param tt A topTable object.
#'   
#' @return Character vector of the gene symbol with the probe ID.
#' 
#' @examples
#' \dontrun{
#' rownames_to_symprobe(esprs(eset), topTable(fit, number=nrow(fit)))
#' }
#'   
#' @export
rownames_to_symprobe <- function (exprset, tt) {
    if(class(tt)!="data.frame") stop("tt isn't a data frame")
    if(class(exprset)!="matrix") stop("exprset isn't a matrix")
    if(nrow(tt)!=nrow(exprset)) stop("nrow(tt) != nrow(exprset)")
    tt_index <- sapply(rownames(exprset), function(pattern) grep(pattern, tt$ID))
    newrownames <- paste(tt$Symbol[tt_index], tt$ID[tt_index], sep="_")
    newrownames
}
