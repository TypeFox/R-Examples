#' Tidy DESeq2 result
#' 
#' Returns a tidy version of a DESeq2 results table.
#' 
#' @author Stephen Turner
#' @keywords NA
#' 
#' @param deseqresult Results from running \code{results(dds)} on a DESeqDataSet object.
#'   
#' @return a tidy version of the DESeq2 results.
#' 
#' @import dplyr
#' 
#' @examples
#' \dontrun{
#' res <- results(dds)
#' res <- deseqresult2tbl
#' }
#'   
#' @export
deseqresult2tbl <- function(deseqresult) {
    gene=baseMean=padj=pvalue=NULL #gets rid of the note on pkg check
    if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
    deseqresult <- as.data.frame(deseqresult)
    deseqresult$gene <- rownames(deseqresult)
    rownames(deseqresult) <- NULL
    deseqresult <- tbl_df(deseqresult)
    deseqresult <- select(deseqresult, gene, baseMean:padj)
    deseqresult %>% arrange(padj,pvalue, baseMean)
}