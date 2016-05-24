#' intersectMatrix
#'
#' This function will count the intersection tally of AVS and genomic regions
#' @param avs A GRanges object which is outputted by loadLd function
#' @param regions A data frame. The data frame contains sample sheet identical to DiffBind or ChIPQC input sample sheets
#' @param ... Heatmap arguments
#' @keywords AVS,Granges
#' @examples
#' \dontrun{
#' intersectMatrix(avs,
#'                 regions=samples,
#'                 col=c("white","grey10"),
#'                 scale="none",
#'                 margins=c(10,5),
#'                 cexRow = 1,
#'                 cexCol = 0.5,
#'                 Rowv=NA,
#'                 Colv=NA)
#' }
#' @import GenomicRanges
#' @export
intersectMatrix <- function(avs, regions, ...){
  no_of_tags <- length(avs)
  no_of_beds <- length(regions$SampleID)
  overlap <- matrix(NA, nrow = no_of_beds, ncol = no_of_tags)
  for (i in 1:no_of_beds){
    bed_path <- as.character(regions$Peaks[i])
    bed.gr <- bedToGRanges(bed_path)
    overlap[i,] <- ifelse(countOverlaps(avs, bed.gr, ignore.strand=TRUE)>0,1,0)
  }
  row.names(overlap) <- regions$SampleID
  avs.ids <- list()
  for (i in 1:length(avs)){
    avs.ids <- c(avs.ids, as.character(elementMetadata(avs[[i]])[1,2]))
  }
  colnames(overlap) <- avs.ids
  heatmap(overlap, ...)
  return(overlap)
}
