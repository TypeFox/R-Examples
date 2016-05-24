#' bedToGRanges
#'
#' This function will convert a bed to GRanges object
#' @param file A bed file. Must contain at least three columns: chr, start and end.
#' @keywords bed,Granges
#' @examples
#' \dontrun{
#' bedToGRanges(file)
#' }
#' @import GenomicRanges
#' @export
bedToGRanges <- function(file){
  df <- read.table(file, sep="\t")
  if (is.numeric(df[1,2])){
    cnames <- c("chr","start","end")
    colnames(df)<-cnames;
  } else {
    cnames <- colnames(df)
  }
  if (!("chr" %in% cnames)){
    stop("No chr column found in df");
  }
  if (!("start" %in% cnames)){
    stop("No start column found in df");
  }
  if (!("end" %in% cnames | "stop" %in% cnames)){
    stop("No end column found in df");
  }
  df <- df[,-c(4:ncol(df))]
  df$chr<- as.factor(unlist(lapply(df$chr, function(x) gsub("chr", "", x))))
  makeGRangesFromDataFrame(df, keep.extra.columns=FALSE)
}
