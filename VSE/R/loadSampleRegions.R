#' loadSampleRegions
#'
#' This function will download sample bed files from www.hansenhelab.org/VSE/sample_regions in /VSE_samples
#' @keywords sample,histone,MCF7
#' @examples
#' \dontrun{
#' loadSampleRegions()
#' }
#' @description The sample bed files are DNAse-seq, ChIP-seq for H3K4me1, H3K4me3, H3K27ac, H3K27me3 and H3K36me3 for MCF7 cells. The data are obtained from ENCODE project. There is also one sampleSheet.csv which is the sample sheet for the bed regions in the format similar to ChIPQC or DiffBind requirement.
#' @return A directory names VSE_samples that will contain 6 bed files and 1 sampleSheet.csv
#' @importFrom utils write.csv download.file
#' @export
loadSampleRegions <- function(){
  tmpdir <- tempdir()
  tmpdir <- paste(tmpdir,"VSE_samples",sep="/")
  dir.create(tmpdir)
  beds <- c("DHS.bed","H3K27ac.bed","H3K27me3.bed","H3K36me3.bed","H3K4me1.bed","H3K4me3.bed","sampleSheet.csv")
  for (i in 1:length(beds)){
    url <- paste0("http://www.hansenhelab.org/VSE/sample_regions/",beds[i])
    download.file(url, destfile = paste(tmpdir,beds[i],sep="/"), method = "curl")
  }
  samp_path <- paste(tmpdir,"sampleSheet.csv",sep="/")
  samp <- read.csv(samp_path, header=TRUE)
  samp$Peaks <- paste(tmpdir,samp$Peaks,sep="/")
  write.csv(samp, file = samp_path, quote = FALSE)
  return(samp_path)
}
