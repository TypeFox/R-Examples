importBlastTab <- function(file){
  out <- read.table(file,stringsAsFactors=FALSE, sep="\t")
  out
}