gethits <-
function(hitsfile){
  data <- read.table(file = hitsfile)
  chromosomes <- unique(data$V1)
  if(length(chromosomes)>1){stop("You have provided hits originating from more than one chromosome. Please parse chromosomes separately")}
  hits <- data$V2
  return(hits)
}
