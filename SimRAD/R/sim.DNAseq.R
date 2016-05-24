sim.DNAseq <-
function(size=10000, GCfreq=0.46){
  ATfreq <- 1-GCfreq
  simseq <- paste(sample(c("A", "C", "G", "T"), size=size, replace=TRUE, prob=c(ATfreq/2, GCfreq/2, GCfreq/2, ATfreq/2)), collapse="")
}
