getProbsRegions <- function(probs, regions, intensities, nclass=3)
 {
  blocks<-regions@featureData@data
  annotation<-intensities[,1:4]
  probsBlock<-lapply(1:nrow(blocks), getProbsRegions.i, blocks=blocks, probs=probs, annotation=annotation, nclass=nclass)
  probsBlock
 }
