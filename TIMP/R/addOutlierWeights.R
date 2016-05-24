"addOutlierWeights" <- function(m, weightListOutliers) {
  for(i in 1:length(m))
    m[[i]]@weightM <- weightListOutliers[[i]]
 
  m


}
