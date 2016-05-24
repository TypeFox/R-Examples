# For finding detections within a detectionList object
# Used when score cutoffs have been updated
# Modified: 2013 JUNE 10

findDetections <-
function(
detection.obj, 
keep.verify=FALSE
) {

   for(i in names(detection.obj@templates)) {

      score.cutoff <- detection.obj@templates[[i]]@score.cutoff
      pks <- detection.obj@peaks[[i]]
      pks$detection <- FALSE
      pks$detection[pks$score>=score.cutoff] <- TRUE
      hits <- pks[pks$score>=score.cutoff, ]
      hits$detection <- NULL
      if(nrow(hits)>0) rownames(hits) <- 1:nrow(hits)
      #if(keep.verify && ncol(detection.obj@detections[[i]]) == 5 && names(detection.obj@detections[[i]])[5] == "true") hits <- merge(hits, detection.obj@detections[[i]], all=TRUE)
      if(keep.verify && ncol(detection.obj@detections[[i]]) == 5) hits <- merge(hits, detection.obj@detections[[i]], all=TRUE)
      detection.obj@detections[[i]] <- hits
      detection.obj@peaks[[i]] <- pks

   }

   return(detection.obj)

}
