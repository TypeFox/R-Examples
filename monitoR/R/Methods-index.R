# S. Hafner
# Modified: 2015 JULY 21

setMethod('[', 'TemplateList', 
   function(x, i=names(x@templates)) {
      x@templates <- x@templates[i]
      return(x)
   }
)

setMethod('[', 'templateScores', 
   function(x, i=names(x@templates)) {
      x@templates <- x@templates[i]
      x@scores <- x@scores[i]
      return(x)
   }
)

setMethod('[', 'detectionList', 
   function(x, i=names(x@templates)) {
      x@templates <- x@templates[i]
      x@scores <- x@scores[i]
      x@peaks <- x@peaks[i]
      x@detections <- x@detections[i]
      return(x)
   }
)


