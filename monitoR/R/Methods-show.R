# Show methods for objects
# S. Hafner
# Modified: 2015 JULY 21

setMethod('show', 'corTemplateList', 
   function(object) {

      templates <- object@templates
      dat <- data.frame(
         'original.recording'=sapply(templates, function(x) x@clip.path), 
         'sample.rate'=       sapply(templates, function(x) x@samp.rate), 
         'lower.frequency'=   sapply(templates, function(x) round(x@frq.lim[1], 3)), 
         'upper.frequency'=   sapply(templates, function(x) round(x@frq.lim[2], 3)), 
         'lower.amp'=         sapply(templates, function(x) round(min(x@pts[, 'amp']), 3)), 
         'upper.amp'=         sapply(templates, function(x) round(max(x@pts[, 'amp']), 3)), 
         'duration'=          sapply(templates, function(x) round(x@duration, 3)), 
         'n.points'=          sapply(templates, function(x) nrow(x@pts)), 
         'score.cutoff'=      sapply(templates, function(x) x@score.cutoff)
      )

      cat('\nObject of class \"', class(object), '\"\n', sep='')
      cat('\tcontaining ', length(object@templates), ' templates\n')
      print(dat)
   }
)

setMethod('show', 'binTemplateList', 
   function(object) {
      templates <- object@templates
      dat <- data.frame(
         'original.recording'=sapply(templates, function(x) x@clip.path), 
         'sample.rate'=sapply(templates, function(x) x@samp.rate), 
         'lower.frequency'=   sapply(templates, function(x) x@frq.lim[1]), 
         'upper.frequency'=   sapply(templates, function(x) x@frq.lim[2]), 
         'duration'=          sapply(templates, function(x) round(x@duration, 2)), 
         'on.points'=         sapply(templates, function(x) nrow(x@pt.on)), 
         'off.points'=        sapply(templates, function(x) nrow(x@pt.off)), 
         'score.cutoff'=      sapply(templates, function(x) x@score.cutoff)
      )

      cat('\nObject of class \"', class(object), '\"\n', sep='')
      cat('\n\tcontaining ', length(object@templates), ' templates\n')
      print(dat)
   }
)

setMethod('show', 'templateScores', 
   function(object) {

      cat('\nA \"templateScores\" object\n')

      cat('\nBased on the survey file: ', object@survey.name, '\n')

      cat('\nAnd ', length(object@templates), ' templates ')

      cat('\nScore information\n')
      print(
         data.frame(
            'min.score'=sapply(object@scores, function(x) min(round(x$score, 2))), 
            'max.score'=sapply(object@scores, function(x) max(round(x$score, 2))), 
            'n.scores'=sapply(object@scores, function(x) nrow(x))
         )
      )
   }
)


setMethod('show', 'detectionList', 
   function(object) {

      cat('\nA \"detectionList\" object\n')

      cat('\nBased on survey file: ', object@survey.name, '\n')

      cat('\nand ', length(object@templates), ' templates\n')

      cat('\nDetection information\n')
      print(
         data.frame(
            'n.peaks'=sapply(object@peaks, function(x) nrow(x)), 
            'n.detections'=sapply(object@detections, function(x) nrow(x)), 
            'min.peak.score'=sapply(object@peaks, function(x) min(x$score)), 
            'max.peak.score'=sapply(object@peaks, function(x) max(x$score)), 
            'min.detection.score'=sapply(object@detections, function(x) if(length(x$score) > 0) min(x$score) else NA), 
            'max.detection.score'=sapply(object@detections, function(x) if(length(x$score) > 0) max(x$score) else NA)
         )
      )

   }
)

