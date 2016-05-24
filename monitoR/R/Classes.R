# Class definitions
# Modified: 2014 MAR 28
# Note: should switch back to slots argument at some point
# representation() is used for compatability with old version of R

setClass('binTemplate', 
#   slots=c(clip.path='character', samp.rate='integer', pt.on='matrix', pt.off='matrix', t.step='numeric', 
#                  frq.step='numeric', n.t.bins='integer', first.t.bin='numeric', n.frq.bins='integer', 
#                  duration='numeric', frq.lim='numeric', wl='integer', ovlp='integer', wn='character', 
#                  score.cutoff='numeric', comment="character"
#   )
   representation=representation(clip.path='character', samp.rate='integer', pt.on='matrix', pt.off='matrix', t.step='numeric', 
                  frq.step='numeric', n.t.bins='integer', first.t.bin='numeric', n.frq.bins='integer', 
                  duration='numeric', frq.lim='numeric', wl='integer', ovlp='integer', wn='character', 
                  score.cutoff='numeric', comment="character"
   )
)

setClass('binTemplateList', 
#   slots=c(templates='list'), 
   representation=representation(templates='list'), 
   validity=function(object) {
      if (any(sapply(object@templates, function(x) !is(x, 'binTemplate'))))
         return('templates slot is not a list of binTemplate objects')
      if (is.null(names(object@templates)))
         return('templates slot is missing names')
      return(TRUE) 
   }
)

setClass('corTemplate', 
#   slots=c(clip.path='character', samp.rate='integer', pts='matrix', t.step='numeric', 
#                  frq.step='numeric', n.t.bins='integer', first.t.bin='numeric', n.frq.bins='integer', 
#                  duration='numeric', frq.lim='numeric', wl='integer', ovlp='integer', wn='character', 
#                  score.cutoff='numeric', comment="character"
#   )
   representation=representation(clip.path='character', samp.rate='integer', pts='matrix', t.step='numeric', 
                  frq.step='numeric', n.t.bins='integer', first.t.bin='numeric', n.frq.bins='integer', 
                  duration='numeric', frq.lim='numeric', wl='integer', ovlp='integer', wn='character', 
                  score.cutoff='numeric', comment="character"
   )

)

setClass('corTemplateList', 
#   slots=c(templates='list'), 
   representation=representation(templates='list'), 
   validity=function(object) {
      if (any(sapply(object@templates, function(x) !is(x, 'corTemplate'))))
         return('templates slot is not a list of corTemplate objects')
      if (is.null(names(object@templates)))
         return('templates slot is missing names')
      return(TRUE) 
   }
)

# Create a class "union" for use in classes that may contain either type of template
setClassUnion('TemplateList', 
   c('binTemplateList', 'corTemplateList')
)
setClassUnion('Template', 
   c('binTemplate', 'corTemplate')
)

setClass('templateScores', 
#   slots=c(survey.name='character', survey='Wave', survey.data='list', templates='list', scores='list', time='character'), 
   representation=representation(survey.name='character', survey='Wave', survey.data='list', templates='list', scores='list', time='character'), 
   validity=function(object) {
      if (any(sapply(object@templates, function(x) !is(x, 'Template'))))
         return('templates slot is not a list of Template objects')
      if (is.null(names(object@templates)))
         return('templates slot is missing names')
      return(TRUE) 
   }
)

# MAY WANT TO ADD A validity CHECK FOR scores HERE
setClass('detectionList', 
#   slots=c(survey.name='character', survey='Wave', survey.data='list', templates='list', scores='list', peaks='list', detections='list'), 
   representation=representation(survey.name='character', survey='Wave', survey.data='list', templates='list', scores='list', peaks='list', detections='list'), 
   validity=function(object) {
      if (any(sapply(object@templates, function(x) !is(x, 'Template'))))
         return('templates slot is not a list of Template objects')
      if (is.null(names(object@templates)))
         return('templates slot is missing names')
      return(TRUE) 
   }
)

