# Summary methods for objects
# Modified: 2013 JUNE 14

setMethod('summary', 'TemplateList', function(object) {show(object)})
setMethod('summary', 'templateScores', function(object) {show(object)})
setMethod('summary', 'detectionList', function(object) {show(object)})
