# For extracting templates from any object that contains them
# Modified: 2013 JUNE 12

getTemplates <-
function(
   object, 
   which.ones=names(object@templates)
) {

   if(class(object@templates[[1]]) == 'binTemplate') return(new('binTemplateList', templates=object@templates[which.ones]))
   if(class(object@templates[[1]]) == 'corTemplate') return(new('corTemplateList', templates=object@templates[which.ones]))

}
