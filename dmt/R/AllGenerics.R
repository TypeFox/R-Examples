setGeneric("getW",           function( model ){ standardGeneric ("getW")}) 
setGeneric("getPhi",         function( model ){ standardGeneric ("getPhi")}) 
setGeneric("getScore",       function( model ){ standardGeneric ("getScore")}) 
setGeneric("getModelMethod", function( model ){ standardGeneric ("getModelMethod")}) 
setGeneric("getParams",      function( model ){ standardGeneric ("getParams")}) 
setGeneric("getZ",           function( model, X, Y){ standardGeneric ("getZ")}) 

setGeneric("getWindowSize",  function( model ){ standardGeneric ("getWindowSize")}) # TODO: remove this?

#setGeneric("isEmpty",       function( model ){ standardGeneric ("isEmpty")})
#setGeneric("orderGenes",    function( model ){ standardGeneric ("orderGenes")})
#setGeneric("findModel",     function( model, name){ standardGeneric ("findModel")})