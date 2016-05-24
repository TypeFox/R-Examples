setGeneric(name = "rrecsys", def = function(data, alg, ...) standardGeneric("rrecsys"))

setGeneric(name = "colRatings", def = function(x) standardGeneric("colRatings"))

setGeneric(name = "rowRatings", def = function(x) standardGeneric("rowRatings"))

setGeneric(name = "numRatings", def = function(x) standardGeneric("numRatings"))

setGeneric(name = "sparsity", def = function(x) standardGeneric("sparsity"))

setGeneric(name = "evalRec", def = function(model, ...) standardGeneric("evalRec"))

setGeneric(name = "evalPred", def = function(model, ...) standardGeneric("evalPred"))

setGeneric(name = "evalModel", def = function(data, folds) standardGeneric("evalModel"))

setGeneric(name = "getAUC", def = function(model, ...) standardGeneric("getAUC"))

setGeneric(name = "predict", def = function(model, ...) standardGeneric("predict")) 

# setGeneric(name = "", def = function() standardGeneric("")) 
