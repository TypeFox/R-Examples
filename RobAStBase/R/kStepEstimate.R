###############################################################################
## Functions and methods for "kStepEstimate" classes and subclasses
###############################################################################

setMethod("pICList", "kStepEstimate", function(object) object@pICList)
setMethod("ICList", "kStepEstimate", function(object) object@ICList)
setMethod("start", "kStepEstimate", function(x) x@start)
setMethod("startval", "kStepEstimate", function(object) object@startval)
setMethod("ustartval", "kStepEstimate", function(object) object@ustartval)
setMethod("ksteps", "kStepEstimate", function(object, diff = FALSE) {
     mm <- cbind(object@startval,object@ksteps)
     rownames(mm) <- rownames(object@ksteps)
     if(diff){
        return(t(apply(mm,1,diff)))
     }
     colnames(mm) <- paste((1:ncol(mm))-1)
     return(mm)
})
setMethod("uksteps", "kStepEstimate", function(object, diff = FALSE) {
     mm <- cbind(object@ustartval,object@uksteps)
     rownames(mm) <- rownames(object@uksteps)
     if(diff){
        return(t(apply(mm,1,diff)))
     }
     colnames(mm) <- paste((1:ncol(mm))-1)
     return(mm)
})
