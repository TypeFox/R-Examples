# This is file ../spam/R/s4coerce.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


# a few coercions that make sense...

#      showMethods(coerce)

setAs("spam","logical", def=function(from) {
    if(.Spam$structurebased) {
        return( as.logical(from@entries))     
    }else{
        inefficiencywarning( gettextf("This operation may be inefficient"), prod(dim(from)))
        return( as.logical(as.matrix(from)))
    }})

setAs("spam","vector", def=function(from) {
    if(.Spam$structurebased) {
        return( as.vector(from@entries))     
    }else{
        inefficiencywarning( gettextf("This operation may be inefficient"), prod(dim(from)))
        return( as.vector(as.matrix(from)))
    }})

setAs("spam","integer", def=function(from) {
    if(.Spam$structurebased) {
        return( as.integer(from@entries))     
    }else{
        inefficiencywarning( gettextf("This operation may be inefficient"), prod(dim(from)))
        return( as.integer(as.matrix(from)))
    }})

setAs("spam","matrix", def=function(from) {
    inefficiencywarning( gettextf("This operation may be inefficient"), prod(dim(from)))
    return( as.logical(as.matrix(from)))
})
setAs("spam","list", def=function(from) {
    return( triplet(from))
})
