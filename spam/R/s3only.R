# This is file ../spam/R/s3only.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

eigen.spam <- function(x, ...) {
    inefficiencywarning( "This 'eigen' operation may be inefficient", prod(dim(x)))
    eigen(as.matrix(x), ...)
}

var.spam <- function(x, ...) {
    inefficiencywarning( "This 'var' operation may be inefficient", prod(dim(x)))
    var(as.matrix(x), ...)
}
