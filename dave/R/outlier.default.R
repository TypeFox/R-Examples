outlier.default <-
function(veg,thresh,y,...) {
    o.outlier<- outly(veg,thresh,y)
    o.outlier$call<- match.call()
    cat("Call:\n")
    class(o.outlier)<- "outlier"
    print(o.outlier$call)
    o.outlier
    cat("\n")
    cat("Threshold value:                     ",o.outlier$threshold,"\n") 
    o.outlier
    cat("Original number of releves, species: ",o.outlier$olddim,"\n") 
    o.outlier
    cat("Remaining number of releves, species:",o.outlier$newdim,"\n") 
    o.outlier
    }
