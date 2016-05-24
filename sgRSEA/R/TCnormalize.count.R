TCnormalize.count <-
function(X, round=TRUE, add.one=TRUE, print=FALSE){
    mat = as.matrix(X)
    cSum = colSums(mat)
    c0 = 1/cSum*mean(cSum)
    tccount = sweep( mat, 2, FUN='*', c0)
    if (add.one) tccount = tccount+1
    if (round) tccount = round(tccount)
    if (print) {
        cat("before","\n")
        print(cSum)
        cat("after","\n")
        print( colSums(tccount) )
        cat("size factor","\n")
        print(c0)
        }
    return(tccount)
    }
