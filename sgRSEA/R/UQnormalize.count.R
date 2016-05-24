UQnormalize.count <-
function(X, round=TRUE, add.one=TRUE, print=FALSE){
    mat = as.matrix(X)
    cQuan = apply( mat, MARGIN=2, quantile, 0.75)
    cq = 1/cQuan*mean(cQuan)
    uqcount = sweep( mat, 2, FUN='*', cq)
    if (add.one) uqcount = uqcount+1
    if (round) uqcount = round(uqcount)
    if (print) {
        cat("before","\n")
        print(cQuan)
        cat("after","\n")
        print( apply( uqcount, MARGIN=2, quantile, 0.75) )
        cat("size factor","\n")
        print(cq)
        }
    return(uqcount)
    }
