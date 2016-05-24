setMethod("sqrt", signature(x = "PosSemDefSymmMatrix"), function(x){
            er <- eigen(x)
            d <- sqrt(er$values)
            D <- if(length(d)==1) diag(1)*d else diag(d)
            return(PosSemDefSymmMatrix(er$vectors %*% D %*% t(er$vectors)))
})

