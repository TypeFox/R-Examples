setMethod("solve", signature(a = "ANY", b = "ANY"), function(a,b, 
             generalized = getdistrOption("use.generalized.inverse.by.default"),
          tol = .Machine$double.eps, ...) {
                 if(!generalized) return(base::solve(a,b, tol = tol, ...))
                 else if(is(try({
                            ab <- base::solve(a,b, tol = tol, ...)
                            if(missing(b))
                                 dimnames(ab) <-  rev(dimnames(a))
                            else names(ab) <-  colnames(a)
                            return(ab)
                            }, silent = TRUE), "try-error")){
             if (!missing(b))
                if(!(length(b)==nrow(a))) stop("non-conformable arguments")
             a.m <- MASS::ginv(a)
             dimnames(a.m) <- rev(dimnames(a))             
             if (missing(b)) return(a.m) 
             else return(a.m %*% b)
             }})

setMethod("solve", signature(a = "PosSemDefSymmMatrix", b = "ANY"), 
           function(a,b, 
               generalized = getdistrOption("use.generalized.inverse.by.default"), 
               tol = .Machine$double.eps, ...){
          if(!generalized) return(base::solve(a,b, tol = tol, ...))
          else{
            er <- eigen(a)
            d1 <- er$values
            d <- 1/d1[d1 > tol]
            ev <- er$vectors[,d1 > tol]
            A <- if (length(d)) ev %*% (t(ev)*d) else 0*a
            dimnames(A) <- dimnames(a)
            if(missing(b)) return(A)
            else return(A%*%b)}   
})

setMethod("solve", signature(a = "PosDefSymmMatrix", b = "ANY"), function(a,b, ...){
base::solve(a,b, ...)})
