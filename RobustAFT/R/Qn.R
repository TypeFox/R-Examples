"Qn" <-
function(y, fact=0, ns=0){
       call <- match.call()
       n    <- length(y)
       if (ns < n) ns <- n
       if (fact <= 0) {
         dn <- 1
         const <- c(1,0.399,0.994,0.512,0.844,0.611,0.857,0.669,0.872)
         if (n <= 9) dn <- const[n]
         else {if (n %% 2 == 1) dn <- n/(n+1.4) else dn <- n/(n+3.8)}
         fact <- dn*2.2219
       }
       y    <- sort(y)
       f.res <- .Fortran("qn",
                   y=to.single(y),
                   n=to.integer(n),
                   scale=single(1),
                   sv=single(ns),
                   siw=integer(ns),
                   sw=single(ns),
                   work=single(ns),
                   left=integer(ns),
                   right=integer(ns),
                   weight=integer(ns),
                   Q=integer(ns),
                   P=integer(ns),
                   ns=to.integer(ns))
       scale <- fact*f.res$scale
       list(scale=scale, call=call)
}


