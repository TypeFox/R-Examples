.construct.loc.est<-function(x, u, df=NA, xi, ui, dfi=rep(NA, length(xi)), 
                u.eff, w, method, method.details) {
        if(missing(w)) w <- rep(1, length(x))
        rv<-list( x=x, u=u, df=df, xi=xi, ui=ui, dfi=dfi, u.eff=u.eff, w=w, 
                        method=method, method.details=method.details)
        class(rv) <- "loc.est"
        return(rv)
}

print.loc.est <- function(x, ...) {
        cat(sprintf("\nLocation estimate: Method=%s\n\n", x$method))
        est<-matrix(c(x$x, x$u), ncol=2)
        dimnames(est) <- list("estimate", c("Value", "u"))
        print(est)
        cat("\n")
}

