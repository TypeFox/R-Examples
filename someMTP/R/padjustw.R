p.adjust.w <- function (p, method = c("bonferroni","holm","BHfwe","BH","BY"), n = length(p),w=NULL){
	
	if(is.null(w)) w=rep(1,length(p)) else w=w/sum(w)*n
	if(!is.null(names(w)) & !is.null(names(p))) w=w[names(p)]

    method <- match.arg(method)
    if (method == "fdr") 
        method <- "BH"  
    p0 <- p
    if (all(nna <- !is.na(p))) 
        nna <- TRUE
    p <- as.vector(p[nna])
    stopifnot(n >= length(p))
    if (n <= 1) 
        return(p0)
    p0[nna] <- switch(method, bonferroni = pmin(1, p*w*n), 
	   holm = { #Holm's (1979) weighted procedure
        o <- order(p/w)
        ro <- order(o)
		i <- cumsum(w[o])
        pmin(1, cummax((n - i + w[o]) * p[o]))[ro]
    },  BHfwe= { #Benjamini-Hochberg's (1997) weighted procedure for familywise error rate control
        o <- order(p)
        ro <- order(o)
		i <- cumsum(w[o])	
        pmin(1, cummax((n - i + w[o]) * p[o]))[ro]
    }, BH = { #Benjamini-Hochberg's (1997) weighted procedure for weighted false discovery rate control
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
		i <- cumsum(w[o[length(w):1]])[length(w):1]
        pmin(1, cummin(n/i * p[o]))[ro]
    }, BY = { #Benjamini and Yekutieli (2001) as extended by finos and salmaso (2007)
		o <- order(p, decreasing = TRUE)
        ro <- order(o)
        q <- sum(1/(cumsum(w[o[length(w):1]])))
        pmin(1, cummin(q * n/i * p[o]))[ro]
    }, none = p)
    attributes(p0)$w=w
	attributes(p0)$method=method
	p0
}
