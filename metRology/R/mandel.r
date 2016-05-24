#Mandel's t
#Mandel's k


#Calculation of mandel's h or k together with Grouped plots

#This version (v1) defines mandel.khas a generic to simplify early reshaping
#Matrix and vector methods reshape to a data frame and 
#pass the result to the data frame method.

#The default takes a matrix or data frame.
#Most other forms package the input into a data frame and then call the 
#data.frame method.


mandel.kh <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, type=c("h", "k"), method=c("classical", "robust"), n=NA, ...) {
        UseMethod("mandel.kh")
}

mandel.kh.default <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, type=c("h", "k"), method=c("classical", "robust"), n=NA, ...) {

        name.g <- .get.mandel.rowname(deparse(substitute(g)), rowname)
        
        #Unstack x if x is a vector and m is present
        if( is.vector(x) ) {
                if( !is.null(m) ) {
                        fm <- factor(m)
                        if( is.null(g) ) {
                                if( !all( (tm<-table(m)) == tm[1] ) ) {
                                        stop("g must be present if x is a vector and group sizes in m are unequal")
                                
                                } else {
                                        #Dummy grouping factor for .to.wide
                                        #assuming one observation per subject
                                        g<-factor( ave(1:length(m), m, FUN=function(x) seq(length(x))) )
                                        xw <- .to.wide(x,g, m)
                                }
                        } else {
                                xw <- .to.wide(x, g, m)
                        }
                        x <- xw[ , 2:ncol(xw), drop=FALSE]
                        g <- xw$g
                } else {
                        x <- data.frame(x=x)
                }
        } else {
                stop(sprintf("mandel.kh does not support objects of type %s", class(x)) )
        }

        mkh<-mandel.kh(x=x, g=g, m=m, na.rm=na.rm, rowname=rowname, type=type, method=method, n=n, ...)

        attr(mkh, "grouped.by") <- name.g
        return(mkh)
}

mandel.kh.array <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, type=c("h", "k"), method=c("classical", "robust"), n=NA, ...) {

        name.g <- .get.mandel.rowname(deparse(substitute(g)), rowname)

        #Convert to vector x if x is a 1-d array and m is present
        if( length(dim(as.array(x))) ==1 ) {
                x  <- as.vector(x)
        } else if( length(dim(as.array(x))) ==2 ){
                #Formally unnecessary, as a 2-d array has class 'matrix'
                x  <- as.data.frame(x)
        } else {
                stop("mandel.kh does not support arrays with more than 2 dimensons")
        }

        mkh<-mandel.kh(x=x, g=g, m=m, na.rm=na.rm, rowname=rowname, type=type, method=method, n=n, ...)

        attr(mkh, "grouped.by") <- name.g

        return(mkh)
}


mandel.kh.matrix <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, type=c("h", "k"), method=c("classical", "robust"), n=NA, ...) {

        mkh<-mandel.kh(x=as.data.frame(x), g=g, m=m, na.rm=na.rm, rowname=rowname, type=type, method=method, n=n, ...)

        name.g <- .get.mandel.rowname(deparse(substitute(g)), rowname)

        attr(mkh, "grouped.by") <- name.g
        
        return(mkh)
}

mandel.kh.data.frame <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, type=c("h", "k"), method=c("classical", "robust"), n=NA, ...) {
        #x is a vector, data frame or matrix
        #If a matrix, x is coerced to a data frame
        #Columns are taken as variables, rows as cases
        #if x is a vector and m NULL, x is coerced to a single column data frame
        #if m and g are present, x is first reshaped to
        #a data frame with columns named for unique(m)

        #if g is present, x is aggregated by g
        #if not, x is taken as a matrix of means or sd's
        #If g is missing for type==k, n MUST be present. 
        #rowname is used as a row name prefix if g is unspecified.
        
        #m is ignored in the data frame method because the colums are expected
        #to correspond to m
        
        #if method=="robust" robust analogues are calculated

        method <- match.arg(method)
        
        h <- function(y, na.rm=T) {
                #y is a vector of means
                return( (y-mean(y, na.rm=na.rm))/sd(y, na.rm=na.rm) )
        }
        
        k <- function(y,  na.rm=T) {
                #y is a vector of sd's
                y.omit<-na.omit(y)
                pooled.sd <- sqrt(sum(y.omit^2)/length(y.omit))
                return( y/pooled.sd )
        }

        #Simple robust h variant based on MASS hubers
        h.robust <- function(y, na.rm=T, ...) {
                #y is a vector of means
                y.omit<- if(na.rm) na.omit(y) else y
                H <- hubers(y.omit, ...)
                return( (y-H$mu)/H$s )
        }
        
        #Robust k variant based on AlgS
        k.robust <- function(y,  degfree, na.rm=T) {
                #y is a vector of sd's
                y.omit<- if(na.rm) na.omit(y) else y
                pooled.sd <- algS(y.omit, degfree)
                return( y/pooled.sd )
        }
        
        #Set the row name if g is missing
        name.g <- .get.mandel.rowname(deparse(substitute(g)), rowname)  
        
        if( is.null(g) ) { #Assume one group per row
                was.null.g <- TRUE
                if(is.null(rownames(x))) {
                        g <- paste(rowname, format(1:nrow(x)), sep="")
                        g<-factor(g, levels=g)
                        print(g)
                } else {
                        g <- factor(rownames(x), levels=rownames(x))
                                #preserves row naming order if any
                }
        } else {
                was.null.g <- FALSE
        }
        
        #Check that n is specified if type=="k" and g is missing or one per group
        if(type[1]=="k" && is.na(n) ) {
                #Get n for each column in x
                if( !was.null.g ) {
                        n.all <- aggregate(x, by=list(g=g), 
                                FUN=function(x) sum(!is.na(x)) )
                        if(ncol(n.all)>2) 
                                n.all <- stack(n.all[,2:ncol(n.all)])
                        else
                                names(n.all) <- c("g", "value") 
                                        #This guarantees presence of n.all$value
                        n <- median(n.all$value[n.all$value>0], na.rm=na.rm)
                }
                if(is.na(n) || n <=1) {
                        stop("n must be specified and >1 with type=='k' and only one value per group" )
                } 
        }

        if(type[1]=="h") {
                x <- aggregate(x, by=list(g=g), FUN=mean, na.rm=na.rm)
                mkh <- if(method=="robust") {
                                as.data.frame( lapply(x[,2:ncol(x), drop=FALSE], h.robust, na.rm=na.rm, ...) )
                        } else {
                                as.data.frame( lapply(x[,2:ncol(x), drop=FALSE], h, na.rm=na.rm) )
                        }
        } else if(type[1]=="k") {
                x <- aggregate(x, by=list(g=g), 
                        FUN=function(x, na.rm) {
                                if(length(x)==1) 
                                        x 
                                else sd(x, na.rm=na.rm)
                            },
                            na.rm=na.rm )
                if(method=="robust") {
                        #NB: uses n-1 as degfree
                        mkh <- as.data.frame( lapply(x[,2:ncol(x), drop=FALSE], k.robust, na.rm=na.rm, degfree=n-1) )
                        
                } else {
                        mkh <- as.data.frame( lapply(x[,2:ncol(x), drop=FALSE], k, na.rm=na.rm) )
                }
                
        } else {
                stop("type must be one of 'h' or 'k'")
        }
        
        row.names(mkh) <- as.character(x[[1]]) 
        
        attr(mkh, "mandel.type") <- type[1]

        attr(mkh, "mandel.method") <- method[1]
        
        attr(mkh, "grouped.by") <- name.g

        attr(mkh, "n") <- n
        
        class(mkh) <- c("mandel.kh", class(mkh))
        
        return(mkh)
        
}

mandel.kh.ilab <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, type=c("h", "k"), method=c("classical", "robust"), n=NA, ...) {

        name.g <- .get.mandel.rowname(deparse(substitute(g)), rowname)
        
        mkh<-mandel.kh(x=x$data$x, g=g, m=x$data$measurand, na.rm=na.rm, rowname=rowname, type=type, method=method, n=n, ...)
        
        attr(mkh, "grouped.by") <- name.g
        
        return(mkh)
        
}

#These convenience functions are generic to allow special handling for ilab grouped.by
mandel.h <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, method=c("classical", "robust"), n=NA, ...) {
        UseMethod("mandel.h")
}

mandel.k <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, method=c("classical", "robust"), n=NA, ...) {
        UseMethod("mandel.k")
}

#Defaults:
mandel.h.default <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, method=c("classical", "robust"), n=NA, ...) {
        
        mkh<-mandel.kh(x=x, g=g, m=m, na.rm=na.rm, rowname=rowname, type="h", method=method, n=n, ...)

        name.g <- .get.mandel.rowname(deparse(substitute(g)), rowname)

        attr(mkh, "grouped.by") <- name.g

        return(mkh)
}

mandel.k.default <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, method=c("classical", "robust"), n=NA, ...) {
        
        mkh<-mandel.kh(x=x, g=g, m=m, na.rm=na.rm, rowname=rowname, type="k", method=method, n=n, ...)

        name.g <- .get.mandel.rowname(deparse(substitute(g)), rowname)

        attr(mkh, "grouped.by") <- name.g

        return(mkh)
}


#ilab methods
#
mandel.h.ilab <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, method=c("classical", "robust"), n=NA, ...) {

        if(missing(g)) {
                g<-x$data$org
                name.g <- if(is.null(rowname)) 
                                "Organisation" 
                        else 
                                rowname
        } else {
                name.g <- if(is.null(rowname)) 
                                deparse(substitute(g))
                        else 
                                rowname
        }
        
        if(missing(m)||is.null(m)) m <- x$data$measurand
        
        mkh<-mandel.kh(x=x$data$x, g=g, m=m, na.rm=na.rm, rowname=rowname, type="h", method=method, n=n, ...)
        
        attr(mkh, "grouped.by") <- name.g
        
        return(mkh)
}

mandel.k.ilab <- function(x, g=NULL, m=NULL, na.rm=T, rowname=NULL, method=c("classical", "robust"), n=NA, ...) {

        if(missing(g)) {
                g<-x$data$org
                name.g <- if(is.null(rowname)) 
                                "Organisation" 
                        else 
                                rowname
        } else {
                name.g <- if(is.null(rowname)) 
                                deparse(substitute(g))
                        else 
                                rowname
        }
        
        if(missing(m)||is.null(m)) m <- x$data$measurand

        mkh<-mandel.kh(x=x$data$x, g=g, m=m, na.rm=na.rm, rowname=rowname, type="k", method=method, n=n, ...)
        
        attr(mkh, "grouped.by") <- name.g
        
        return(mkh)
}



#
# Plot method
#

plot.mandel.kh <- function(x, probs=c(0.95, 0.99), main, 
                                xlab=attr(x, "grouped.by"), ylab= attr(x, "mandel.type") ,
                                ylim=NULL, las=1, axes=TRUE, cex.axis=1, frame.plot = axes,
                                lwd=1, lty=1,col=par("col"), 
                                col.ind=1, lty.ind=c(2,1), lwd.ind=1,
                                separators=TRUE, col.sep="lightgrey", lwd.sep=1, lty.sep=1,
                                zero.line=TRUE, lwd.zero=1, col.zero=1, lty.zero=1,
                                p.adjust="none", ...) {
                
                if(missing(main) ) 
                        main <- paste(  deparse(substitute(x)), " - Mandel's", 
                                        attr(x, "mandel.type"), 
                                        if(attr(x, "mandel.method") == "robust") "(Robust variant)" 
                                     )
                
                ni<-ncol(x)
                #       #Number of items 
                
                ng <- nrow(x)
                #       #Number of groups

                mids <- gplot(x, main=main, xlab=xlab, ylab=ylab,
                                ylim=ylim, las=las, axes=axes, cex.axis=cex.axis, 
                                frame.plot=frame.plot, lwd=lwd, lty=lty, col=col,
                                separators=separators, col.sep=col.sep, 
                                lwd.sep=lwd.sep, lty.sep=lty.sep,
                                zero.line=zero.line, lwd.zero=lwd.zero, 
                                col.zero=col.zero, lty.zero=lty.zero, ...)
                
                
                if( !is.na(probs[1]) ) {
                        if(p.adjust != "none" ) {
                                probs <- 1 - p.adjust(1-probs, method=p.adjust, n = ng * ncol(x))
                        }
                        if(attr(x, "mandel.type") == "h" ) {
                                #Mandel's h
                                #Use 2-sided intervals
                                probs <- 1 - (1 - probs)/2
                                probs <- c(probs, 1-probs)
                                if(attr(x, "mandel.method") == "classical" ) {
                                        abline(h=qmandelh(probs, ng), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                } else {
                                        #Robust: use qnorm indicators
                                        abline(h=qnorm(probs), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                }
                        } else {
                                #Mandel's k
                                if(attr(x, "mandel.method") == "classical" ) {
                                        abline(h=qmandelk(probs, ng, attr(x, "n")), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                } else {
                                        #Robust: use f(n-1, Inf)
                                        abline(h=sqrt(qf(probs, attr(x, "n")-1, Inf)), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                }
                        }
                
                }
                
                return(invisible(mids))
}


barplot.mandel.kh <- function(height, probs=c(0.95, 0.99), main, 
                                xlab=attr(height, "grouped.by"), ylab=attr(height, "mandel.type"),
                                separators=TRUE, zero.line=TRUE, ylim,  p.adjust="none", frame.plot = TRUE,
                                ... , 
                                col.ind=1, lty.ind=c(2,1), lwd.ind=1, 
                                col.sep="lightgrey", lwd.sep=1, lty.sep=1,
                                lwd.zero=1, col.zero=1, lty.zero=1) {
                
                if(missing(main) ) 
                        main <- paste(  deparse(substitute(height)), " - Mandel's", 
                                        attr(height, "mandel.type"), 
                                        if(attr(height, "mandel.method") == "robust") "(Robust variant)" 
                                     )
                
                ng <- nrow(height)
                        #Number of groups
                        
                
                if(missing(ylim)) ylim <- range(pretty(c(0, na.omit(stack(height))$values)))
                
                mids <- barplot(t(as.matrix(height)), beside=TRUE, 
                        ylim=ylim, main=main, xlab=xlab, ylab=ylab,   ...)
                        
                if(separators) {
                        mid.max<-mids[nrow(mids), ]
                        abline(v=c(0.5, mid.max+1), col=col.sep, lty=lty.sep, lwd=lwd.sep)
                }
                if(zero.line) abline(h=0, col=col.zero, lwd=lwd.zero, lty=lty.zero)
                
                if(frame.plot) box()
                
                if( !is.na(probs[1]) ) {
                        if(p.adjust != "none" ) {
                                probs <- 1 - p.adjust(1-probs, method=p.adjust, n = ng * ncol(height))
                        }
                        if(attr(height, "mandel.type") == "h" ) {
                                #Mandel's h
                                #Use 2-sided intervals
                                probs <- 1 - (1 - probs)/2
                                probs <- c(probs, 1-probs)
                                if(attr(height, "mandel.method") == "classical" ) {
                                        abline(h=qmandelh(probs, ng), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                } else {
                                        #Robust: use qnorm indicators
                                        abline(h=qnorm(probs), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                }
                        } else {
                                #Mandel's k
                                if(attr(height, "mandel.method") == "classical" ) {
                                        abline(h=qmandelk(probs, ng, attr(height, "n")), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                } else {
                                        #Robust: use f(n-1, Inf)
                                        abline(h=qf(probs, attr(height, "n")-1, Inf), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                }
                        }
                
                }
                
                return(invisible(mids))
}

boxplot.mandel.kh <- function(x, probs=c(0.95, 0.99), main,  
                                xlab=attr(x, "grouped.by"), ylab=attr(x, "mandel.type"),
                                separators=FALSE, zero.line=TRUE, ylim,  p.adjust="none", 
                                frame.plot = TRUE, horizontal=FALSE, at,
                                ... , 
                                col.ind=1, lty.ind=c(2,1), lwd.ind=1, 
                                col.sep="lightgrey", lwd.sep=1, lty.sep=1,
                                lwd.zero=1, col.zero=1, lty.zero=1,
                                outlier.labels=row.names(x), cex.lab=0.7, col.lab=1, 
                                adj=NULL, pos=NULL, srt=0 ) {
                
                if(missing(main) ) 
                        main <- paste(  deparse(substitute(x)), " - Mandel's", 
                                        attr(x, "mandel.type"), 
                                        if(attr(x, "mandel.method") == "robust") "(Robust variant)" 
                                     )
                
                ng <- nrow(x)
                        #Number of groups
                        
                if(missing(at)) at <- 1:ncol(x)
                if(missing(ylim)) ylim <- range(pretty(c(0, na.omit(stack(x))$values)))
                
                bx <- boxplot(as.matrix(x),  horizontal=horizontal, at=at,
                        ylim=ylim, main=main, xlab=xlab, ylab=ylab,   ...)
                        
                if(separators) {
                        if(length(at)> 1 ) {
                                offset.at <- diff(at[1:2])/2
                                sep.at <-c(at[1]-offset.at, at[1]+offset.at, at[-1]+diff(at)/2)
                        } else {
                                sep.at <- at+c(-0.5,0.5)
                        }
                        if(horizontal) 
                                abline(h=sep.at, col=col.sep, lty=lty.sep, lwd=lwd.sep)
                        else
                                abline(v=sep.at, col=col.sep, lty=lty.sep, lwd=lwd.sep)
                }
                if(zero.line) abline(h=0, col=col.zero, lwd=lwd.zero, lty=lty.zero)
                
                if(frame.plot) box()
                
                if( !is.na(probs[1]) ) {
                        if(p.adjust != "none" ) {
                                probs <- 1 - p.adjust(1-probs, method=p.adjust, n = ng * ncol(x))
                        }
                        if(attr(x, "mandel.type") == "h" ) {
                                #Mandel's h
                                #Use 2-sided intervals
                                probs <- 1 - (1 - probs)/2
                                probs <- c(probs, 1-probs)
                                if(attr(x, "mandel.method") == "classical" ) {
                                        if(horizontal) 
                                                abline(v=qmandelh(probs, ng), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                        else
                                                abline(h=qmandelh(probs, ng), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                } else {
                                        #Robust: use qnorm indicators
                                        if(horizontal) 
                                                abline(v=qnorm(probs), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                        else
                                                abline(h=qnorm(probs), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                }
                        } else {
                                #Mandel's k
                                if(attr(x, "mandel.method") == "classical" ) {
                                        if(horizontal) 
                                                abline(v=qmandelk(probs, ng, attr(x, "n")), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                        else
                                                abline(h=qmandelk(probs, ng, attr(x, "n")), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                } else {
                                        #Robust: use f(n-1, Inf)
                                        if(horizontal) 
                                                abline(v=qf(probs, attr(x, "n")-1, Inf), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                        else
                                                abline(h=qf(probs, attr(x, "n")-1, Inf), lty=lty.ind, col=col.ind, lwd=lwd.ind)
                                }
                        }
                
                }
                
                if( ifelse(is.logical(outlier.labels[1]),outlier.labels[1], !is.na(outlier.labels[1])  ) ) {
                        if(is.logical(outlier.labels[1])) {
                                #Names not specified and labelling is required by TRUE: get labels
                                outlier.labels <- row.names(x) 
                        }
                        #Now go find all those outlier locations in the grouped data:
                        out.index <- rep(NA, length(bx$out))
                        for(i in 1:length(bx$out)) {
                                out.index[i] <- which.min( abs( x[,bx$group[i]] - bx$out[i] ) ) 
                        }
                        if(is.null(pos) && is.null(adj)) pos <- 4
                        
                        if(horizontal) 
                                text(bx$out, at[bx$group], outlier.labels[out.index], 
                                        cex=cex.lab, col=col.lab, pos=pos, adj=adj, srt=srt)
                        else
                                text(at[bx$group], bx$out, outlier.labels[out.index], 
                                        cex=cex.lab, col=col.lab, pos=pos, adj=adj, srt=srt)
                }
                
                return(invisible(bx))
}



#Mandel statistic quantiles and probabilities
#Mandel's k
#Mandel's k has minimum 0 and maximum sqrt(g), where g is the number of groups

qmandelk <- function(p, g, n, lower.tail = TRUE, log.p = FALSE) {
        # p: Probability
        # g: Number of groups (labs (p), in ISO 17025),
        # n: number of replicates per group
        # Principle: k^2 / g is distributed as Beta((n-1)/2, (g-1)(n-1)/2)
        # This gives qk as sqrt( g * qbeta( (n-1)/2, (g-1)*(n-1)/2))

        sqrt( g * qbeta( p, (n-1)/2, (g-1)*(n-1)/2, lower.tail=lower.tail, log.p=log.p) )

}

pmandelk <- function(q, g, n, lower.tail = TRUE, log.p = FALSE) {
        # q: Quantile
        # g: Number of groups (labs (p), in ISO 17025)
        # n: number of replicates per group
        # Principle: k^2 / g is distributed as Beta((n-1)/2, (g-1)(n-1)/2)
        # This gives pk as pbeta( q^2 / g, (n-1)/2, (g-1)*(n-1)/2)

        pbeta( q^2 / g, (n-1)/2, (g-1)*(n-1)/2, lower.tail=lower.tail, log.p=log.p)     

}

dmandelk <- function(x, g, n, log = FALSE) {
        # x: Quantile
        # g: Number of groups (labs (p), in ISO 17025)
        # n: number of replicates per group
        # Principle: k^2 / g is distributed as Beta((n-1)/2, (g-1)(n-1)/2)
        # This gives dk as 2k/g dbeta( q^2 / g, (n-1)/2, (g-1)*(n-1)/2)
        2 * x * dbeta( x^2 / g, (n-1)/2, (g-1)*(n-1)/2, log = FALSE) / g        
}


rmandelk <- function(B, g, n) {
        # B: number required (note B and not n, as n is usd in mandelk)
        # g: Number of groups (labs (p), in ISO 17025),
        # n: number of replicates per group
        # Principle: k^2 / g is distributed as Beta((n-1)/2, (g-1)(n-1)/2)
        # This gives qk as sqrt( g * qbeta( (n-1)/2, (g-1)*(n-1)/2))

        sqrt( g * rbeta( B, (n-1)/2, (g-1)*(n-1)/2) )

}


#Mandel's h distribution functions

qmandelh <- function(p, g, lower.tail = TRUE, log.p = FALSE) {
        # p: Probability
        # g: Number of groups (labs (p), in ISO 17025)
        #Principle: (1+h*sqrt(g)/(g-1))/2 is distributed as Beta((g-2)/2, (g-2)/2)
        
        ((g-1)/sqrt(g))*(2*qbeta(p, (g-2)/2, (g-2)/2, lower.tail = TRUE, log.p = FALSE)-1)      
}

pmandelh <- function(q, g, lower.tail = TRUE, log.p = FALSE) {
        # q: Quantile
        # g: Number of groups (labs (p), in ISO 17025)
        #Principle: (1+h*sqrt(g)/(g-1))/2 is distributed as Beta((g-2)/2, (g-2)/2)
        #This gives  ph as pbeta( (1+h*sqrt(g)/(g-1))/2, (g-2)/2, (g-2)/2 )
        pbeta( (1+q*sqrt(g)/(g-1))/2, (g-2)/2, (g-2)/2, lower.tail = TRUE, log.p = FALSE)       
}


dmandelh <- function(x, g, log = FALSE) {
        # x: xuantile
        # g: Number of groups (labs (p), in ISO 17025)
        #Principle: (1+h*sqrt(g)/(g-1))/2 is distributed as Beta((g-2)/2, (g-2)/2)
        #This gives  dh as 
        dbeta( (1+x*sqrt(g)/(g-1))/2, (g-2)/2, (g-2)/2, log = FALSE) / (2*(g-1)/sqrt(g))        
}

rmandelh <- function(B, g) {
        # B: number required (note B and not n, as n is usd in mandelk)
        # g: Number of groups (labs (p), in ISO 17025)
        #Principle: (1+h*sqrt(g)/(g-1))/2 is distributed as Beta((g-2)/2, (g-2)/2)
        
        ((g-1)/sqrt(g))*(2*rbeta(B, (g-2)/2, (g-2)/2)-1)        
}



#
# Utility functions
#

.get.mandel.rowname <- function(g, rowname=NULL) {
        rv <- if(g=="NULL" | g=="") {
                if(is.null(rowname)) 
                        "Row"
                else
                        rowname
        } else {
                g
        }
        
        return(rv)
}


.to.wide <- function(x, g, m) {
        #x is a vector grouped by g and m
        #where mandel's k will be calculated for 
        #all m and grouped by g. In interlab studies, g is usually the lab.
        
        #Generally there are n>=1 for all g:m
        #but the group size may not be constant
        
        #to.wide returns a data frame with colums named by m and a further grouping column g
        #that is, in a form usable by mandel.hk
        
        #reshape ought to be able to do this, but writing 
        #a custom routine was quicker than working out why 
        #reshape wasn't doing what I wanted. 
        
        if(!is.factor(g)) g <- factor(g)
        if(!is.factor(m)) g <- factor(m)
        
        #Count size of each group:
        n.per.g <- tapply(x, g:m, length)
        n <- max(n.per.g, na.rm=TRUE)
        
        d <- data.frame(g = factor(rep(levels(g), each=n), levels=levels(g)))
        row.names(d) <- paste(d$g, rep(1:n, along=d$g), sep=":")
        
        n.x <- paste(g, ave(x, g, m, FUN=function(x) 1:length(x)), sep=":")

        for(nn in levels(m) ) {
                subx <- x[m==nn]
                names(subx) <- n.x[m==nn]
                d[[nn]] <- subx[row.names(d)] 
        }
        
        return(d)
}

