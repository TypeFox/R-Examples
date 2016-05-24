#Interlab object

#Wish list: An 'assigned' function to extract or update assigned values 
#for a study

construct.ilab<-function(org, item, measurand, x, u, df, k, U, U.lower, U.upper, distrib=NULL, distrib.pars=NULL, study=NA, title=NA, p=0.95, ...) {

        rv<-list()
        
        rv$title <- title
        
        rv$subset <- NA
        
        L <- length(x)
        org<-rep(org, length.out=L)
        if(is.character(org)) org <- factor(org)
        item<-rep(item, length.out=L)
        measurand<-rep(measurand, length.out=L)
        study<-rep(study, length.out=L)
        
        l. <- as.data.frame(list(...))
        
        if(missing(df)) df <- rep(NA, L)

        if( !missing(U) ) {
                if( is.factor(U) ) U <- as.character(U)
                if( is.character(U) ) {
                        U.l <- U.r <- U.lower <- U.upper <- rep(NA, L)
                        #Form 1: range +-a - +-b
                        AtoB <- grep("[-+.0-9]+ *- *[-+.0-9]+",U)
                        U.l[AtoB] <- as.numeric(gsub("([-+]?[.0-9]+) *- *[-+]?[.0-9]+","\\1", U[AtoB]))
                        U.r[AtoB] <- as.numeric(gsub("[-+]?[.0-9]+ *- *([-+]?[.0-9]+)","\\1", U[AtoB]))
                        U.lower[AtoB] <- x[AtoB] - pmin(U.l[AtoB], U.r[AtoB])
                        U.upper[AtoB] <-  pmax(U.l[AtoB], U.r[AtoB]) - x[AtoB]
                        if( any( c(U.lower[AtoB], U.upper[AtoB]) <0 ) ) 
                                stop("Some x values outside range given by U=\"a-b\"")
                        
                        #Form 2: +-a/-+b
                        AslashB <- grep("[-+.0-9]+ */ *[-+.0-9]+",U)
                        U.l[AslashB] <- as.numeric(gsub("([-+]?[.0-9]+) */ *[-+]?[.0-9]+","\\1", U[AslashB]))
                        U.r[AslashB] <- as.numeric(gsub("[-+]?[.0-9]+ */ *([-+]?[.0-9]+)","\\1", U[AslashB]))
                        U.lower[AslashB] <- - pmin(U.l[AslashB], U.r[AslashB]) #unary - because lower must be -ve
                        U.upper[AslashB] <-  pmax(U.l[AslashB], U.r[AslashB])
                        
                        #Simple:
                        simple <- (1:L)[- c(AtoB, AslashB)]
                        U.l[simple] <- U.lower[simple] <- U.upper[simple] <- as.numeric(U[simple])
                        
                        U <- rep(NA, L)
                        U[simple] <- U.l[simple]
                        
                } else {
                        U.lower <- U.upper <- U
                }
        } else {
                if(!missing(u) && !missing(k) && missing(U.lower) && missing(U.upper) ) {
                        #No simple form for U if upper and lower bounds present
                        U <- k * u
                } else
                        U <- rep(NA, L)

        } 
        
        if(missing(U.lower)) U.lower <- if(!missing(U)) U else rep(NA, L)
        if(missing(U.upper)) U.upper <- if(!missing(U)) U else rep(NA, L)
        
        if(missing(u)) {
                if(!missing(U) && !missing(k) ) {
                        u <- U / k 
                } else
                        u <- rep(NA, L)

        }

        if(missing(k)) {
                if(!missing(U) && !missing(u) ) {
                        k <- U / u 
                } else
                        k <- rep(NA, L)

        }
        
        
        
        rv$data <- data.frame(
                        org=org, item=item, measurand=measurand, x=x, u=u, df=df, k=k, 
                        U=U, U.lower=U.lower, U.upper=U.upper, study=study
                    )
        
        l. <- list(...)
        if( length(l.) > 0) rv$data <- cbind(rv$data, as.data.frame(l.))
        
        if(!is.null(distrib) ) { #Omitted if not provided
                #Add distrib (can be vector)
                if(is.list(distrib)) {
                        rv$distrib<-distrib
                } else {
                        if(length(distrib) < L ) distrib <- rep(distrib, length.out=L)
                        rv$distrib<-as.list(distrib)
                }
                #Add any df deducible from k, p and distrib
                for(n in 1:L) {
                        if(is.na(rv$data$df[n]) && rv$distrib[[n]] %in% c("t", "t.scaled")) {
                                if(!is.na(rv$data$k[n])) rv$data$df[n] <- .get.df(rv$data$k[n], p)
                        }
                }
                #Add distrib.pars
                if(!is.null(distrib.pars)) { 
                        rv$distrib.pars<-as.list(distrib.pars)
                } else {
                        #Estimate distrib.pars from distrib etc
                        rv$distrib.pars<-list()
                        for( n in 1:L ) {
                                rv$distrib.pars[[n]]<-.get.pars(distrib[[n]], rv$data$x[n], rv$data$u[n], rv$data$df[n])
                        }
                }

        } else {
                rv$distrib<-as.list(rep(NA, L))
                rv$distrib.pars<-as.list(rep(NA, L))
        }
        
        class(rv) <- "ilab"
        return(rv)
}

#
# Adding information to ilab objects
#

# None, yet


#
# Print and plot functions for ilab objects
#

print.ilab <- function(x, ..., digits=NULL, right=FALSE) {
        maxwidth<-12L
        if(!is.na(x$title[1])) {
                for(s in x$title) cat(sprintf("%s\n", s))
        } else {
                cat("Interlaboratory study:\n")
        }
        
        if(!is.na(x$subset)) {
                cat(sprintf("Subset: %s\n", x$subset))
        }
        
        dp<-x$data
        if(!is.null(x[["distrib", exact=TRUE]]) ) {
                fdp<-function(x) { if(is.function(x)) deparse(x)[1] else paste(x) }
                distrib.labels<- as.vector( sapply(x$distrib, fdp ) )
                dp$distrib<-sub(paste("(.{",maxwidth,",",maxwidth,"})(.+)", sep=""), 
                        "\\1...",distrib.labels)
        }
        if(!is.null(x$distrib.pars)) {
                dp$distrib.pars <- vector("character", length=nrow(x$data) )
                for(nn in 1:nrow(x$data) ) {
                        dp[nn,"distrib.pars"]<-
                                paste(names(x$distrib.pars[[nn]]), 
                                        format(x$distrib.pars[[nn]], digits=digits),
                                        sep="=", collapse=", ")
                }
                
        }
        print.data.frame(dp,digits=digits, right=right, ...)

}

plot.ilab <- function(x, ...) {
        pars<-c(list(x=x), list(...) )
        do.call("kplot", pars)
}

#
# Functions to extract from ilab objects
#

subset.ilab <- function(x, subset, drop=FALSE, ...) {
    if (!missing(subset)) {
        e <- substitute(subset)
        r <- eval(e, x$data, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
        
        x$subset <- sprintf("subset(%s, %s)", deparse(substitute(x)), deparse(substitute(subset)))
                
        x$data<-x$data[r, ,drop=drop]
        if(!is.null(x$distrib)) x$distrib <- x$distrib[r]
        if(!is.null(x$distrib.pars)) x$distrib.pars <- x$distrib.pars[r]
        
    }

    return(x) 
}

'[.ilab' <- function(x, i, j) {
        
        x$subset <- sprintf("%s[%s, %s]", deparse(substitute(x)), 
                        deparse(substitute(i)), deparse(substitute(j)))
        
        x$data <- x$data[i,j, drop=FALSE]

        if( !is.null(x$distrib) ) x$distrib <- x$distrib[i]
        if( !is.null(x$distrib.pars) ) x$distrib.pars <- x$distrib.pars[i]
                
        return(x)
}

#
# Functions to combine or extend ilab objects
#

#base::rbind behaves oddly, so mask with a new generic

rbind<-function(..., deparse.level = 1) UseMethod("rbind")

# rbind.default<-base::rbind
        #This documented construct causes CRAN checks to fail as cbind uses .Internal

rbind.default <- function(..., deparse.level=1) base::rbind(..., deparse.level=deparse.level)

rbind.ilab<-function(..., deparse.level = 1) {
        ilab.list <- list(...)
        
        #check classes:
        il.classes <- sapply(ilab.list, function(x) class(x)[1])
        if(any(il.classes != "ilab")) 
                stop("All objects must be of class 'ilab'", call.=TRUE)
        
        if(length(ilab.list) == 0 ) {
                return(NULL)
        } else if(length(ilab.list) == 1) {
                return(ilab.list[[1]])
        } else {
                rv <- ilab.list[[1]]
                if(is.null(rv$distrib)) rv$distrib<-rep(NA, nrow(rv$data))
                if(is.null(rv$distrib.pars)) rv$distrib.pars<-as.list(rep(NA, nrow(rv$data)))
                
                for( i in 2:length(ilab.list) ) {
                        if(!isTRUE(all.equal(sort(names(rv)), sort(names(ilab.list[[i]])) ))) {
                                stop(sprintf("Names in %s do not match previous names.", names(ilab.list)[i]), call.=TRUE)
                        } else {
                                print(paste("Binding ", i, "\n"))
                                rv$data<-rbind(rv$data, ilab.list[[i]]$data, deparse.level=deparse.level)
                                if(is.null(ilab.list[[i]]$distrib)) ilab.list[[i]]$distrib<-rep(NA, nrow(ilab.list[[i]]$data))
                                if(is.null(ilab.list[[i]]$distrib.pars)) ilab.list[[i]]$distrib.pars<-as.list(rep(NA, nrow(ilab.list[[i]]$data)))
                                
                                rv$distrib<-c(rv$distrib, ilab.list[[i]]$distrib)
                                rv$distrib.pars<-c(rv$distrib.pars, ilab.list[[i]]$distrib.pars)
                                
                        }
                }
        }
        return(rv)
}

c.ilab<-function(..., recursive=FALSE) {
        rbind.ilab(...)
}

#cbind also behaves oddly, so mask with a new generic

cbind<-function(..., deparse.level = 1) UseMethod("cbind")

# cbind.default<-base::cbind
        #This documented construct causes CRAN checks to fail as cbind uses .Internal

cbind.default <- function(..., deparse.level=1) base::cbind(..., deparse.level=deparse.level)

cbind.ilab<-function(..., deparse.level = 1) {
        l<-list(...)
        L<-length(l)
        #Find the first ilab object (usually first because of dispatch,
        #but may be called directly)
        i.ilab <- which( sapply( l, function(x) class(x)[1] ) =="ilab")
        if(length(i.ilab) == 0) 
                        stop("Only one ilab object permitted in cbind.ilab", call.=TRUE)
        
        if(length(i.ilab) > 1) 
                        stop("cbind.ilab requires one ilab object", call.=TRUE)
        
        i.args <- (1:L)[-i.ilab]        
        
        #Check that args are atomic or data frame or list:
        args.ok <- sapply(l[i.args], is.atomic) | sapply(l[i.args], is.data.frame) 
        
        if( any(!args.ok) )
                stop("Arguments to cbind.ilab must be atomic, data frame or class 'ilab'", call.=TRUE)

        ilab<-l[[i.ilab]]
        #Check for over-sized data frames or vectors..
        for(i in i.args) {
                nm <- names(l)[i]
                if(is.null(nm)) nm <- sprintf("argument %d", i+1)
                if(length(dim(l[[i]]))>2) 
                        stop(sprintf("Number of dimensions of parameter %s exceeds 2", nm), call.=TRUE)
                if(length(dim(l[[i]]))==2) {
                        if( nrow(l[[i]])>nrow(ilab$data) ) 
                                stop(sprintf("Number of rows in %s exceeds rows in %s", nm, deparse(substitute(ilab))), call.=TRUE)
                } else {
                        if(length(l[[i]])>nrow(ilab$data)) 
                                stop(sprintf("Length of %s exceeds rows in %s", nm, deparse(substitute(ilab))), call.=TRUE)
                }
        }
        ilab$data <- do.call(base::cbind, c(list(ilab$data), l[i.args], list(deparse.level = deparse.level)))
        return(ilab)
}
