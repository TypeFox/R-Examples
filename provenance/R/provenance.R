#' Read a .csv file with continuous (detrital zircon) data
#'
#' Reads a data table containing continuous data (e.g. detrital zircon
#' ages)
#' @param fname the path of a .csv file with the input data,
#' arranged in columns.
#' @param errorfile the (optional) path of a .csv file with the
#' standard errors of the input data, arranged by column in the same
#' order as \code{fname}. Must be specified if the data are to be
#' compared with the Sircombe-Hazelton dissimilarity.
#' @param method an optional string specifying the dissimilarity
#' measure which should be used for comparing this with other
#' datasets. Should be one of either \code{"KS"} (for
#' Kolmogorov-Smirnov) or \code{"SH"} (for Sircombe and Hazelton). If
#' \code{method = "SH"}, then \code{errorfile} should be specified. If
#' \code{method = "SH"} and \code{errorfile} is unspecified, then the
#' program will default back to the Kolmogorov-Smirnov dissimilarity.
#' @param xlab an optional string specifying the nature and units of
#' the data.  This string is used to label kernel density estimates.
#' @param colmap an optional string with the name of one of R's
#' built-in colour palettes (e.g., heat.colors, terrain.colors,
#' topo.colors, cm.colors), which are to be used for plotting the data.
#' @return an object of class \code{distributional}, i.e. a list with the
#' following items:
#' 
#' \code{x}: a named list of vectors containing the numerical data for each sample
#' 
#' \code{err}: an (optional) named list of vectors containing the standard errors of \code{x}
#'
#' \code{method}: either "KS" (for Kolmogorov-Smirnov) or "SH" (for Sircombe Hazelton)
#' 
#' \code{breaks}: a vector with the locations of the histogram bin edges
#' 
#' \code{xlab}: a string containing the label to be given to the x-axis on all plots
#' @examples
#' agefile <- system.file("DZ.csv",package="provenance")
#' errfile <- system.file("DZerr.csv",package="provenance")
#' DZ <- read.distributional(agefile,errfile)
#' plot(KDE(DZ$x$N1))
#' @export
read.distributional <- function(fname,errorfile=NA,method="KS",xlab="age [Ma]",colmap='rainbow') {
    out <- list()
    out$name <- basename(substr(fname,1,nchar(fname)-4))
    if (method=="SH" & is.na(errorfile)) method <- "KS"
    class(out) <- "distributional"
    out$method <- method
    out$x <- list()
    out$err <- list()
    out$colmap <- colmap
    dat <- utils::read.csv(fname,header=TRUE)
    ns = length(dat)
    for (i in 1:ns){
        out$x[[names(dat)[i]]] = dat[!is.na(dat[,i]),i]
    }
    if (!is.na(errorfile)){
        err <- utils::read.csv(errorfile,header=TRUE)
        for (i in 1:ns) {
            out$err[[names(dat)[i]]] = dat[!is.na(err[,i]),i]
        }
    }
    d <- unlist(out$x)
    ng <- length(d) # number of grains
    nb <- log(ng/ns,base=2)+1
    out$breaks <- seq(min(d),max(d),length.out=nb+1)
    out$xlab <- xlab
    return(out)
}

#' Read a .csv file with categorical data
#'
#' Reads a data table containing categorical data (e.g. petrographic,
#' heavy mineral or geochemical data)
#'
#' @param fname a string with the path to the .csv file
#' @param method either "bray" (for the Bray-Curtis distance) or
#' "aitchison" (for Aitchison's central logratio distance). If
#' omitted, the function defaults to 'aitchison', unless there are
#' zeros present in the data.
#' @param colmap an optional string with the name of one of R's
#' built-in colour palettes (e.g., heat.colors, terrain.colors,
#' topo.colors, cm.colors), which are to be used for plotting the data.
#' @return an object of class \code{compositional}, i.e. a list with the
#' following items:
#' 
#' \code{x}: a data frame with the samples as rows and the categories as columns
#'
#' \code{method}: either "aitchison" (for Aitchison's centred logratio
#' distance) or "bray" (for the Bray-Curtis distance)
#' @examples
#' fname <- system.file("Major.csv",package="provenance")
#' Major <- read.compositional(fname)
#' plot(PCA(Major))
#' @export
read.compositional <- function(fname,method=NULL,colmap='rainbow') {
    out <- list()
    out$name <- basename(substr(fname,1,nchar(fname)-4))
    class(out) <- "compositional"
    out$x <- utils::read.csv(fname,header=TRUE,row.names=1)
    if (is.null(method)){
        if (any(out$x==0)) { method <- "bray" }
        else { method <- "aitchison" }
    }
    out$method <- method
    out$colmap <- colmap
    if (any(out$x==0) & method=="aitchison"){
        stop(paste("This dataset contains zeros and is",
                   "incompatible with the 'aitchison' distance"))
    }
    return(out)
}
#' Read a .csv file with mineral and rock densities
#'
#' Reads a data table containing densities to be used for
#' hydraulic sorting corrections (minsorting and srd functions)
#'
#' @param fname a string with the path to the .csv file
#' @return a vector with mineral and rock densities
#' @examples
#' data(Namib,densities)
#' N8 <- subset(Namib$HM,select="N8")
#' distribution <- minsorting(N8,densities,phi=2,sigmaphi=1,medium="air",by=0.05)
#' plot(distribution)
#' @export
read.densities <- function(fname){
    return(utils::read.csv(fname,header=TRUE))
}

#' create a \code{data.frame} object
#'
#' Convert an object of class \code{compositional} to a
#' \code{data.frame} for use in the \code{robCompositions} package
#'
#' @param x an object of class \code{compositional}
#' @param ... optional arguments to be passed on to the generic function
#' @return a \code{data.frame}
#' @examples
#' data(Namib)
#' qfl <- ternary(Namib$PT,c('Q'),c('KF','P'),c('Lm','Lv','Ls'))
#' plot(qfl,type="QFL.dickinson")
#' qfl.frame <- as.data.frame(qfl)
#' ## uncomment the next two lines to plot an error
#' ## ellipse using the robCompositions package:
#' # library(robCompositions)
#' # pca <- pcaCoDa(qfl.frame)
#' # plot(pca,xlabs=rownames(qfl.frame))
#' @export
as.data.frame.compositional <- function(x,...){
    nc <- ncol(as.matrix(x$x))
    if (nc==3) out <- data.frame(x$x[,c(2,3,1)],...)
    if (nc>3) out <- data.frame(x$x[,c(2,3,1,4:nc)],...)
    if (nc<3) out <- data.frame(x$x,...)
    return(out)
}

#' create an \code{acomp} object
#'
#' Convert an object of class \code{compositional} to an object of
#' class \code{acomp} for use in the \code{compositions} package
#'
#' @param x an object of class \code{compositional}
#' @return a \code{data.frame}
#' @examples
#' data(Namib)
#' qfl <- ternary(Namib$PT,c('Q'),c('KF','P'),c('Lm','Lv','Ls'))
#' plot(qfl,type="QFL.dickinson")
#' qfl.acomp <- as.acomp(qfl)
#' ## uncomment the next two lines to plot an error
#' ## ellipse using the compositions package: 
#' # library(compositions)
#' # ellipses(mean(qfl.acomp),var(qfl.acomp),r=2)
#' @export
as.acomp <- function(x){
    if (!methods::is(x,"compositional")){
        stop("not an object of class compositional or ternary")
    }
    dat <- as.matrix(as.data.frame(x))
    out <- structure(dat)
    attributes(out) <- list(dim = dim(dat), dimnames=dimnames(dat), class="acomp")
    return(out)
}

as.compositional.matrix <- function(x,method=NULL,colmap='rainbow'){
    out <- list(x=NULL,method=method,colmap=colmap)
    class(out) <- "compositional"
    out$x <- as.matrix(x)
    nc <- ncol(out$x)
    if (nc==3) out$x <- out$x[,c(3,1,2)]
    if (nc>3) out$x <- out$x[,c(3,1,2,4:nc)]
    if (nc<3) out$x <- out$x
    return(out)
}

#' create a \code{compositional} object
#'
#' Convert an object of class \code{matrix}, \code{data.fram} or
#' \code{acomp} to an object of class \code{compositional}
#'
#' @param x an object of class \code{matrix}, \code{data.fram} or
#' \code{acomp}
#' @param method dissimilarity measure, either 'aitchison' for
#' Aitchison's CLR-distance or 'bray' for the Bray-Curtis distance.
#' @param colmap the colour map to be used in pie charts.
#' @return an object of class \code{compositional}
#' @examples
#' data(Namib)
#' PT.acomp <- as.acomp(Namib$PT)
#' PT.compositional <- as.compositional(PT.acomp)
#' print(Namib$PT$x - PT.compositional$x)
#' ## uncomment the following lines for an illustration of using this 
#' ## function to integrate the \code{provenance} package with \code{compositions}
#' # library(compositions)
#' # data(Glacial)
#' # a.glac <- acomp(Glacial)
#' # c.glac <- as.compositional(a.glac)
#' # summaryplot(c.glac,ncol=8)
#' @export
as.compositional <- function(x,method=NULL,colmap='rainbow'){
    if (methods::is(x,"acomp")){
        attr <- attributes(x)
        attributes(x) <- NULL
        x[!is.numeric(x)] <- NA
        y <- matrix(x,nrow=attr$dim[[1]],ncol=attr$dim[[2]],dimnames=attr$dimnames)
        print(print(attr$dim[[1]]))
        return(as.compositional.matrix(y))
    } else if (methods::is(x,"data.frame") | methods::is(x,"matrix")){
        y <- as.matrix(x)
        dimnames(y) <- dimnames(x)
        return(as.compositional.matrix(y,method,colmap))
    } else {
        stop(paste("cannot convert an object of class",class(x),
                   "into an object of class compositional"))
    }
}

#' Kolmogorov-Smirnov dissimilarity
#'
#' Returns the Kolmogorov-Smirnov dissimilarity between two samples
#'
#' @param x the first sample as a vector
#' @param y the second sample as a vector
#' @return a scalar value representing the maximum vertical distance
#' between the two cumulative distributions
#' @examples
#' data(Namib)
#' print(KS.diss(Namib$DZ$x[['N1']],Namib$DZ$x[['T8']]))
#' @export
KS.diss <- function(x,y) {
    xx = sort(x)
    cdftmp = stats::ecdf(xx)
    cdf1 = cdftmp(xx)
    xy = sort(y)
    cdftmp = stats::ecdf(xy)
    cdfEstim = cdftmp(xy)
    cdfRef = stats::approx(xx, cdf1, xy, yleft = 0, yright = 1, ties = "mean")
    dif = cdfRef$y - cdfEstim
    dif = abs(dif)
    out = max(dif)
    return(out)
}

#' Calculate the dissimilarity matrix between two \code{distributional} or
#' \code{compositional} datasets
#'
#' Calculate the dissimilarity matrix between two datasets of class
#' \code{distributional} or \code{compositional} using the Kolmogorov-Smirnov,
#' Sircombe-Hazelton, Aitchison or Bray Curtis distance
#' 
#' @param x an object of class \code{distributional} or \code{compositional}
#' @param method (optional) either "KS", "SH", "aitchison" or "bray"
#' @examples
#' data(Namib)
#' print(round(100*diss(Namib$DZ)))
#' @return an object of class \code{diss}
#' @rdname diss
#' @export
diss <- function(x,method){ UseMethod("diss",x) }
#' @rdname diss
#' @export
diss.distributional <- function(x,method=NULL) {
    if (!is.null(method)) x$method <- method
    n <- length(x$x)
    d <- mat.or.vec(n,n)
    rownames(d) <- names(x$x)
    colnames(d) <- names(x$x)
    if (x$method=="SH") c2 <- getc2(x)
    for (i in 1:n){
        for (j in 1:n){
            if (x$method=="SH"){
                d[i,j] <- SH.diss(x,i,j,c.con=c2)
            }
            if (x$method=="KS"){
                d[i,j] <- KS.diss(x$x[[i]],x$x[[j]])
            }
        }
    }
    out <- stats::as.dist(d)
    class(out) <- append("diss",class(out))
    return(out)
}
#' @rdname diss
#' @export
diss.compositional <- function(x,method=NULL){
    if (!is.null(method)) x$method <- method
    if (x$method=="aitchison"){
        out <- stats::dist(CLR(x)$x)
    } else {
        snames <- names(x)
        ns <- length(snames)
        d <- mat.or.vec(ns,ns)
        rownames(d) <- snames
        colnames(d) <- snames
        for (i in 1:ns){
            for (j in 1:ns){
                d[i,j] <- bray.diss(x$x[i,],x$x[j,])
            }
        }
        out <- stats::as.dist(d)
    }
    class(out) <- append("diss",class(out))
    return(out)
}

#' Bray-Curtis dissimilarity
#'
#' Calculates the Bray-Curtis dissimilarity between two samples
#' @param x a vector containing the first compositional sample
#' @param y a vector of length(x) containing the second compositional sample
#' @return a scalar value
#' @examples
#' data(Namib)
#' print(bray.diss(Namib$HM$x["N1",],Namib$HM$x["N2",]))
#' @export
bray.diss <- function(x,y){
    return(as.numeric(sum(abs(x-y))/sum(x+y)))
}

#' Multidimensional Scaling
#'
#' Performs classical or nonmetric Multidimensional Scaling analysis
#' @param x an object of class \code{distributional}, \code{compositional} or \code{diss}
#' @param classical boolean flag indicating whether classical (TRUE)
#' or nonmetric (FALSE) MDS should be used
#' @param ... optional arguments to be passed onto \code{diss} (if
#' \code{x} is of class \code{compositional} or \code{distributional})
#' or onto \code{cmdscale} or \code{isoMDS} (if \code{x} is of class
#' \code{dist}).
#' @return an object of class \code{MDS}, i.e. a list containing the
#' following items:
#'
#' \code{points}: a two column vector of the fitted configuration
#'
#' \code{classical}: a boolean flag indicating whether the MDS
#' configuration was obtained by classical (\code{TRUE}) or nonmetric
#' (\code{FALSE}) MDS.
#'
#' \code{diss}: the dissimilarity matrix used for the MDS analysis
#' 
#' \code{stress}: (only if \code{classical=TRUE}) the final stress
#' achieved (in percent)
#' @examples
#' data(Namib)
#' plot(MDS(Namib$Major,classical=TRUE))
#' @rdname MDS
#' @importFrom MASS isoMDS
#' @export
MDS <- function(x,...){ UseMethod("MDS",x) }
#' Multidimensional Scaling of compositional data
#'
#' @rdname MDS
#' @export
MDS.compositional <- function(x,classical=FALSE,...){
    d <- diss.compositional(x,...)
    return(MDS.diss(d,classical=classical))
}
#' Multidimensional Scaling of distributional data
#'
#' @rdname MDS
#' @export
MDS.distributional <- function(x,classical=FALSE,...){
    d <- diss.distributional(x,...)
    return(MDS.diss(d,classical=classical))
}
#' Multidimensional Scaling of a dissimilarity matrix
#'
#' @rdname MDS
#' @export
MDS.diss <- function(x,classical=FALSE,...){
    out <- list() 
    if (classical){
        out$points <- stats::cmdscale(x)
    } else {
        out <- MASS::isoMDS(d=x,...)
    }
    out$classical <- classical
    out$diss <- x
    class(out) <- "MDS"
    return(out)
}

#' Centred logratio transformation
#'
#' Calculates Aitchison's centered logratio transformation for a
#' dataset of class \code{compositional}
#' @param x an object of class \code{compositional}
#' @return a matrix of CLR coordinates
#' @examples
#' # The following code shows that applying provenance's PCA function
#' # to compositional data is equivalent to applying R's built-in
#' # princomp function to the CLR transformed data.
#' data(Namib)
#' plot(PCA(Namib$Major))
#' dev.new()
#' clrdat <- CLR(Namib$Major)$x
#' biplot(princomp(clrdat))
#' @export
CLR <- function(x){
    if (!methods::is(x,'compositional')){stop('CLR(x): x is not of class compositional.')}
    out <- x
    g <- apply(log(x$x),1,mean)
    nc <- ncol(x$x)
    gg <- matrix(rep(g,nc),ncol=nc,byrow=FALSE)
    out$x <- log(x$x) - gg
    return(out)
}

#' Principal Component Analysis
#'
#' Performs PCA of compositional data using a centred logratio distance
#' @param x an object of class \code{compositional}
#' @param ... optional arguments to R's \code{princomp function}
#' @return an object of classes \code{PCA}, which is synonymous to
#' the stats packages' \code{princomp} class.
#' @examples
#' data(Namib)
#' plot(MDS(Namib$Major,classical=TRUE))
#' dev.new()
#' plot(PCA(Namib$Major),asp=1)
#' print("This example demonstrates the equivalence of classical MDS and PCA")
#' @export
PCA <- function(x,...){
    if (!methods::is(x,'compositional')){stop('x is not of class compositional in PCA(x)')}
    clrdat <- CLR(x)
    pc <- stats::princomp(clrdat$x,...)
    class(pc) <- append("PCA",class(pc))
    return(pc)
}

#' Get a subset of distributional data
#'
#' Return a subset of provenance data according to some specified indices
#' @param x an object of class \code{distributional}
#' @param subset logical expression indicating elements or rows to keep:
#' missing values are taken as false.
#' @param select a vector of sample names
#' @param ... optional arguments for the generic subset function
#' @return an object of class \code{distributional}
#' @seealso read.distributional
#' @examples
#' data(Namib)
#' coast <- subset(Namib$HM,select=c("N1","N2","T8","T13","N12","N13"))
#' summaryplot(coast,ncol=2)
#' @export
subset.distributional <- function(x,subset=NULL,select=NULL,...){
    out <- x
    if (!is.null(subset)){
        i <- which(subset,arr.ind=TRUE)
    } else if (!is.null(select)){
        i <- which(names(x) %in% select)
    } else {
        return(out)
    }
    if (length(x$err)==length(x$x)) out$err <- x$err[i]
    out$x <- x$x[i]
    return(out)
}
#' Get a subset of compositional data
#'
#' Return a subset of provenance data according to some specified indices
#' @param x an object of class \code{compositional}
#' @param subset logical expression indicating elements or rows to keep:
#' missing values are taken as false.
#' @param select a vector of sample names.
#' @param components a vector specifying a subcomposition
#' @param ... optional arguments for the generic subset function
#' @return an object of class \code{compositional}
#' @seealso read.compositional
#' @export
subset.compositional <- function(x,subset=NULL,select=NULL,components=NULL,...){
    out <- x
    if (!is.null(subset)){
        i <- which(subset,arr.ind=TRUE)
    } else if (!is.null(select)){
        i <- which(names(x) %in% select)
    } else {
        i <- 1:length(names(x))
    }
    if (!is.null(components)){
        j <- which(colnames(x$x) %in% components,arr.ind=TRUE)
    } else {
        j <- 1:length(colnames(x$x))
    }
    out$x <- x$x[i,j]
    if (methods::is(x,"SRDcorrected")){
        out$restoration <- x$restoration[i]
        for (sname in rownames(out$x)){
            out$restoration[[sname]] <- subset(x$restoration[[sname]],select=j)
        }
    }
    return(out)
}

# returns list of dissimilarities between common items
getdisslist <- function(slist){
    dnames <- names(slist)
    lablist <- lapply(slist,function(x) names(x))
    commonlabels <- Reduce(intersect,lablist)
    for (name in dnames){
        slist[[name]] <- subset(slist[[name]],select=commonlabels)
    }
    disslist <- slist
    for (name in dnames){
        disslist[[name]] <- diss(slist[[name]])
    }
    return(disslist)
}

#' Generalised Procrustes Analysis of provenance data
#'
#' Given a number of input datasets, this function performs an MDS
#' analysis on each of these and the feeds the resulting
#' configurations into the \code{GPA()} function.
#'
#' @param ... a sequence of datasets of classes \code{distributional}
#' and \code{compositional}
#' @return an object of class \code{GPA}, i.e. a list containing the
#' following items:
#' 
#' \code{points}: a two column vector with the coordinates of the
#' group configuration
#'
#' \code{labels}: a list with the sample names
#' @author Pieter Vermeesch
#' @references Gower, J.C. (1975). Generalized Procrustes analysis,
#' Psychometrika, 40, 33-50.
#' @examples
#' data(Namib)
#' gpa <- procrustes(Namib$DZ,Namib$HM)
#' plot(gpa)
#' @importFrom MASS isoMDS
#' @seealso GPA
#' @export
procrustes <- function(...) {
    dnames <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
    slist <- list(...)
    names(slist) <- dnames
    disslist <- getdisslist(slist)
    n <- length(labels(disslist[[1]]))
    m <- length(disslist)
    X <- array(dim=c(n,2,m))
    for (i in 1:m){
        md <- MDS(disslist[[i]],FALSE)
        if (md$stress < 0.05) md <- MDS(disslist[[i]],TRUE)
        X[,,i] <- md$points
    }
    result <- GPA(X)
    out <- list()
    out$points <- result
    out$labels <- labels(disslist[[1]])
    class(out) <- "GPA"
    return(out)
}

#  based on a Wikipedia algorithm
#' Generalised Procrustes Analysis of configurations
#'
#' Given a number of (2D) configurations, this function uses a
#' combination of transformations (reflections, rotations,
#' translations and scaling) to find a 'consensus' configuration which
#' best matches all the component configurations in a least-squares
#' sense.
#' 
#' @param X a list of dissimilarity matrices
#' @param scale boolean flag indicating if the transformation should include the scaling operation
#' @return a two column vector with the coordinates of the
#' group configuration
#' @seealso procrustes
#' @export
GPA <- function(X,scale=TRUE){
    if (length(dim(X))<3) {
        return(X)
    } else if (dim(X)[3]<3){
        return(procfit(X[,,1],X[,,2])$Yrot)
    } else {
        Y <- X # initialise fitted configurations
        refconf <- X[,,1] # reference configuration
        for (j in 1:100){
            for (i in 1:dim(X)[3]){
                Y[,,i] <- procfit(refconf,X[,,i])$Yrot
            }
            meanconf <- apply(Y,c(1,2),'mean')
            misfit <- sum((refconf-meanconf)^2)
            if (misfit < 1e-10){
                break
            } else {
                refconf <- meanconf
            }
        }
        return(refconf)
    }
}

# Procrustes analysis of two configurations
# based on the 'procrustes' function of the 'vegan' package
procfit <- function (X, Y, scale=TRUE, symmetric=FALSE, ...) {
    if (nrow(X) != nrow(Y)) 
        stop("Matrices have different number of rows: ", nrow(X), 
            " and ", nrow(Y))
    if (ncol(X) < ncol(Y)) {
        warning("X has fewer axes than Y: X adjusted to comform Y\n")
        addcols <- ncol(Y) - ncol(X)
        for (i in 1:addcols) X <- cbind(X, 0)
    }
    ctrace <- function(MAT) sum(MAT^2)
    c <- 1
    if (symmetric) {
        X <- scale(X, scale = FALSE)
        Y <- scale(Y, scale = FALSE)
        X <- X/sqrt(ctrace(X))
        Y <- Y/sqrt(ctrace(Y))
    }
    xmean <- apply(X, 2, mean)
    ymean <- apply(Y, 2, mean)
    if (!symmetric) {
        X <- scale(X, scale = FALSE)
        Y <- scale(Y, scale = FALSE)
    }
    XY <- crossprod(X, Y)
    sol <- svd(XY)
    A <- sol$v %*% t(sol$u)
    if (scale) {
        c <- sum(sol$d)/ctrace(Y)
    }
    Yrot <- c * Y %*% A
    b <- xmean - c * ymean %*% A
    R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
    reslt <- list(Yrot = Yrot, X = X, ss = R2, rotation = A, 
        translation = b, scale = c, xmean = xmean, symmetric = symmetric, 
        call = match.call())
    reslt$svd <- sol
    class(reslt) <- "procrustes"
    reslt
}

# calculate the trace of a matrix
tr <- function (m){
    if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) 
        stop("m must be a square matrix")
    return(sum(diag(m)))
}

get.data.names <- function(dlist){
    out <- c()
    for (d in dlist){
        out <- c(out,d$name)
    }
    out
}

# set minimum and maximum values of a dataset
setmM <- function(x,from=NA,to=NA,log=FALSE){
    if (is.na(from)) { from <- min(x); setm <- TRUE }
    else { setm <- FALSE }
    if (is.na(to)) { to <- max(x); setM <- TRUE }
    else { setM <- FALSE }
    if (setm) {
        if (log) { from <- from/2 }
        else {
            if (2*from-to<0) {from <- 0}
            else {from <- from-(to-from)/10}
        }
    }
    if (setM) {
        if (log) { to <- 2*to }
        else { to <- to+(to-from)/10 }
    }
    return(list(m=from,M=to))
}

#' @export
names.distributional <- function(x){
    return(names(x$x))
}
#' @export
names.compositional <- function(x){
    return(rownames(x$x))
}
#' @export
names.KDEs <- function(x){
    return(names(x$kdes))
}
#' @export
names.ternary <- function(x){
    return(rownames(x$x))
}

#' Calculate the number of grains required to achieve a desired level of sampling resolution
#'
#' Returns the number of grains that need to be analysed to decrease
#' the likelihood of missing any fraction greater than a given size
#' below a given level.
#' @param f the size of the smallest resolvable fraction (0<f<1)
#' @param p the probability that all n grains in the sample have missed
#' at least one fraction of size f
#' @param n, the number of grains in the sample
#' @return the number of grains needed to reduce the chance of missing
#' at least one fraction f of the total population to less than p
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' # number of grains required to be 99% that no fraction greater than 5% was missed:
#' print(get.n(0.01))
#' # number of grains required to be 90% that no fraction greater than 10% was missed:
#' print(get.n(p=0.1,f=0.1))
#' @export
get.n <- function(p=0.05,f=0.05){
    n <- 1
    while(T){
        pp <- get.p(n,f)
        if (pp<p){ break }
        else {n <- n+1}
    }
    return(n)
}

#' Calculate the probability of missing a given population fraction
#'
#' For a given sample size, returns the likelihood of missing any
#' fraction greater than a given size
#' @param n the number of grains in the detrital sample
#' @param f the size of the smallest resolvable fraction (0<f<1)
#' @return the probability that all n grains in the sample have missed
#' at least one fraction of size f
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' print(get.p(60))
#' print(get.p(117))
#' @export
get.p <- function(n,f=0.05){
    p <- 0
    M <- 1/f
    for (i in 1:M){
        p <- p + (-1)^(i-1) * choose(M,i)*(1-i*f)^n
    }
    return(p)
}

#' Calculate the largest fraction that is likely to be missed
#'
#' For a given sample size, returns the largest fraction which has
#' been sampled with p x 100 % likelihood.
#' @param n the number of grains in the detrital sample
#' @param p the required level of confidence
#' @return the largest fraction that is sampled with at least 100 x p%
#' certainty
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' print(get.f(60))
#' print(get.f(117))
#' @export
get.f <- function(n,p=0.05){
    fmin <- 0
    fmax <- 1
    for (i in 1:100){
        f <- (fmax+fmin)/2
        if (get.p(n,f)<p) { fmax <- f }
        else { fmin <- f }
    }
    return((fmin+fmax)/2)
}

get.densities <- function(X,dtable){
    if (!methods::is(X,"compositional")) stop("input is not of class compositional")
    minerals <- colnames(X$x)
    i <- which(colnames(dtable) %in% colnames(X$x), arr.ind=TRUE)
    return(dtable[i])
}

ndim <- function(X){
    return(length(dim(X)))
}

sumcols <- function(X,x){
    if (length(x)>1 & ndim(X[,x])>0) # >1 class, >1 sample
        out <- apply(X[,x],1,sum)
    if (length(x)>1 & ndim(X[,x])==0) # >1 class, 1 sample
        out <- sum(X[,x])
    if (length(x)==1 & ndim(X[,x])>0) # 1 class, >1 sample
        out <- sum(X[,x])
    if (length(x)==1 & ndim(X[,x])==0) # 1 class, 1 sample
        out <- X[,x]
    names(out) <- rownames(X)
    return(out)
}

# X is a vector of strings
sumlabels <- function(X){
    out <- X[1]
    n <- length(X)
    if (n==1) return(out)
    for (i in 2:length(X)){
        out <- paste(out,X[i],sep='+')
    }
    return(out)
}

#' Group components of a composition
#'
#' Adds several components of a composition together into a single component
#' @param X a compositional dataset
#' @param ... a series of new labels assigned to strings or vectors of strings
#' denoting the components that need amalgamating
#' @return an object of the same class as X with fewer components
#' @examples
#' data(Namib)
#' HMcomponents <- c("zr","tm","rt","TiOx","sph","ap","ep",
#'                   "gt","st","amp","cpx","opx")
#' am <- amalgamate(Namib$PTHM,feldspars=c("KF","P"),
#'                  lithics=c("Lm","Lv","Ls"),heavies=HMcomponents)
#' plot(ternary(am))
#' @rdname amalgamate
#' @export
amalgamate <- function(X,...){ UseMethod("amalgamate",X) }
#' @rdname amalgamate
#' @export
amalgamate.default <- function(X,...){
    groups <- list(...)
    ng <- length(groups)
    labels <- names(groups)
    out <- NULL
    for (i in 1:ng){
        colsum <- sumcols(X,groups[[i]])
        out <- cbind(out,colsum)
    }
    colnames(out) <- labels
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.compositional <- function(X,...){
    out <- X
    out$x <- amalgamate(X$x,...)
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.SRDcorrected <- function(X,...){
    out <- X
    out$x <- amalgamate.default(X$x,...)
    for (sname in names(X$restoration)){
        out$restoration[[sname]] <-
            amalgamate.default(X$restoration[[sname]],...)
    }
    return(out)
}

#' Combine samples of distributional data
#'
#' Lumps all single grain analyses of several samples together under a new name
#' @param X a distributional dataset
#' @param ... a series of new labels assigned to strings or vectors of strings
#' denoting the samples that need amalgamating
#' @return a distributional data object with fewer samples than X
#' @examples
#' data(Namib)
#' combined <- combine(Namib$DZ,east=c('N3','N4','N5','N6','N7','N8','N9','N10'),
#'                        west=c('N1','N2','N11','N12','T8','T13'))
#' summaryplot(KDEs(combined))
#' @export
combine <- function(X,...){
    out <- X
    groups <- list(...)
    ng <- length(groups)
    labels <- names(groups)
    out$x <- list()
    out$err <- list()
    loadErr <- (length(X$err)>0)
    for (i in 1:ng){
        out$x[[labels[i]]] <- NULL
        for (g in groups[[i]]){
            out$x[[labels[i]]] <- c(out$x[[labels[i]]],X$x[[g]])
            if (loadErr) { out$err[[labels[i]]] <- c(out$err[[labels[i]]],X$err[[g]]) }
        }
    }
    return(out)    
}

ternaryclosure <- function(X,x,y,z){ 
    xlab <- sumlabels(x)
    ylab <- sumlabels(y)
    zlab <- sumlabels(z)
    out <- cbind(sumcols(X,x),sumcols(X,y),sumcols(X,z))
    den <- rowSums(out)
    out <- apply(out,2,'/',den)
    if (methods::is(out,"matrix")) {
        colnames(out) <- c(xlab,ylab,zlab)
    } else {
        names(out) <- c(xlab,ylab,zlab)
    }
    return(out)
}

#' Define a ternary composition
#'
#' Create an object of class \code{ternary}
#' @param X an object of class \code{compositional}
#' @param x string or a vector of strings indicating the variables making up
#' the first subcomposition of the ternary system. If omitted, the first
#' component of X is used instead.
#' @param y second (set of) variables
#' @param z third (set of) variables
#' @return an object of class \code{ternary}, i.e. a
#' list containing:
#'
#' x: a three column matrix (or vector) of ternary compositions.
#'
#' and (if X is of class \code{SRDcorrected})
#'
#' restoration: a list of intermediate ternary compositions inherited
#' from the SRD correction
#'
#' @seealso restore
#' @examples
#' data(Namib)
#' tern <- ternary(Namib$PT,c('Q'),c('KF','P'),c('Lm','Lv','Ls'))
#' plot(tern,type="QFL")
#' @export
ternary <- function(X,x=NULL,y=NULL,z=NULL){
    if (!methods::is(X,"compositional")) stop("X is not of class compositional")
    out <- list()
    class(out) <- append("ternary","compositional")
    if (is.null(x)) x <- colnames(X$x)[1]
    if (is.null(y)) y <- colnames(X$x)[2]
    if (is.null(z)) z <- colnames(X$x)[3]
    arg <- deparse(substitute(x))
    out$x <- ternaryclosure(X$x,x,y,z)
    if (methods::is(X,"SRDcorrected")){
        out$restoration <- list()
        snames <- names(X$restoration)
        for (sname in snames){
            out$restoration[[sname]] <-
                ternaryclosure(X$restoration[[sname]],x,y,z)
        }
        class(out) <- append("SRDcorrected",class(out))
    }
    return(out)
}
