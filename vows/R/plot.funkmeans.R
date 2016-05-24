#' Plotting of k-means clustering results for massively parallel smooths
#' 
#' Visualization of functional k-means clustering as implemented by
#' \code{\link{funkmeans}}.
#' 
#' 
#' @param x a functional k-means clustering object obtained from
#' \code{\link{funkmeans}}.
#' @param fdobj a functional data object, of class \code{"\link[fda]{fd}"},
#' defining the set of curves being clustered. By default, this is taken to be
#' \code{x$fdobj}; but if the latter is \code{NULL}, \code{fdobj} must be
#' specified. See the two cases in the example.
#' @param deriv which derivative to display in the plots, which show 30
#' randomly selected curves, along with the cluster center, from each cluster.
#' By default, the "0th derivative" is used (i.e., the curves themselves).
#' @param ncluster number of clusters to display.  By default, all are
#' displayed.
#' @param new.array logical: if \code{TRUE}, plots will be displayed in an 
#' array whose dimensions are set by the \code{mfrow} argument.
#' @param mfrow a vector of length 2 giving the numbers of rows and columns for
#' the array of plots. By default, the number of rows will exceed the number of
#' columns by \code{0} or \code{1}, depending on \code{ncluster}.
#' @param colvec a vector of colors for the clusters. By default, this is set
#' to the first \code{ncluster} elements of \code{c("dodgerblue", "green",
#' "red", "orange", "yellow", "orchid",} \code{ "brown", "grey", "purple")}, if
#' \code{ncluster <= 9}.
#' @param cex.mtext magnification for mtext command to display the size of each
#' cluster above the corresponding subfigure.
#' @param xlabs,ylabs,titles ?????\code{NULL} or a character vector of length 1 
#' or \code{ncluster}, specifying titles (x axis, y axis, overall titles) for each
#'  cluster. If vector's length equals 1, each cluster plot has the same title. 
#'  By default, it's \code{NULL}
#' @param ... arguments passed to \code{\link[graphics]{plot}.}
#' @author Yin-Hsiu Chen \email{enjoychen0701@@gmail.com}, Philip Reiss
#' \email{phil.reiss@@nyumc.org}, Lan Huo, and Ruixin Tan
#' @seealso \code{\link{funkmeans}}, \code{\link{funkpanel}}
#' @examples
#' 
#' data(test)
#' d4 = test$d4
#' x = test$x
#' semi.obj = semipar4d(d4, formula = ~sf(x), data = data.frame(x = x), lsp=-5:5)
#' myfdobj = extract.fd(semi.obj)
#' 
#' # Case 1: fd object is stored in funkmeans object...
#' fkmobj = funkmeans(myfdobj, ncomp = 8, centers = 6)
#' plot(fkmobj)
#' 
#' # Case 2: fd object is not stored...
#' fkmobj = funkmeans(myfdobj, ncomp = 8, centers = 6, store.fdobj=FALSE)
#' plot(fkmobj, myfdobj)
#' @export
plot.funkmeans = function(x, fdobj=NULL, deriv=0, ncluster = nrow(x$centers), new.array=TRUE, mfrow = NULL,
            colvec = NULL, cex.mtext=.7, xlabs="", ylabs="", titles = "", ...) {
    if (is.null(fdobj)) {
    	if (is.null(x$fdobj)) stop("Must specify 'fdobj'")
    	fdobj = x$fdobj
    }
    if (new.array) {
        if (!is.null(mfrow)) par(mfrow=mfrow)
        else {
    	      nro = ceiling(sqrt(ncluster))
    	      nco = ceiling(ncluster/nro)
    	      par(mfrow=c(nro,nco))
        }
    }
    if (is.null(colvec)) {
    	if (ncluster<=9) colvec = c("dodgerblue", "green", "red", "orange", "yellow", "orchid", "brown", "grey", "purple")[1:ncluster]
        else stop("Please set 'colvec' to a numeric or character vector of colors")
    }
    for (ii in 1:ncluster) {
        ee = sample(which(x$cluster==ii), min(30, sum(x$cluster==ii)))
        ee.fd = deriv.fd(fd(coef = fdobj$coef[ , ee], basisobj = fdobj$basis), deriv)
        
        if (length(xlabs) <= 1) xlab = xlabs
          else if (length(xlabs) == ncluster) xlab = xlabs[ii]
          else stop("The length of xlabs should be 1 or equal to the number of clusters")
        if (length(ylabs) <= 1) ylab = ylabs
          else if (length(ylabs) == ncluster) ylab = ylabs[ii]
          else stop("The length of ylabs should be 1 or equal to the number of clusters")
        if (length(titles) <= 1) main = titles
          else if (length(titles) == ncluster) main = titles[ii]
          else stop("The length of titles should be 1 or equal to the number of clusters")
        
        plot(ee.fd, lty=1, col=colvec[x$cluster[ee]], xlab=xlab, ylab=ylab, main = main, ...)   
        mtext(sum(x$cluster==ii), cex=cex.mtext)        
        idx = which(x$cluster==ii)
        mean.cluster = deriv.fd(fd(coef = apply(fdobj$coef[ , idx],1,mean), basisobj = fdobj$basis), deriv)
        lines(mean.cluster, lty=1, col="black", lwd=2)
    }
}
