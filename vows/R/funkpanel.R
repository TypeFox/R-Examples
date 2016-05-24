#' Interactive visualization of massively parallel smoothing results
#' 
#' This function uses interactive graphics tools, provided by the \pkg{rpanel}
#' package (Bowman et al., 2007), to visualize the results of functional
#' k-means clustering as implemented by \code{\link{funkmeans}}.
#' 
#' The clusters defined by the input object are shown for a cross-section
#' (slice), and a slider allows viewing of different slices.  Clicking on a
#' voxel produces a scatterplot of the data at that voxel, along with the
#' fitted scatterplot.  The "legend", produced by \code{\link{plot.funkmeans}},
#' consists of plots of 30 randomly selected curves, along with the cluster
#' center, from each cluster.
#' 
#' @param fkmobj a functional k-means clustering object obtained from
#' \code{\link{funkmeans4d}}.
#' @param semiobj the massively parallel smoothing object on which the
#' clustering was based; ordinarily produced by \code{\link{semipar.mp}} or
#' \code{\link{semipar4d}}.
#' @param arr4d a 4-dimensional array containing the raw data that were
#' smoothed at each point.  The first 3 dimensions refer to x, y, and z
#' coordinates and the last dimension corresponds to different images.
#' @param predictor a vector or matrix of covariates.
#' @param titl title of the panel.
#' @param xlab,ylab x- and y-axis labels.
#' @param ncluster number of clusters to display.  By default, all are
#' displayed.
#' @param slice index of the slice to be shown initially in the panel.
#' @param ylim.scatter the y limits of the voxelwise scatterplots.
#' @param deriv.legend which derivative to plot in the "legend"; see Details.
#' By default, the curves themselves are used.
#' @param ylim.legend the y limits used in the "legend"; see Details.
#' @param scattermain title for the scatter plots.
#' @param colvec a vector of colors for the clusters. By default, this is set
#' to the first \code{ncluster} elements of \code{c("dodgerblue", "green",
#' "red", "orange", "yellow", "orchid",} \code{ "brown", "grey", "purple")}, if
#' \code{ncluster < 9}.
#' @author Lei Huang \email{huangracer@@gmail.com}, Yin-Hsiu Chen
#' \email{enjoychen0701@@gmail.com}, and Lan Huo
#' @seealso \code{\link{funkmeans}}, \code{\link{funkmeans4d}},
#' \code{\link{plot.funkmeans}}
#' @references Bowman, A., Crawford, E., Alexander, G., and Bowman, R. (2007).
#' rpanel: Simple interactive controls for R functions using the tcltk package.
#' \emph{Journal of Statistical Software}, 17(9).
#' @examples
#' 
#' data(test)
#' d4 = test$d4
#' x = test$x
#' semi.obj = semipar4d(d4, ~sf(x), -5:5, data.frame(x = x))
#' fdobj = extract.fd(semi.obj)
#' fkmobj = funkmeans4d(fdobj, d4, ncomp=6, centers=3)
#' funkpanel(fkmobj, semi.obj, d4, x)
#' @export
funkpanel <- function(fkmobj, semiobj, arr4d, predictor, titl="", xlab="", ylab="", ncluster = nrow(fkmobj$centers), 
                      slice=dim(fkmobj$arr.cluster)[3]%/%2, ylim.scatter=NULL, deriv.legend=0, ylim.legend=NULL, 
                      scattermain=NULL, colvec = NULL)  {
    if (!inherits(fkmobj, "funkmeans4d")) stop("'fkmobj' must be an object of class 'funkmeans4d'")
    if (is.null(colvec)) {
    	if (ncluster<=9) colvec = c("dodgerblue", "green", "red", "orange", "yellow", "orchid", "brown", "grey", "purple")[1:ncluster]
        else stop("Please set 'colvec' to a numeric or character vector of colors")
    }
    tkrp = tkrp1 = NULL
    screened = min(fkmobj$arr.cluster, na.rm=TRUE) == 0
    if (!is.null(semiobj$incl.inds)) attr(arr4d,"has.data")= semiobj$include
    data.inds = which(attributes(arr4d)$has.data==TRUE, arr.ind=TRUE)
    x.coord = attributes(arr4d)$coord[[1]]
    y.coord = attributes(arr4d)$coord[[2]]
    z.coord = attributes(arr4d)$coord[[3]]
    axis.flag = TRUE
    ttl = "z"
   
    draw <- function(panel) {
        image(x = x.coord, y = y.coord, z = fkmobj$arr.cluster[ , , panel$slice], 
              col = if (screened) c("grey", colvec) else colvec, breaks = .5 + (-screened):ncluster, 
              main = paste(ttl, "=", z.coord[panel$slice]), xlab="x", ylab="y", axes = axis.flag)
        panel
    }
    
    redraw <- function(panel) {
        rpanel::rp.tkrreplot(panel, tkrp)
        panel
    }
    
    scatterplot <- function(panel,x,y)  {
        arr.ind = c(which.min(abs(x.coord-x)), which.min(abs(y.coord-y)), panel$slice)
        incl.inds = if (is.null(semiobj$incl.inds)) data.inds else semiobj$incl.inds
        which.vox = which(incl.inds[,1]==arr.ind[1] & incl.inds[,2]==arr.ind[2] & incl.inds[,3]==arr.ind[3])
        if (length(which.vox)>0) plot(semiobj, semiobj$Y, which.vox=which.vox, ylab=ylab, main=scattermain, ylim=ylim.scatter)
        else cat("Scatterplots available only for voxels within the displayed clusters!\n")
        return(panel)
   }
 
    legend.draw <- function(panel){
    	fdobj = list()
    	fdobj$coef <- fkmobj$coef
    	fdobj$basis <- fkmobj$basis
    	plot(fkmobj, fdobj = fdobj, xlab = xlab, ylab = ylab, deriv = deriv.legend, 
    	     ylim = ylim.legend, ncluster = ncluster, colvec = colvec)
    	panel
    }
    
    imgplot <- rpanel::rp.control(title = titl, slice=slice)
    rpanel::rp.tkrplot(imgplot, tkrp, draw, action=scatterplot, hscale = 0.8, vscale = 1.2, pos=list(row=0, column=0))
    rpanel::rp.tkrplot(imgplot, tkrp1, legend.draw, hscale = 0.8, vscale = 1.2, pos=list(row=0, column=1))
    rpanel::rp.slider(panel = imgplot, variable=slice, action = redraw, from = 1, to = dim(fkmobj$arr.cluster)[3], resolution = 1, title=ttl, showvalue=is.null(attributes(arr4d)$coord), pos=list(row=1, column=1))
}

