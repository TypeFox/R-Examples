#
#  EM algorithm for Mixture of Unrestricted Multivariate Skew t-distributioins
#  Package: EMMIX-uskew
#  Version: 0.11-5
#
#  Code by S. X. Lee
#  Updated on 01 Nov 2013
#
# Lee S. and Mclachlan, G.J. (2013) On the fitting of finite mixtures of
#   multivariate skew t-distributions via the EM algorithm
#

################################################################################
#  SECTION 4
#                           Contour Plots
#
################################################################################
#2D Contours
fmmst.contour.2d <- function(dat, model, grid=50, drawpoints=TRUE, clusters=NULL,
    levels=10, map=c("scatter", "heat", "cluster"), 
    component=NULL, xlim, ylim, xlab, ylab, main, tmethod=1, ...) {  
    if(length(levels)==1) {
        if(levels<0) stop("'levels' cannot have negative value")
        if(is.whole(levels)) lflag<- T
        else if(levels>1) stop("'levels' must be either an integer or a quantile value") else lflag<-F
    }else if (all(levels > 0) && all(levels<1)) lflag <- F
    else stop("'levels' must be either an integer or a vector of quantile values")   
    if(ncol(dat)<2) stop("Data must be 2D!")
    if(ncol(dat)!=2) warning("Data should be 2-D. Only the first two columns will be used.")   
    if(is.null(clusters)) clusters <- rep(1, nrow(dat))
    dat<-as.matrix(dat)[,1:2]
    rx<-range(dat[,1])+c(-1,1)
    ry<-range(dat[,2])+c(-1,1)
    if (missing(xlim))  xlim<-rx
    if (missing(ylim))  ylim<-ry 
    g<-length(model$dof) ;p<-2
    if (missing(xlab))  xlab<-"X"
    if (missing(ylab))  ylab<-"Y" 
    if (missing(main)) {
        msg <- bquote(paste("Multivariate Skew ", italic(t), "-mixture Model (g = ", .(g), ")", sep=""))
        #msg <- paste("Multivariate Skew t-mixture Model (g = ",g,")", sep="")
    } else msg=main
    x <- seq(xlim[1], xlim[2], length=grid)
    y <- seq(ylim[1], ylim[2], length=grid)
    nx <- length(x)
    ny <- length(y)
    xoy <- cbind(rep(x,ny), as.vector(matrix(y,nx,ny,byrow=TRUE)))
    X <- matrix(xoy, nx*ny, 2, byrow=FALSE)
    
    ColClust <- c("lightyellow", "lightgreen", "lightblue3", "lightpink", "lightgrey", 
        "lightgreen", "lightsalmon", colors()[sample(1:length(colors()),g)])
    ColPoint <- c("blue", "red", "black", "Brown", "green", colors()[sample(1:length(colors()),g)])         
    
    if(map[1]=="heat") {    
        smoothScatter(dat, nrpoints=Inf, colramp=colorRampPalette(heatColor()),xlab=xlab,ylab=ylab,main=msg, ...)
        drawpoints <- FALSE
    }else if(map[1]=="cluster") {
        cat("creating cluster map\n")
        clustergrid <- fmmstDA(g, X, model, tmethod)
        plot(x = 0, y = 0,type = "n", xlim = xlim, ylim = ylim,xlab=xlab,ylab=ylab, main=msg)        
        for(i in 1:g) points(X[clustergrid==i,], pch=22, lwd=grid/5, col=ColClust[i], bg=ColClust[i])
    }else plot(x = 0, y = 0,type = "n", xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab, main=msg)
    if(drawpoints) {for(i in 1:g) points(dat[clusters==i,],pch=19,col=ColPoint[i])}

    if(is.null(component)) {
        dens <- dfmmst(X, model$mu, model$sigma, model$delta, model$dof, model$pro, tmethod=tmethod)
        densMat <- matrix(dens,nx,ny)
        if(lflag) contour(x, y, densMat, nlevels=levels, add=TRUE, drawlabels=FALSE, lty=1, col = 4,...)
        else {
            levels <- quantile(densMat, probs=levels)
            contour(x, y, densMat, levels=levels, add=TRUE, drawlabels=FALSE, lty=1, col = 4,...)         
        }
    } else {
        ncomp <- length(as.numeric(component))
        if(!all(component %in% c(1:g))) stop("argument 'component' is incorrectly specified!")
        if(!lflag) levels <- quantile(densMat, probs=levels)
        for(i in 1:ncomp) {
            dens <- dmst(X, model$mu[[component[i]]], 
                model$sigma[[component[i]]], model$delta[[component[i]]], 
                model$dof[component[i]], tmethod=tmethod) 
            densMat <- matrix(dens,nx,ny)     
            if(lflag) contour(x, y, densMat, nlevels=levels, add=TRUE, 
                drawlabels=FALSE, lty=1,col = ColPoint[i],...) 
            else contour(x, y, densMat, levels=levels, add=TRUE, 
                    drawlabels=FALSE, lty=1,col = ColPoint[i],...)    
        }
    }  
} 

#perform discriminant analysis (DA) for a specified FM-MST distribution
fmmstDA <- function(g, dat, model, tmethod=1) {
    if (g == 1) return(rep(1, nrow(dat)))
    probMat <- matrix(0,g,nrow(dat))
    for(i in 1:g) probMat[i,] <- model$pro[i]*dmst(dat, model$mu[[i]], model$sigma[[i]], model$delta[[i]], model$dof[i], tmethod=tmethod)
    clust <- apply(probMat,2,which.max)  
    return(clust)
}

#3D Contours
fmmst.contour.3d <- function(dat, model, grid=20, drawpoints=TRUE, levels=0.9, 
    clusters=NULL, xlim, ylim, zlim, xlab, ylab, zlab, main, component=NULL, ...) {
         
    if(length(levels)==1) {
        if(levels<0) stop("'levels' cannot have negative value")
        if(is.whole(levels)) lflag<- T
        else {if(levels>1) stop("'levels' must be either an integer or a quantile value") else lflag<-F}
    }else if (all(levels > 0) && all(levels<1)) lflag <- F
    else stop("'levels' must be either an integer or a vector of quantile values") 
    if(missing(dat)) {
        if(drawpoints)  stop("dat must be provided if drawpoints=TRUE.")
        if(missing(xlim) || missing(ylim) || missing(zlim)) stop("x, y, z limits must be provided if dat is not provided.")
    } else {
        if(ncol(dat)<3) stop("Data must be 3D!")
        if(ncol(dat)!=3) warning("Data should be 3-D. Only the first three columns will be used.")  
        if(is.null(clusters)) clusters <- rep(1, nrow(dat))
        dat<-as.matrix(dat)[,1:3]   
        rx<-range(dat[,1])+c(-1,1)
        ry<-range(dat[,2])+c(-1,1) 
        rz<-range(dat[,3])+c(-1,1)  
        if (missing(xlim))  xlim<-rx
        if (missing(ylim))  ylim<-ry  
        if (missing(zlim))  zlim<-rz 
    }
    g<-length(model$dof) ;p<-3
    if (missing(xlab))  xlab<-"X"
    if (missing(ylab))  ylab<-"Y" 
    if (missing(zlab))  zlab<-"Z" 
    library(rgl)
    ColPoint <- c("magenta", "green", "red", "blue", "grey", colors()[sample(1:length(colors()),g)]) 
       
    if(drawpoints) {
        plot3d(dat[clusters==1,1], dat[clusters==1,2], dat[clusters==1,3], 
            col="magenta", xlab=xlab, ylab=ylab, zlab=zlab, xlim=xlim, ylim=ylim, zlim=zlim)
        if(g>=2){ 
            for(i in 2:g)            
            points3d(dat[clusters==i,1], dat[clusters==i,2], dat[clusters==i,3], col=ColPoint[i])        
        }
    } else {
        plot3d(1,1,1,col="white", xlab=xlab, ylab=ylab, zlab=zlab, xlim=xlim, ylim=ylim, zlim=zlim)        
    } 
    XVec <- seq(xlim[1], xlim[2], length=grid)       
    YVec <- seq(ylim[1], ylim[2], length=grid)
    ZVec <- seq(zlim[1], zlim[2], length=grid) 
    if(is.null(component)) {
        mst.contour <- function(x, y, z) dfmmst(cbind(x,y,z), model$mu, 
              model$sigma, model$delta, model$dof, model$pro)
        contour3d(mst.contour, levels, XVec, YVec, ZVec, add=T, color="blue", fill=T, alpha=0.25, nflag=lflag)
    } else{
        ncomp <- length(as.numeric(component))
        if(!all(component %in% c(1:g))) stop("argument 'component' is incorrectly specified!")
        for(i in 1:ncomp){
          mst.contour <- function(x, y, z) dmst(cbind(x,y,z), model$mu[[component[i]]], 
              model$sigma[[component[i]]], model$delta[[component[i]]], model$dof[component[i]])
          cat("Plotting component ", component[i] ,"\n")
          contour3d(mst.contour, levels, XVec, YVec, ZVec, add=T, color=ColPoint[component[i]], fill=T, alpha=0.25, nflag=lflag)
          }
    }  
}


################################################################################
#other functions
heatColor <- function(){
    rgb(c(255,255,254,254,253,252,227,189,128),
        c(255,237,217,178,141,78,26,0,0),
        c(204,160,118,76,60,42,28,38,38),maxColorValue=255)                     
}

#The following code is adapted from the package misc3d (with minor changes)
contour3d <- function(f, level, x = 1:dim(f)[1], y = 1:dim(f)[2], z = 1:dim(f)[3],
    mask = NULL, color = "white", color2 = NA, alpha = 1, fill = TRUE, 
    col.mesh = if (fill) NA else color, material = "default", smooth = 4,
    add = FALSE, draw = TRUE, engine = "rgl", separate=FALSE, nflag=F, ...){

    if (! all(is.finite(x), is.finite(y), is.finite(z)))
        stop("'x', 'y', and 'z' values must be finite and non-missing")
    if (is.function(f) || is.array(f) && length(dim(f)) == 3){
        if (is.function(f)){
            if (length(formals(f)) < 3)
                stop("The function must have at least 3 arguments.")
            vol <- fgrid(f, x, y, z)
        }else{
          if (dim(f)[1] != length(x) || dim(f)[2] != length(y) ||  dim(f)[3] != length(z))
            stop("dimensions of f do not match x, y, or z")
          vol <- f
        }              
        maxvol <- max(vol)
        minvol <- min(vol)
        if(missing(level)) {quan<-quantile(vol); level <- quan[4];}
        else {
            if(nflag) level <- quantile(vol, probs=seq(0,1,length.out=level))
            else  level <- quantile(vol, probs=level)
        }
        #cat("lflag is ", nflag, "\n")
        #cat("Your quantile levels are ", level, "\n")
        con <- which(! level <= maxvol & level >= minvol)
        if (length(con) == length(level))
            stop(paste("The 'level' has to be within the range of 'f' (between ", minvol, " and ", maxvol,").\n", sep=""))
        else if (length(con) > 0){
            warning(paste("The 'level' outside the range of 'f' (between ", minvol, " and ", maxvol, ") has been removed. \n", sep=""))
            level <- level[-con]
            if (is.list(mask)) mask <- mask[-con]
            if (length(color) > 1) color <- color[-con]
            if (length(color2) > 1) color2 <- color2[-con]
            if (length(alpha) > 1) alpha <- alpha[-con]
            if (length(fill) > 1) fill <- fill[-con]
            if (length(col.mesh) > 1) col.mesh <- col.mesh[-con]
            if (length(material) > 1) material <- material[-con]
            if (length(smooth) > 1) smooth <- smooth[-con]
        }

      }

    else stop("vol has to be a function or a 3-dimensional array")

    scene <- contourTriangles(vol, maxvol, level, x, y, z, mask, color, color2,
                              alpha, fill, col.mesh, material, smooth)
    if (! draw || engine == "none"){
        if (! any(separate))
            scene
        else{
            if (length(level)==1){
                newScene <- separateTriangles(scene)
                cat("Triangles are separated into ", length(newScene),
                          " chunks.", "\n", sep="")
            }
            else{
                if (length(separate) < length(level))
                    separate <- c(separate, rep(FALSE, length(level)-length(separate)))

                newScene <- NULL
                for (i in 1:length(level)){
                    if (separate[i]){
                        new <- separateTriangles(scene[[i]])
                        newScene <- c(newScene, new)
                        cat("Triangles from level ", level[i],
                                  " are separated into ", length(new), " chunks.",
                                  "\n", sep="")
                }
                    else
                        newScene <- c(newScene, list(scene[[i]]))
                }
            }
            newScene
        }
    }
    else {
        scene <- colorScene(scene)
        if (engine == "rgl")
            drawScene.rgl(scene, add = add, ...)
        else if (engine %in% c("standard", "grid"))
            drawScene(scene, add = add, engine = engine, ...)
        else stop(paste("unknown rendering engine:", engine))
    }
}

PreProcessing <- local({
    explode <- function(x)
        floor(((x - 1) %% 2^(1:8))/2^(0:7))

    BasicRotation <-
        matrix(c(1,2,3,4,5,6,7,8,5,6,2,1,8,7,3,4,8,7,6,5,4,3,2,1,
                 4,3,7,8,1,2,6,5,2,6,7,3,1,5,8,4,6,5,8,7,2,1,4,3,
                 5,1,4,8,6,2,3,7,4,1,2,3,8,5,6,7,3,4,1,2,7,8,5,6,
                 2,3,4,1,6,7,8,5,6,7,3,2,5,8,4,1,7,8,4,3,6,5,1,2,
                 8,5,1,4,7,6,2,3,7,3,2,6,8,4,1,5,4,8,5,1,3,7,6,2,
                 3,2,6,7,4,1,5,8,2,1,5,6,3,4,8,7,1,4,8,5,2,3,7,6,
                 1,5,6,2,4,8,7,3,5,8,7,6,1,4,3,2,8,4,3,7,5,1,2,6,
                 3,7,8,4,2,6,5,1,7,6,5,8,3,2,1,4,6,2,1,5,7,3,4,8),
               ncol=8, byrow=TRUE)

    CaseRotation <-
        matrix(c(1,24,2,19,2,17,3,17,2,24,4,24,3,24,6,10,2,15,3,19,
                 4,17,6,9,3,9,6,8,6,1,9,23,2,20,3,18,4,7,6,16,5,24,
                 7,5,7,24,12,9,4,20,6,22,8,24,10,24,7,9,15,24,13,20,
                 6,20,2,21,4,6,3,16,6,4,4,16,8,23,6,14,10,23,5,21,
                 7,10,7,16,15,9,7,2,13,8,12,23,6,6,3,6,6,17,6,18,
                 9,18,7,4,13,17,15,18,6,13,7,6,12,16,13,18,6,2,11,24,
                 7,3,7,12,3,12,2,23,5,23,4,23,7,1,3,14,7,14,6,21,
                 15,23,4,15,7,19,8,19,13,23,6,11,12,17,10,19,6,23,4,12,
                 7,18,8,22,13,16,7,13,11,23,13,21,7,15,8,21,13,22,14,24,
                 8,15,13,11,7,7,8,12,4,22,3,23,7,23,6,24,12,18,6,7,
                 13,19,9,24,6,19,7,21,11,18,13,24,7,20,15,16,7,22,6,15,
                 3,22,6,3,15,17,10,22,6,12,12,24,7,11,6,5,3,15,13,10,
                 7,8,8,20,4,9,7,17,5,22,4,18,2,22,2,22,4,18,5,22,
                 7,17,4,9,8,20,7,8,13,10,3,15,6,5,7,11,12,24,6,12,
                 10,22,15,17,6,3,3,22,6,15,7,22,15,16,7,20,13,24,11,18,
                 7,21,6,19,9,24,13,19,6,7,12,18,6,24,7,23,3,23,4,22,
                 8,12,7,7,13,11,8,15,14,24,13,22,8,21,7,15,13,21,11,23,
                 7,13,13,16,8,22,7,18,4,12,6,23,10,19,12,17,6,11,13,23,
                 8,19,7,19,4,15,15,23,6,21,7,14,3,14,7,1,4,23,5,23,
                 2,23,3,12,7,12,7,3,11,24,6,2,13,18,12,16,7,6,6,13,
                 15,18,13,17,7,4,9,18,6,18,6,17,3,6,6,6,12,23,13,8,
                 7,2,15,9,7,16,7,10,5,21,10,23,6,14,8,23,4,16,6,4,
                 3,16,4,6,2,21,6,20,13,20,15,24,7,9,10,24,8,24,6,22,
                 4,20,12,9,7,24,7,5,5,24,6,16,4,7,3,18,2,20,9,23,
                 6,1,6,8,3,9,6,9,4,17,3,19,2,15,6,10,3,24,4,24,
                 2,24,3,17,2,17,2,19,1,24),
               ncol=2,byrow=TRUE)

    CaseRotationFlip <-
        cbind(CaseRotation,
              c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,
                1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,-1,-1,1,1,
                1,1,1,1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,1,
                -1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,1,1,1,1,
                1,-1,1,-1,1,1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,-1,-1,-1,-1,-1,
                -1,-1,
                -1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,1,1,1,1,-1,
                -1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,-1,
                -1,-1,
                -1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,
                -1,-1,
                1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,
                -1,-1,1,
                1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,
                -1,-1,-1,
                -1,-1,-1,-1,-1,-1))
    EdgePoints <-
        matrix(c(1,1,2,2,2,3,3,3,4,4,4,1,5,5,6,6,6,7,
                 7,7,8,8,8,5,9,1,5,10,2,6,11,3,7,12,4,8,13,9,9),
               ncol=3, byrow=TRUE)
    BasicEdges <- list(c(1,4,9),
                       c(2,4,9,10),
                       c(1,4,5,6,9,10),
                       c(1,4,6,7,9,11),
                       c(1,4,10,11,12),
                       c(2,4,6,7,9,10,11),
                       c(1,2,5,6,7,8,9,10,11,13),
                       c(9,10,11,12),
                       c(1,2,7,8,9,11),
                       c(1,2,3,4,5,6,7,8,13),
                       c(1,2,6,7,9,12),
                       c(1,4,5,8,9,10,11,12,13),
                       c(1,2,3,4,5,6,7,8,9,10,11,12,13),
                       c(1,4,7,8,10,11))
    EdgeSequence1 <-
        list(c(1,2,3),
             c(1,2,4,2,3,4),
             list(c(1,2,5,3,4,6),
                  c(1,2,6,2,4,6,2,3,4,2,5,3)),
             list(c(1,2,5,3,4,6),
                  c(1,3,5,3,4,5,2,5,4,1,2,6,1,6,3,2,4,6)),
             c(1,3,2,2,3,5,3,4,5),
             list(c(1,2,6,2,5,6,3,4,7),
                  c(1,2,7,2,4,7,2,5,4,3,4,5,3,5,6,1,6,3,1,3,7,7,3,4),
                  c(1,2,7,2,4,7,2,5,4,3,4,5,3,5,6)),
             list(c(1,8,2,3,7,6,4,5,9),
                  c(1,8,2,3,9,4,3,7,9,5,9,7,5,7,6),
                  c(3,7,6,1,8,4,1,4,5,1,9,2,1,5,9),
                  c(4,5,9,1,7,6,1,6,2,2,6,3,2,3,8),
                  c(1,10,2,2,10,9,5,9,10,5,10,6,6,10,7,3,7,10,3,10,4,4,10,8,1,
                    8,10),
                  c(1,10,2,1,10,7,6,10,7,5,10,6,5,9,10,4,9,10,3,10,4,3,8,10,2,
                    10,8),
                  c(5,9,10,2,10,9,1,10,2,1,7,10,6,10,7,3,10,6,3,8,10,4,10,8,4,
                    5,10),
                  c(1,7,6,1,6,9,1,9,2,5,9,6,3,8,4),
                  c(1,7,8,3,8,7,3,7,6,3,6,5,3,5,4,4,5,9,4,9,8,2,8,9,1,8,2)),
             c(1,2,3,1,3,4),
             c(1,2,6,1,3,5,1,6,3,3,4,5),
             list(c(1,4,5,2,6,3,3,6,7,4,8,5),
                  c(1,4,3,1,3,2,3,4,8,3,8,7,5,8,6,6,8,7,1,2,6,1,6,5),
                  c(1,2,9,1,9,5,5,9,8,4,8,9,3,4,9,3,9,7,6,7,9,2,6,9),
                  c(5,9,6,1,9,5,1,4,9,4,8,9,7,9,8,3,9,7,2,9,3,2,6,9),
                  c(1,2,5,2,6,5,3,4,8,3,8,7)),
             c(1,2,3,1,3,6,1,6,5,3,4,6),
             list(c(1,6,2,2,6,8,3,5,4,6,7,8),
                  c(1,5,3,1,3,6,3,6,7,3,7,4,2,4,8,2,5,4,1,5,2,4,7,8),
                  c(1,9,2,2,9,5,3,5,9,4,9,8,3,9,4,7,8,9,6,7,9,1,9,6),
                  c(4,9,5,1,5,9,1,9,2,2,9,8,7,8,9,6,9,7,3,6,9,3,9,4),
                  c(1,5,2,3,8,4,3,6,8,6,7,8)),
             list(
                  ##13.1
                  c(7,8,12,2,3,11,1,4,9,5,6,10),
                  ##13.2
                  c(2,  3, 11,  7,  8, 12,  9,  5,  4,  5,  4,  6,  4,  6,  1,
                    6,  1, 10),
                  c(1,  4,  9,  7,  8, 12, 10,  2,  5,  2,  5,  3,  5,  3,  6,
                    3,  6, 11),
                  c(5,  6, 10,  1,  4,  9, 11,  7,  2,  7,  2,  8,  2,  8,  3,
                    8,  3, 12),
                  c(5,  6, 10,  2,  3, 11, 12,  4,  7,  4,  7,  1,  7,  1,  8,
                    1,  8,  9),
                  c(5,  6, 10,  7,  8, 12,  2, 11,  1, 11,  1,  9, 11,  9,  3,
                    9,  3,  4),
                  c(2,  3, 11,  4,  1,  9,  5, 10,  8, 10,  8, 12, 10, 12,  6,
                    12,  6,  7),
                  ##13.3
                  c(7, 8, 12, 13, 3, 11, 13, 11, 6, 13,  6, 5, 13, 5,  9, 13,
                    9,  4, 13,  4,  1, 13, 1,10,13,10, 2,13,2,3),
                  c(2, 3, 11, 13, 6, 10, 13, 10, 1, 13,  1, 4, 13, 4, 12, 13,
                    12,  7, 13,  7,  8, 13, 8, 9,13,9, 5, 13,5,6),
                  c(7, 8, 12, 13, 6,  5, 13,  5, 9, 13,  9, 4, 13, 4,  3, 13,
                    3, 11, 13, 11,  2, 13, 2, 1,13,1,10,13,10,6),
                  c(2, 3, 11, 13,  4,  1, 13,  1, 10, 13, 10, 6, 13, 6, 7,13,
                    7, 12, 13, 12,  8, 13, 8, 5,13,5, 9,13, 9,4),
                  c(1, 4, 9,  13,  8, 12, 13, 12,  3, 13,  3, 2, 13, 2, 10, 13,
                    10, 5, 13,  5, 6, 13, 6,11,13,11,7,13, 7,8),
                  c(7, 8, 12, 13,  5,  6, 13,  6, 11, 13, 11, 3, 13, 3,  4, 13,
                    4, 9, 13,  9, 1, 13, 1, 2,13,2,10,13,10,5),
                  c(1, 4,  9, 13,  3,  2, 13,  2, 10, 13, 10, 5, 13, 5,  8, 13,
                    8, 12,13, 12, 7, 13, 7, 6,13,6,11,13,11,3),
                  c(5, 6, 10, 13,  1,  9, 13,  9,  8, 13,  8, 7, 13, 7, 11, 13,
                    11,  2,13,  2,  3, 13, 3, 12, 13, 12, 4,13,  4, 1),
                  c(5, 6, 10, 13,  8,  7, 13,  7, 11, 13, 11, 2, 13, 2,  1, 13,
                    1,  9,13,  9,  4, 13, 4,  3, 13, 3, 12,13, 12, 8),
                  c(1, 4,  9, 13,  2,  3, 13,  3, 12, 13, 12, 8, 13, 8,  5, 13,
                    5, 10, 13, 10, 6, 13, 6,  7, 13, 7, 11,13, 11, 2),
                  c(5, 6, 10, 13,  7,  8, 13,  8,  9, 13,  9, 1, 13, 1,  2, 13,
                    2, 11, 13, 11, 3, 13, 3,  4, 13, 4, 12,13, 12, 7),
                  c(2, 3, 11, 13,  1,  4, 13,  4, 12, 13, 12, 7, 13, 7,  6, 13,
                    6, 10, 13, 10, 5, 13, 5,  8, 13, 8, 9, 13,  9, 1),
                  ##13.4
                  c(13, 3, 11, 13, 11, 6, 13, 6, 7, 13, 7, 12, 13, 12, 8, 13,
                    8, 5, 13, 5, 9, 13, 9, 4,
                    13, 4, 1,  13, 1, 10, 13,10, 2, 13, 2,  3),
                  c(13, 4, 12, 13, 12, 7, 13, 7, 8, 13, 8,  9, 13,  9, 5, 13,
                    5, 6, 13, 6, 10, 13, 10, 1, 13,
                    1,  2, 13,  2, 11, 13, 11,  3, 13,  3,  4),
                  c(13, 2, 10, 13, 10, 5, 13, 5, 6, 13, 6, 11, 13, 11, 7, 13,
                    7, 8, 13, 8, 12, 13,12, 3, 13,
                    3,  4, 13,  4,  9, 13,  9,  1, 13,  1,  2),
                  c(13, 1,  9, 13,  9, 8, 13, 8, 5, 13, 5, 10, 13, 10, 6, 13,
                    6, 7, 13, 7, 11, 13,11, 2, 13,
                    2,  3, 13,  3, 12, 13, 12,  4, 13,  4,  1),
                  ##13.5.1
                  c(7,  8, 12,  2,  1, 10,  3,  4, 11,  4, 11,  6,  4,  6,  9,
                    6,  9,  5),
                  c(3,  2, 11,  8,  5,  9,  4,  1, 12,  1, 12,  7,  1,  7, 10,
                    7, 10,  6),
                  c(1,  4,  9,  6,  7, 11,  2,  3, 10,  3, 10,  5,  3,  5, 12,
                    5, 12,  8),
                  c(5,  6, 10,  4,  3, 12,  1,  2,  9,  2,  9,  8,  2,  8, 11,
                    8, 11,  7),
                  ##13.5.2
                  c(1, 2, 10, 8, 5, 9,  8, 9,  4, 8, 4,12, 4,12, 3,12, 3,11,
                    12,11, 7,11, 7, 6, 7, 6, 8, 6, 8, 5),
                  c(8, 5,  9, 3, 4, 12, 3, 12, 7, 3, 7,11, 7,11, 6,11, 6,10,
                    11,10, 2,10, 2, 1, 2, 1, 3, 1, 3, 4),
                  c(6, 7, 11, 1, 2, 10, 1, 10, 5, 1, 5, 9, 5, 9, 8, 9, 8,12,
                    9,12, 4,12, 4, 3, 4, 3, 1, 3, 1, 2),
                  c(3, 4, 12, 6, 7, 11, 6, 11, 2, 6, 2, 10,2, 10, 1,10,1, 9,
                    10,9, 5, 5, 9, 8, 5, 8, 6, 8, 6,  7)),
             c(1,4,2,1,6,4,1,5,6,3,6,4))

    switch23 <- function(x){
        num <- length(x) / 3
        temp <- x[0: (num-1)*3 + 2]
        x[0: (num-1)*3 + 2] <- x[0: (num-1)*3+ 3]
        x[0: (num-1)*3 + 3] <- temp
        x
    }

    SwitchSeq <- function(ed){
        if (is.list(ed)){
            lapply(1:length(ed), function(x) SwitchSeq(ed[[x]]))
        }
        else
            switch23(ed)
    }

    EdgeSequence2 <- SwitchSeq(EdgeSequence1)

    getedge <- function(x){
        case <- x[1]
        rotation <- x[2]
        map <- rep(0,8)
        for(i in 1:8){
            temp <- as.integer(BasicRotation[rotation,][i])
            map[temp] <- i
        }
        sapply(BasicEdges[[case-1]], function(x){
            if (x!=13){
                EndP1 <- EdgePoints[x,2]
                EndP2 <- EdgePoints[x,3]
                newEdge <- EdgePoints[(EdgePoints[,2]==map[EndP1]
                                       &EdgePoints[,3]==map[EndP2])|
                                      (EdgePoints[,3]==map[EndP1]
                                       &EdgePoints[,2]==map[EndP2]),][1]
             }
            else  newEdge <- 13
            newEdge})
    }

    GetEdges <- local({
        Edges <- apply(CaseRotationFlip[-c(1,256),], 1, getedge)
        Case <- cbind(seq(1:256), CaseRotationFlip[,c(1,3)])
        Edges <- apply(Case[-c(1,256),], 1, function(x){
            case <- x[2]-1
            EdgeNum <- x[1]-1
            if (x[3]==1)
                sapply(EdgeSequence1[[case]], function(x) Edges[[EdgeNum]][x])
            else sapply(EdgeSequence2[[case]], function(x) Edges[[EdgeNum]][x])
        })
        Edges
    })

    BasicFace <- list(c(0),c(0),c(1),c(7),c(0),c(2,7),c(1,2,6,7),c(0),c(0),
                      c(5,6,7),c(0), c(1,4,7),c(1,2,3,4,5,6,7),c(0))
    FacePoints <- matrix(c(seq(1,6),1,2,4,1,1,5,6,7,7,8,3,7,2,3,3,4,2,
                           6,5,6,8,5,4,8),
                         ncol=5)
    FacePoints <- cbind(FacePoints, apply(FacePoints[,2:5],1,prod))

    getface <- function(x) {
        case <- x[1]
        rotation <- x[2]
        map <- rep(0,8)
        for(i in 1:8){
            temp <- as.integer(BasicRotation[rotation,][i])
            map[temp] <- i
        }
        sapply(BasicFace[[case-1]], function(x){
            EndP <- rep(0,4)
            if (x==0) newFace <- 0
            else if (x==7) newFace <- 7
            else {
                for (i in 1:4){
                     point <- FacePoints[x,i+1]
                     EndP[i] <- map[point]
                }
                newFace<- FacePoints[FacePoints[,6]==prod(EndP[1:4]),][1]
            }
            newFace})
    }

    flipface <- function(case, face){
        if (face!=0){
            index <- explode(case+1)
            if (sum(index) > 4)
                index <- ifelse(index==0,1,0)
            if (face!=7 && index[FacePoints[face,2]]==0)
                face <- -face
            else if (face==7){
                tcase <- CaseRotationFlip[case+1,1]-1
                if ((tcase == 4 || tcase==6 || tcase==10 ||tcase==12)
                    && !(index[1]+index[7]==2) && !(index[3]+index[5]==2))
                    face <- -face
                else if (tcase==7
                         && !(index[1]+index[7]==0) && !(index[3]+index[5]==0))
                    face <- -face
            }
        }
        face
    }

    GetFaces <- local({
        Faces <- apply(CaseRotationFlip[-c(1,256),], 1, getface)
        for (i in 1:254)
            for(j in 1:length(Faces[[i]]))
                Faces[[i]][j] <- flipface(i, Faces[[i]][j])
        Faces
    })

    special <- list(name =  c(3, 4, 6, 7,  10,12,13),
                    nface=  c(1, 1, 2, 4,  3, 3, 7),
                    sev =   c(0, 1, 1, 1,  1, 1, 1),
                    nedge = c(18,24,48,177,96,96,816),
                    ind = list(c(0,1), c(0,1), c(0,2,1,3),
                               c(0,8,4,12,2,10,1,9,6,14,5,13,3,11,15,7),
                               c(0,4,1,5,2,6,3,7),c(0,4,2,6,1,5,3,7),
                               c(0,1,2,4,8,16,32,3,9,17,33,6,18,34,12,20,36,24,
                                 40,35,25,22,44,19,41,38,28,83,105,102,92)),
                    position=list(list(c(1:6),c(7:18)),
                                  list(c(1:6),c(7:24)),
                                  list(c(1:9),c(10:33),c(34:48),c(34:48)),
                                  list(c(1:9),c(1:9),c(10:24),c(10:24),
                                       c(25:39),c(25:39),c(40:54), c(40:54),
                                       c(55:81), c(55:81),c(82:108),c(82:108),
                                       c(109:135),c(109:135),c(136:150),
                                       c(151:177)),
                                  list(c(1:12),c(13:36),c(37:60),c(37:60),
                                       c(61:84), c(61:84),c(85:96),c(85:96)),
                                  list(c(1:12), c(13:36), c(37:60), c(37:60),
                                       c(61:84), c(61:84),c(85:96),c(85:96)),
                                  list(c(1:12),
                                       c(13:30),c(31:48),c(49:66),c(67:84),
                                       c(85:102),c(103:120),
                                       c(121:150),c(151:180),c(181:210),
                                       c(211:240),c(241:270),c(271:300),
                                       c(301:330),
                                       c(331:360),c(361:390),c(391:420),
                                       c(421:450),c(451:480),
                                       c(481:516),c(517:552),c(553:588),
                                       c(589:624),
                                       c(625:642),c(643:660),c(661:678),
                                       c(679:696),
                                       c(697:726),c(727:756),c(757:786),
                                       c(787:816))
                                 )
                    )

    list(Edges = GetEdges, Faces = GetFaces, EdgePoints = EdgePoints,
         FacePoints = FacePoints, CaseRotationFlip = CaseRotationFlip,
         special = special)
})

Faces <- PreProcessing$Faces
Edges <- PreProcessing$Edges
EdgePoints <- PreProcessing$EdgePoints
FacePoints <- PreProcessing$FacePoints
CaseRotationFlip <- PreProcessing$CaseRotationFlip
special <- PreProcessing$special

fgrid <- function(fun, x, y, z) {
    g <- expand.grid(x = x, y = y, z = z)
    array(fun(g$x, g$y, g$z), c(length(x), length(y), length(z)))
}

faceType <- function(v, nx, ny, level, maxvol) {
    if(level==maxvol)
      p <- v >= level
    else p <- v > level
    v[p] <- 1; v[! p] <- 0
    v[-nx, -ny] + 2 * v[-1, -ny] + 4 * v[-1, -1] + 8 * v[-nx, -1]
}

levCells <- function(v, level, maxvol) {
    nx <- dim(v)[1]
    ny <- dim(v)[2]
    nz <- dim(v)[3]
    cells <- vector("list", nz - 1)
    types <- vector("list", nz - 1)

    bottomTypes <- faceType(v[,,1], nx, ny, level, maxvol)
    for (k in 1 : (nz - 1)) {
        topTypes <- faceType(v[,, k + 1], nx, ny, level, maxvol)
        cellTypes <- bottomTypes + 16 * topTypes
        contourCells <- which(cellTypes > 0 & cellTypes < 255)
        cells[[k]] <- contourCells + (nx - 1) * (ny - 1) * (k - 1)
        types[[k]] <- as.integer(cellTypes[contourCells])
        bottomTypes <- topTypes
    }
    cells <- unlist(cells)
    i <- as.integer((cells - 1) %% (nx - 1) + 1)
    j <- as.integer(((cells - 1) %/% (nx - 1)) %% (ny - 1) + 1)
    k <- as.integer((cells - 1) %/% ((nx - 1) * (ny - 1)) + 1)
    t <- unlist(types)
    list(i = i, j = j, k = k, t = t)
}

CalPoint <- function(x1,x2,y1,y2,z1,z2,v1,v2){
    s <- v1 / (v1-v2)
    x <- x1+s*(x2-x1)
    y <- y1+s*(y2-y1)
    z <- z1+s*(z2-z1)
    c(x,y,z)
}

GetPoints<-function(edge, p1, info){
    ##**** need better name than info
    ## info is the output from GetBasic()
    x1 <- EdgePoints[edge,2]
    x2 <- EdgePoints[edge,3]
    c((1-floor(x1/9))*info[p1+x1-1,1]+floor(x1/9)*info[p1,1],
      (1-floor(x1/9))*info[p1+x2-1,1]+floor(x1/9)*info[p1+1,1],
      (1-floor(x1/9))*info[p1+x1-1,2]+floor(x1/9)*info[p1+1,2],
      (1-floor(x1/9))*info[p1+x2-1,2]+floor(x1/9)*info[p1+2,2],
      (1-floor(x1/9))*info[p1+x1-1,3]+floor(x1/9)*info[p1+1,3],
      (1-floor(x1/9))*info[p1+x2-1,3]+floor(x1/9)*info[p1+5,3],
      (1-floor(x1/9))*info[p1+x1-1,4]+floor(x1/9)*(0*info[p1+1,3]+1),
      (1-floor(x1/9))*info[p1+x2-1,4]+floor(x1/9)*(0*info[p1+1,3]-1))
}

FaceNo7 <- function(faces, p1, info){
    ##**** need better name than info
    ## info is the output from GetBasic()
    index <- ifelse(faces > 0, 1, -1)
    faces <- abs(faces)
    e1 <- FacePoints[faces,2]
    e2 <- FacePoints[faces,3]
    e3 <- FacePoints[faces,4]
    e4 <- FacePoints[faces,5]
    A <- info[p1+e1-1,4]
    B <- info[p1+e2-1,4]
    C <- info[p1+e3-1,4]
    D <- info[p1+e4-1,4]
    index <- index*ifelse (A*B-C*D > 0, 1, -1)
    ifelse(index==1, 1, 0)
}

Face7 <- function(faces, p1, info){
    ## info is the output from GetBasic()
    index <- ifelse(faces > 0, 1, -1)
    A0 <- info[p1,4];   B0 <- info[p1+3,4]
    C0 <- info[p1+2,4]; D0 <- info[p1+1,4]
    A1 <- info[p1+4,4]; B1 <- info[p1+7,4]
    C1 <- info[p1+6,4]; D1 <- info[p1+5,4]
    a <- (A1 - A0)*(C1 - C0) - (B1 - B0)*(D1 - D0)
    b <- C0*(A1 - A0) + A0*(C1 - C0) - D0*(B1 - B0) - B0*(D1 - D0)
    c <- A0*C0 - B0*D0
    tmax <- -b/(2*a)
    maximum <- a*tmax^2 + b*tmax + c
    maximum <- ifelse(maximum=="NaN",-1,maximum)
    cond1 <- ifelse (a < 0, 1 ,0)
    cond2 <- ifelse (tmax > 0, 1 ,0)
    cond3 <- ifelse (tmax < 1, 1, 0)
    cond4 <- ifelse (maximum >0, 1, 0)
    totalcond <- cond1 * cond2 * cond3 * cond4
    index <- index*ifelse(totalcond==1, 1, -1)
    ifelse(index==1, 1, 0)
}


## GetBasic()--- The output matrix "information" consists of 4 columns
## and #-of-cubes*8 rows.  The first 3 columns tell the
## coordinate(x,y,z) of each vertex of the cube, the 4th gives the
## intensity minus the threshold, which actually makes the threshold
## eaqual to 0. This is convenient for further judgment of subcases

GetBasic <- function(R, vol, level, v) {
    cube.1 <- cbind(v$i[R], v$j[R], v$k[R])
    index <- matrix(c(0,1,1,0,0,1,1,0,
                      0,0,1,1,0,0,1,1,
                      0,0,0,0,1,1,1,1),
                    nrow=8)
    ax.inc <- c(1,1,1)

    ver.inc <- t(apply(index,1, function(x) x*ax.inc))
    cube.co <-
        kronecker(rep(1,nrow(cube.1)),ver.inc) + kronecker(cube.1,rep(1,8))

    value <- vol[cube.co] - level
    information <- cbind(cube.co, value)
    information <- rbind(information, rep(0, 4))
    p1 <- (1:length(R) - 1) * 8 + 1
    cases <- v$t[R]
    list(information=information, p1 = p1, cases=cases)
}

PreRender <- function(edges, p1, type, info) {
    if(type==1){
        if (typeof(edges)=="list"){
            count <- sapply(edges, function(x) length(x))
            edges <- cbind(unlist(edges), rep(p1,count))
        }
        else{
            count <- nrow(edges)
            edges <- cbind(as.vector(t(edges)), rep(p1,each=count))
        }
    }
    else{
        if (is.vector(edges))
            edges <- matrix(edges, ncol = length(edges))
        p1 <- edges[, 1]
        count <- ncol(edges) - 1
        edges <- cbind(as.vector(t(edges[, -1])), rep(p1, each = count))
    }
    ##The output of GetPoints() are coordinates of cubes.
    info <- GetPoints(edges[,1],edges[,2], info)
    info <- matrix(info,ncol=8)
    ##The output of CalPoint() are coordinates of triangles.
    info <- CalPoint(info[,1],info[,2],info[,3],info[,4],
                            info[,5],info[,6],info[,7],info[,8])
    matrix(info,ncol=3)
}

rescale <- function(i, x) {
    nx <- length(x)
    low <- pmin(pmax(1, floor(i)), nx - 1)
    x[low] + (i - low) * (x[low + 1] - x[low])
}

computeContour3d <- function (vol, maxvol, level,
                              x = 1:dim(vol)[1],
                              y = 1:dim(vol)[2],
                              z = 1:dim(vol)[3], mask) {

    nx <- length(x)
    ny <- length(y)
    nz <- length(z)

    if (is.function(mask)) mask <- fgrid(mask, x, y, z)
    if (! all(mask)) vol[! mask] <- NA

    v <- levCells(vol, level, maxvol)
    tcase <- CaseRotationFlip[v$t+1,1]-1

    R <- which(tcase %in% c(1,2,5,8,9,11,14))
    if (length(R) > 0){
        Basics <- GetBasic(R, vol, level, v)
        information <- Basics$information
        p1 <- Basics$p1
        cases <- Basics$cases
        edges <- Edges[cases]
        triangles <- PreRender(edges, p1,type=1, information)
    }
    else triangles <- matrix(0, nrow=0,ncol=3) # emty contour, e.g.

    for (i in 1:length(special$name)){
        R <- which(tcase == special$name[i])
        if (length(R) > 0) {
            Basics <- GetBasic(R, vol, level, v)
            information <- Basics$information
            p1 <- Basics$p1
            cases <- Basics$cases

            nface <- special$nface[i]
            nedge <- special$nedge[i]
            faces <- matrix(unlist(Faces[cases]), ncol = nface, byrow = TRUE)

            if (i==1)
                index <- FaceNo7(faces[, 1], p1, information)
            else if (i==2)
                index <- Face7(faces[, 1], p1, information)
            else{
                index <-  Face7(faces[, nface], p1, information)*2^(nface-1)
                for(j in 1:(nface-1)){
                    temp <-  FaceNo7(faces[, j], p1, information)
                    index <- index + temp * 2^(j-1)
                }
            }
            edges <- matrix(unlist(Edges[cases]), ncol = nedge, byrow = TRUE)
            edges <- cbind(edges, p1, index)
            ind <- special$ind[[i]]
            position <- special$position[[i]]

            for (j in 1:length(ind)){
                ed <- edges[which(index == ind[j]), c(nedge+1, position[[j]])]
                if (length(ed) > 0) {
                    prtri <- PreRender(ed,nedge+1,type=2, information)
                    triangles <- rbind(triangles, prtri)
                }
            }
        }
    }

    if (! identical(x, 1 : nx)) triangles[,1] <- rescale(triangles[,1], x)
    if (! identical(y, 1 : ny)) triangles[,2] <- rescale(triangles[,2], y)
    if (! identical(z, 1 : nz)) triangles[,3] <- rescale(triangles[,3], z)

    triangles
}

################################################################################
contourTriangles <- function(vol, maxvol, level,
                             x = 1:dim(vol)[1],
                             y = 1:dim(vol)[2],
                             z = 1:dim(vol)[3],
                             mask = NULL, color = "white", color2 = NA,
                             alpha = 1, fill = TRUE,
                             col.mesh = if (fill) NA else color,
                             material = "default", smooth = 0) {
    if (length(level) > 1) {
        val <- vector("list", length(level))
        for (i in seq(along = level)) {
            m <- if (is.list(mask)) mask[[i]] else mask
            col <- if (length(color) > 1) color[[i]] else color
            col2 <- if (length(color2) > 1) color2[[i]] else color2
            a <- if (length(alpha) > 1) alpha[[i]] else alpha
            fl <- if (length(fill) > 1) fill[[i]] else fill
            cm <- if (length(col.mesh) > 1) col.mesh[[i]] else col.mesh
            mat <- if (length(material) > 1) material[[1]] else material
            sm <- if (length(smooth) > 1) smooth[[1]] else smooth
            val[[i]] <- contourTriangles(vol, maxvol, level[i], x, y, z, m,
                                         col, col2, a, fl, cm, mat, sm)
        }
        val
    }
    else makeTriangles(computeContour3d(vol, maxvol, level, x, y, z, mask),
                       color = color, color2 = color2, alpha = alpha,
                       fill = fill, col.mesh = col.mesh,
                       material = material, smooth = smooth)
}

################################################################################
getlevel <- function(f, x, y, z) {
    if (! all(is.finite(x), is.finite(y), is.finite(z)))
        stop("'x', 'y', and 'z' values must be finite and non-missing")
    if (is.function(f) || is.array(f) && length(dim(f)) == 3){
        if (is.function(f)){
            if (length(formals(f)) < 3)
                stop("The function must have at least 3 arguments.")
            vol <- fgrid(f, x, y, z)
          }
        else{
          if (dim(f)[1] != length(x) || dim(f)[2] != length(y) ||  dim(f)[3] != length(z))
            stop("dimensions of f do not match x, y, or z")
          vol <- f
        }

        maxvol <- max(vol)
        minvol <- min(vol)
        print(paste("Your level is betwen ", minvol, " and ", maxvol, sep=""),quote=F)
    }
}


################################################################################
makeTriangles <- function(v1, v2, v3,
                          color = "red", color2 = NA, alpha = 1,
                          fill = TRUE, col.mesh = if (fill) NA else color,
	                  smooth = 0,  material = "default") {
    if (missing(v2) || missing(v3)) {
        if (missing(v2) && missing(v3))
  	    v <- unzipTriangleMatrix(v1)
        else if (missing(v3))
            v <- ve2t(list(vb = v1, ib = v2))
        else stop("unknown form of triangle specification")
	v1 <- v$v1
	v2 <- v$v2
	v3 <- v$v3
    }
    tris <- structure(list(v1 = v1, v2 = v2, v3 = v3,
                           color = color, color2 = color2, fill = fill,
                           material = material, col.mesh = col.mesh,
                           alpha = alpha, smooth = smooth),
                      class = "Triangles3D")
    colorTriangles(tris)
}

is.Triangles3D <- function(x) identical(class(x), "Triangles3D")

updateTriangles <- function(triangles, color, color2, alpha, fill, col.mesh,
                            material, smooth) {
    if (! missing(color)) triangles$color <- color
    if (! missing(color2)) triangles$color2 <- color2
    if (! missing(fill)) triangles$fill <- fill
    if (! missing(col.mesh)) triangles$col.mesh <- col.mesh
    if (! missing(material)) triangles$material <- material
    if (! missing(alpha)) triangles$alpha <- alpha
    if (! missing(smooth)) triangles$smooth <- smooth
    colorTriangles(triangles)
}

#**** This assumes comparable scaling of dimensions
#**** 5 is the largest exponent for S that will work; smaller is OK
t2ve <- function (triangles)
{
    vb <- rbind(triangles$v1, triangles$v2, triangles$v3)
    vbmin <- min(vb)
    vbmax <- max(vb)
    S <- 10^5
    score <- function(v, d) floor(as.vector(v %*% d))
    scale <- function(v) (1 - 1 / S) * (v - vbmin) / (vbmax - vbmin)
    d <- c(S, S^2, S^3)
    scores <- score(scale(vb), d)
    vb <- vb[! duplicated(scores),]
    scores <- score(scale(vb), d)
    ib <- rbind(match(score(scale(triangles$v1), d), scores),
                match(score(scale(triangles$v2), d), scores),
                match(score(scale(triangles$v3), d), scores))
    list(vb = t(vb), ib = ib)
}

ve2t <- function(ve) {
    list (v1 = t(ve$vb[,ve$ib[1,]]),
          v2 = t(ve$vb[,ve$ib[2,]]),
          v3 = t(ve$vb[,ve$ib[3,]]))
}

unzipTriangleMatrix <- function(tris) {
    if (ncol(tris) != 3)
        stop("triangle matrix must have three columns.")
    if (nrow(tris) %% 3 != 0)
        stop("number of rows in triangle matrix must be divisible by 3")
    n <- nrow(tris) / 3
    list(v1 = tris[3 * (1 : n) - 2,],
         v2 = tris[3 * (1 : n) - 1,],
         v3 = tris[3 * (1 : n),])
}

zipTriangles <- function(tris) {
    n <- nrow(tris$v1)
    if(is.null(n)) n <- 1
    else if (nrow(tris$v2) != n || nrow(tris$v3) != n)
        stop("vertex arrays must have the same number of rows")
    v <- matrix(0, nrow = 3 * n, ncol = 3)
    v[3 * (1 : n) - 2,] <- tris$v1
    v[3 * (1 : n) - 1,] <- tris$v2
    v[3 * (1 : n),] <- tris$v3
    v
}

colorTriangles <- function(triangles) {
    if (is.function(triangles$color) || is.function(triangles$color2)) {
        v <- (triangles$v1 + triangles$v2 + triangles$v3) / 3
        if (is.function(triangles$color))
            triangles$color <- triangles$color(v[,1], v[,2], v[,3])
        if (is.function(triangles$color2))
            triangles$color2 <- triangles$color2(v[,1], v[,2], v[,3])
        if (is.function(triangles$col.mesh))
            triangles$col.mesh <- triangles$col.mesh(v[,1], v[,2], v[,3])
    }
    triangles
}

colorScene <- function(scene) {
    if (is.Triangles3D(scene))
        colorTriangles(scene)
    else lapply(scene, colorTriangles)
}

## **** better to make new triangles including only requested components?
canonicalizeAndMergeScene <- function(scene, ...) {
    which <- list(...)
    if (is.Triangles3D(scene)) {
        #cat("scene v1 is ", print(scene$v1), "\n")
        n.tri <- nrow(scene$v1)
        if(is.null(n.tri))  n.tri<-1
        for (n in which)
            if (length(scene[[n]]) != n.tri)
                scene[[n]] <- rep(scene[[n]], length = n.tri)
        scene
    }
    else {
        scene <- lapply(scene, canonicalizeAndMergeScene, ...)
        x <- scene[[1]]
        x$v1 <- do.call(rbind, lapply(scene, function(x) x$v1))
        x$v2 <- do.call(rbind, lapply(scene, function(x) x$v2))
        x$v3 <- do.call(rbind, lapply(scene, function(x) x$v3))
        for (n in which)
            x[[n]] <- do.call(c, lapply(scene, function(x) x[[n]]))
        x
    }
}

expandTriangleGrid <- function(x, y) {
    nx <- length(x) - 1
    ny <- length(y) - 1
    A <- c(0, 0)
    B <- c(1, 0)
    C <- c(1, 1)
    D <- c(0, 1)
    g <- expand.grid(x = 1 : nx, y = 1 : ny)
    even <- (g$x + g$y) %% 2 == 0
    gx11 <- ifelse(even, g$x + A[1], g$x + A[1])
    gy11 <- ifelse(even, g$y + A[2], g$y + A[2])
    gx12 <- ifelse(even, g$x + A[1], g$x + B[1])
    gy12 <- ifelse(even, g$y + A[2], g$y + B[2])
    i1 <- rbind(cbind(gx11, gy11), cbind(gx12, gy12))
    gx21 <- ifelse(even, g$x + B[1], g$x + B[1])
    gy21 <- ifelse(even, g$y + B[2], g$y + B[2])
    gx22 <- ifelse(even, g$x + C[1], g$x + C[1])
    gy22 <- ifelse(even, g$y + C[2], g$y + C[2])
    i2 <- rbind(cbind(gx21, gy21), cbind(gx22, gy22))
    gx31 <- ifelse(even, g$x + C[1], g$x + D[1])
    gy31 <- ifelse(even, g$y + C[2], g$y + D[2])
    gx32 <- ifelse(even, g$x + D[1], g$x + D[1])
    gy32 <- ifelse(even, g$y + D[2], g$y + D[2])
    i3 <- rbind(cbind(gx31, gy31), cbind(gx32, gy32))
    v1 <- cbind(x[i1[,1]], y[i1[,2]])
    v2 <- cbind(x[i2[,1]], y[i2[,2]])
    v3 <- cbind(x[i3[,1]], y[i3[,2]])
    list(v1 = v1, v2 = v2, v3 = v3)
}

## adapted from lattice ltransform3dto3d
trans3dto3d <- function (x, R.mat) {
    if (length(x) == 0)
        return(x)
    val <- R.mat %*% rbind(t(x), 1)
    val[1, ] <- val[1, ]/val[4, ]
    val[2, ] <- val[2, ]/val[4, ]
    val[3, ] <- val[3, ]/val[4, ]
    t(val[1:3, , drop = FALSE])
}

transformTriangles <- function(triangles, R) {
    tr <- function(v) trans3dto3d(v, R)
    triangles$v1 <- tr(triangles$v1)
    triangles$v2 <- tr(triangles$v2)
    triangles$v3 <- tr(triangles$v3)
    triangles
}

transformScene <- function(scene, rot.mat) {
    if (is.Triangles3D(scene))
        transformTriangles(scene, rot.mat)
    else lapply(scene, transformTriangles, rot.mat)
}

translateTriangles <- function(triangles, x = 0, y = 0, z = 0) {
    M <- diag(4)
    M[1:3,4] <- c(x, y, z)
    transformTriangles(triangles, M)
}

scaleTriangles <- function(triangles, x = 1, y = x, z = x) {
    M <- diag(c(x, y, z, 1))
    transformTriangles(triangles, M)
}

## triangleNormals computes the normal vectors to a collection of
## triangles as the vector crossprocuct of the direction from v1 to v2
## and the direction from v2 to v3.  The result is an n by 3 matrix of
## unit representing the n unit normal vectors.

triangleNormals <- function(triangles) {
   x <- triangles$v2 - triangles$v1
   y <- triangles$v3 - triangles$v2
   z <- cbind(x[,2]*y[,3] - x[,3]*y[,2],
              x[,3]*y[,1] - x[,1]*y[,3],
              x[,1]*y[,2] - x[,2]*y[,1])
   z / sqrt(rowSums(z^2))
}

# adapted from lattice ltransform3dMatrix
trans3dMat <- function (screen, P = diag(4)) {
    givens4 <- function(i, j, gamma) {
        T <- diag(4)
        cgamma <- cos(gamma)
        sgamma <- sin(gamma)
        T[c(i,j),c(i,j)] <- matrix(c(cgamma, sgamma, -sgamma, cgamma), 2, 2)
        T
    }
    screen.names <- names(screen)
    for (i in seq(along = screen.names)) {
        if (screen.names[i] == "x")
            P <- givens4(2, 3, screen[[i]] * pi/180) %*% P
        else if (screen.names[i] == "y")
            P <- givens4(1, 3, -screen[[i]] * pi/180) %*% P #**** whi negative?
        else if (screen.names[i] == "z")
            P <- givens4(1, 2, screen[[i]] * pi/180) %*% P
    }
    P
}

makeViewTransform <- function(ranges, scale, aspect, screen, R.mat) {
    m <- c(mean(ranges$xlim), mean(ranges$ylim), mean(ranges$zlim))
    s <- 0.5 * c(diff(ranges$xlim), diff(ranges$ylim), diff(ranges$zlim))
    if (! scale) s <- rep(max(s), 3)
    else s <- s / c(1, aspect)
    A <- diag(1 / c(s, 1))
    A[1:3, 4] <- -m / s
    trans3dMat(screen, R.mat %*% A)
}

trianglesRanges <- function(triangles, xlim, ylim, zlim) {
    v1 <- triangles$v1
    v2 <- triangles$v2
    v3 <- triangles$v3
    if (is.null(xlim)) xlim <- range(v1[,1], v2[,1], v3[,1], na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(v1[,2], v2[,2], v3[,2], na.rm = TRUE)
    if (is.null(zlim)) zlim <- range(v1[,3], v2[,3], v3[,3], na.rm = TRUE)
    list(xlim = xlim, ylim = ylim, zlim = zlim)
}

sceneRanges <- function(scene, xlim, ylim, zlim) {
    if (is.Triangles3D(scene))
        trianglesRanges(scene, xlim, ylim, zlim)
    else {
        ranges <- lapply(scene, trianglesRanges, xlim, ylim, zlim)
        list(xlim = range(sapply(ranges,function(x) x$xlim)),
             ylim = range(sapply(ranges,function(x) x$ylim)),
             zlim = range(sapply(ranges,function(x) x$zlim)))
    }
}

addTrianglesPerspective <- function(triangles, distance) {
    pt <- function(v) {
        v[, 1] <- v[, 1] / (1 / distance - v[, 3])
        v[, 2] <- v[, 2] / (1 / distance - v[, 3])
        v
    }
    triangles$v1 <- pt(triangles$v1)
    triangles$v2 <- pt(triangles$v2)
    triangles$v3 <- pt(triangles$v3)
    triangles
}

addPerspective <- function(scene, distance) {
    if (is.Triangles3D(scene))
        addTrianglesPerspective(scene, distance)
    else lapply(scene, addTrianglesPerspective, distance)
}

screenRange <- function(v1, v2, v3)
    range(v1[,1:2], v2[,1:2], v3[,1:2], na.rm = TRUE)

vertexTriangles <- function(ve) {
    n.vert <- ncol(ve$vb)
    ib <- ve$ib
    vt <- function(i) which(ib[1,] == i | ib[2,] == i | ib[3,] == i)
    lapply(1 : n.vert, vt)
}

# faster version
vertexTriangles <- function(ve) {
    n.vert <- ncol(ve$vb)
    val <- vector("list", n.vert)
    ib <- ve$ib
    for (i in 1 : ncol(ib)) {
        val[[ib[1,i]]] <- c(val[[ib[1,i]]], i)
        val[[ib[2,i]]] <- c(val[[ib[2,i]]], i)
        val[[ib[3,i]]] <- c(val[[ib[3,i]]], i)
    }
    val
}


vertexNormals <- function(vt, N) {
    vn <- function(tris) {
        z <- apply(N[tris,,drop = FALSE], 2, mean, na.rm = TRUE);
        z <- z / sqrt(sum(z^2))
        if (any(is.na(z))) c(1,0,0) else z
    }
    t(sapply(vt, vn))
}

# faster version
vertexNormals <- function(vt, N) {
    val <- matrix(0, nrow = length(vt), ncol = 3)
    for (i in seq(along = vt)) {
        Ni <- N[vt[[i]],,drop = FALSE]
        Ni1 <- Ni[,1]
        Ni2 <- Ni[,2]
        Ni3 <- Ni[,3]
        z1 <- if (any(is.na(Ni1))) mean(Ni1, na.rm = TRUE)
              else sum(Ni1) / length(Ni1)
        z2 <- if (any(is.na(Ni2))) mean(Ni2, na.rm = TRUE)
              else sum(Ni2) / length(Ni2)
        z3 <- if (any(is.na(Ni3))) mean(Ni3, na.rm = TRUE)
              else sum(Ni3) / length(Ni3)
        z <- c(z1, z2, z3)
        z <- z / sqrt(sum(z^2))
        val[i,] <- if (any(is.na(z))) c(1,0,0) else z
    }
    val
}

interpolateVertexNormals <- function(VN, ib) {
    z <- (VN[ib[1,],] + VN[ib[2,],] + VN[ib[3,],]) / 3
    z / sqrt(rowSums(z^2))
}

## triangleVertexNormals computes the normals at the vertices by
## averaging the normals of the incident triangles.  This is used by
## the rgl engine.  The result form is chosen so zipTriangles can be
## used on it.
triangleVertexNormals <- function(v) {
    N <- triangleNormals(v)
    ve <- t2ve(v)
    vt <- vertexTriangles(ve)
    VN <- vertexNormals(vt, N)
    list(v1 = VN[ve$ib[1,],], v2 = VN[ve$ib[2,],], v3 = VN[ve$ib[3,],])
}

vertexColors <- function(vt, col) {
    C <- t(col2rgb(col))
    val <- matrix(0, nrow = length(vt), ncol = 3)
    for (i in seq(along = vt)) {
        vti <- vt[[i]]
        nti <- length(vti)
        Ci <- C[vti,,drop = FALSE]
        Ci1 <- Ci[,1]
        Ci2 <- Ci[,2]
        Ci3 <- Ci[,3]
        val[i,] <- c(sum(Ci1), sum(Ci2), sum(Ci3)) / nti
    }
    val
}

interpolateVertexColors <- function(VC, ib) {
    TC <- (VC[ib[1,],] + VC[ib[2,],] + VC[ib[3,],]) / 3
    rgb(TC[,1], TC[,2], TC[,3], maxColorValue=255)
}

triangleEdges <- function(vb, ib) {
    edges <- cbind(ib[c(1,2),], ib[c(2,3),], ib[c(3,1),])
    swap <- edges[1,] > edges[2,]
    edges[,swap] <- edges[2:1,swap]
    edges[,! duplicated(edges, MARGIN = 2)]
}

# faster version
triangleEdges <- function(vb, ib) {
    n.vert <- ncol(vb)
    edges <- cbind(ib[c(1,2),], ib[c(2,3),], ib[c(3,1),])
    swap <- edges[1,] > edges[2,]
    edges[,swap] <- edges[2:1,swap]
    score <- as.vector(c(1 + n.vert, 1) %*% edges)
    edges[,! duplicated(score)]
}

triangleMidTriangles <- function(vb, ib, VN) {
    n.vert <- ncol(vb)
    edges <- triangleEdges(vb, ib)
    vb <- (vb[,edges[1,]] + vb[,edges[2,]]) / 2
    d <- c(1 + n.vert, 1)
    scores <- as.vector(d %*% edges)
    mpi <- function(a, b) {
        s <- d[1] * pmin(a, b) + d[2] * pmax(a, b)
        match(s, scores)
    }
    mpi1 <- mpi(ib[1,], ib[2,])
    mpi2 <- mpi(ib[2,], ib[3,])
    mpi3 <- mpi(ib[3,], ib[1,])
    ib <- rbind(mpi1, mpi2, mpi3)
    z <- VN[edges[1,],] + VN[edges[2,],]
    z <- z / sqrt(rowSums(z^2))
    list(vb = vb, ib = ib, VN = z)
}

## surfaceTriangles creates a set of triangles for a grid specified by x,
## y and function falues computed with f if f is a function or taken
## from f if f is a matrix.

surfaceTriangles <- function(x, y, f,
                             color = "red", color2 = NA,  alpha = 1,
                             fill = TRUE, col.mesh = if (fill) NA else color,
                             smooth = 0, material = "default") {
    if (is.function(f))
        ff <- function(ix, iy) f(x[ix], y[iy])
    else
        ff <- function(ix, iy) f[ix + length(x) * (iy - 1)]
    i <- expandTriangleGrid(1 : length(x), 1 : length(y))
    i1 <- i$v1
    i2 <- i$v2
    i3 <- i$v3
    v1 <- cbind(x[i1[,1]], y[i1[,2]], ff(i1[,1], i1[,2]))
    v2 <- cbind(x[i2[,1]], y[i2[,2]], ff(i2[,1], i2[,2]))
    v3 <- cbind(x[i3[,1]], y[i3[,2]], ff(i3[,1], i3[,2]))
    na1 <- is.na(v1[,1]) | is.na(v1[,2]) | is.na(v1[,3])
    na2 <- is.na(v2[,1]) | is.na(v2[,2]) | is.na(v2[,3])
    na3 <- is.na(v3[,1]) | is.na(v3[,2]) | is.na(v3[,3])
    nna <- ! (na1 | na2 | na3)
    makeTriangles(v1[nna,], v2[nna,], v3[nna,],
                  color = color, color2 = color2, fill = fill, smooth = smooth,
                  material = material, col.mesh = col.mesh, alpha = alpha)
}

## pointsTetrahedra computes a collection of tetrahedra centered at
## the specified point locations.  This is useful, for example, for
## displaying raw data along with a density contour in a scene
## rendered with standard or grid graphics. Random orientation might
## be useful to avoid strange results at certain lighting angles.

pointsTetrahedra <- function(x, y, z, size = 0.01, color = "black", ...) {
    n <- length(x)
    if (length(y) != n || length(z) != n)
        stop("coordinate vectors must be the same length.")

    ## Create a basic tetrahedron centered at the origin
    a <- sqrt(3) / 2
    b <- 1 / (2 * sqrt(3))
    h <- sqrt(2 / 3)

    mx <- 1 / 2
    my <- (a + b) / 4
    mz <- h / 4

    A <- c(       -mx,    -my,    -mz)
    B <- c(    1 - mx,    -my,    -mz)
    C <- c(1 / 2 - mx, a - my,    -mz)
    D <- c(1 / 2 - mx, b - my, h - mz)

    v1 <- rbind(B, A, B, C)
    v2 <- rbind(A, B, C, A)
    v3 <- rbind(C, D, D, D)

    ## Scale the tetrahedron
    if (length(size) < 3) size <- rep(size, len = 3)
    if (n == 1) s <- diag(size)
    else s <- diag(size * c(diff(range(x)), diff(range(y)), diff(range(z))))
    sv1 <- v1 %*% s
    sv2 <- v2 %*% s
    sv3 <- v3 %*% s

    ## Compute the tetrahedra for the points, taking advantage of recycling
    x4 <- rep(x, each = 4)
    y4 <- rep(y, each = 4)
    z4 <- rep(z, each = 4)

    V1 <- cbind(x4 + sv1[,1], y4 + sv1[,2], z4 + sv1[,3])
    V2 <- cbind(x4 + sv2[,1], y4 + sv2[,2], z4 + sv2[,3])
    V3 <- cbind(x4 + sv3[,1], y4 + sv3[,2], z4 + sv3[,3])

    makeTriangles(V1, V2, V3, color = color, ...)
}


## Compute for each triangle the indices of triangles that share an
## edge with it.  This could be done more efficiently.
triangleNeighbors <- function(tris) {
   ve <- t2ve(tris)
   vt <- vertexTriangles(ve)
   ib <- ve$ib
   n.tri <- ncol(ib)
   tn <- vector("list", n.tri)
   for (i in 1 : n.tri) {
       v1 <- unique(vt[[ib[1, i]]])
       v2 <- unique(vt[[ib[2, i]]])
       v3 <- unique(vt[[ib[3, i]]])
       i12 <- intersect(v1, v2)
       i23 <- intersect(v2, v3)
       i31 <- intersect(v3, v1)
       u <- union(union(i12, i23), i31)
       tn[[i]] <- u[u != i]
   }
   tn
}
## 'unique' in unique(vt[[ib[1, i]]]) seems to be unnecessary
## unless a triangle has essentially two vertices or one vertex
triangleNeighbors <- function(tris) {
   ve <- t2ve(tris)
   vt <- vertexTriangles(ve)
   ib <- ve$ib
   n.tri <- ncol(ib)
   tn <- vector("list", n.tri)
   for (i in 1 : n.tri) {
       v1 <- vt[[ib[1, i]]]
       v2 <- vt[[ib[2, i]]]
       v3 <- vt[[ib[3, i]]]
       i12 <- intersect(v1, v2)
       i23 <- intersect(v2, v3)
       i31 <- intersect(v3, v1)
       u <- union(union(i12, i23), i31)
       tn[[i]] <- u[u != i]
   }
   tn
}

## Dijkstra's version of Rem's algorithm for computing equivalence
## classes based on a number of vertices 1:nvert and a set of N edges
## provided as an N x 2 matrix.
GetPatches <- function(nvert, edges) {

   f <- 1:nvert

   if (!(is.vector(edges)) && dim(edges)[1] != 0){
       nedge <- nrow(edges)

       for (e in 1:nedge) {
           p0 <- edges[e, 1]
           q0 <- edges[e, 2]
           p1 <- f[p0]
           q1 <- f[q0]
           while (p1 != q1) {
               if (q1 < p1) {
                   f[p0] <- q1
                   p0 <- p1
                   p1 <- f[p1]
               }
               else {
                   f[q0] <- p1
                   q0 <- q1
                   q1 <- f[q1]
               }
           }
       }
   }
   if(is.vector(edges)){
       if(edges[1] < edges[2])
           f[edges[2]] <- edges[1]
       else  f[edges[1]] <- edges[2]
   }

   for (v in 1:nvert)
       f[v] <- f[f[v]]

   split(1:nvert,f)
}

## compute the edges to indicate which triangles share an edge -- this
## needs more error checking
triangleNeighborEdges <- function(tn) {
   edges <- function(i) {
       v <- tn[[i]]
       if (length(v) > 0) cbind(i,v)
       else numeric(0)
   }
   do.call(rbind, lapply(1:length(tn), edges))
}

## separate triangles into disconnected chunks
separateTriangles <- function(contour3dObj){
    tn <- triangleNeighbors(contour3dObj)
    edges <- triangleNeighborEdges(tn)
    edges <- edges[edges[,1] < edges[,2],]
    p <- GetPatches(length(tn), edges)
    newContour3dObj <- vector("list", length(p))
    for(i in 1:length(newContour3dObj)){
        newContour3dObj[[i]] <- contour3dObj
        newContour3dObj[[i]]$v1 <- contour3dObj$v1[p[[i]],]
        newContour3dObj[[i]]$v2 <- contour3dObj$v2[p[[i]],]
        newContour3dObj[[i]]$v3 <- contour3dObj$v3[p[[i]],]
    }
    newContour3dObj
}

drawScene.rgl <- function(scene, add = FALSE, ...) {
    loadRGL()
    if (! rgl.cur())
        open3d()
    if (!add)
        clear3d()

    scene <- colorScene(scene)
    triangles <- canonicalizeAndMergeScene(scene, "color", "color2", "alpha",
                                           "col.mesh", "fill", "smooth")
    use.col2 <- ! all(is.na(triangles$color2))
    if (use.col2 && any(is.na(triangles$color2)))
        warning(paste("mixing surfaces with and without color2 may not",
                      "work properly in the rgl engine"))

    col <- rep(triangles$color, each = 3)
    if (use.col2)
    	col2 <- rep(triangles$color2, each=3)
    alpha <- rep(triangles$alpha, each = 3)
    fill <- rep(triangles$fill, each = 3)
    col.mesh <- rep(triangles$col.mesh, each = 3)
    data <- zipTriangles(triangles)
    if (all(fill)) {
        front <- "filled"
        back <- "filled"
    }
    else if (any(fill))
        ##**** handle these by splitting; OK if no alpha < 1
        stop(paste("for now rgl engine cannot handle mixed fill/wire",
                   "frame surfaces"))
    else {
        front <- "lines"
        back <- "lines"
        col <- col.mesh
    }
    oldstyle = getOption("old.misc3d.orientation")
    if (! is.null(oldstyle) && oldstyle) {
        data <- data[,c(1, 3, 2)]
        data[,3] <- -data[,3]
    }
    if (any(triangles$smooth > 0)) {
        if (any(triangles$smooth == 0))
            stop(paste("for now for the rgl engine cannot handle mixed",
                       "smooth/non-smooth surfaces"))
        normals <- zipTriangles(triangleVertexNormals(triangles))
    }
    else normals <- NULL
    if (nrow(data) > 0) # to avoid a segfault in rgl
    {
        if (! use.col2)
    	    triangles3d(data[,1], data[,2], data[,3],
                      col = col, alpha = alpha, normals = normals,
                      front = front, back = back, ...)
        else {
            triangles3d(data[,1], data[,2], data[,3],
	              col = col, alpha = alpha, normals = normals,
                      front = front, back = "cull", ...)
            triangles3d(data[,1], data[,2], data[,3],
	              col = col2, alpha = alpha, normals = normals,
                      front = "cull", back = back, ...)
        }
    }
        
}

renderScene <- function(scene, box, fill, col.mesh, add, engine, polynum,
                        col.bg, depth, newpage) {
    triangles <- canonicalizeAndMergeScene(scene, "color", "col.light",
                                           "col.mesh", "fill")
    v1 <- triangles$v1
    v2 <- triangles$v2
    v3 <- triangles$v3
    n.tri <- nrow(v1)
    if (fill) {
        fill <- rep(triangles$fill, length = n.tri)
        col.mesh <- rep(triangles$col.mesh, length = n.tri)
    }
    else col.mesh <- rep(col.mesh, length = n.tri)
    col.fill <- ifelse(fill, triangles$col.light, NA)
    z <- (v1[,3] + v2[,3] + v3[,3]) / 3
    if (depth > 0) {
        rgbcol <- col2rgb(col.fill, alpha = TRUE)/255
        rgbcol.bg <- col2rgb(col.bg, alpha = TRUE)/255
        s <- (1 + depth * z) / (1 + depth * max(z))
        col.fill <- rgb(rgbcol[1,] * s + rgbcol.bg[1,] * (1 - s),
                        rgbcol[2,] * s + rgbcol.bg[2,] * (1 - s),
                        rgbcol[3,] * s + rgbcol.bg[3,] * (1 - s),
                        rgbcol[4,])
    }
    col.mesh <- ifelse(is.na(col.mesh), col.fill, col.mesh)
    i <- order(z, na.last = NA)
    if (engine == "grid") 
        render.grid(v1[i,], v2[i,], v3[i,], box, fill[i], col.fill[i],
                    col.mesh[i], add, polynum, newpage)
    else render.standard(v1[i,], v2[i,], v3[i,], box, fill[i], col.fill[i],
                         col.mesh[i], add)
}

render.standard <- function(v1, v2, v3, box, fill, col.fill, col.mesh, add) {
    if (! add) {
        # rr <- screenRange(v1, v2, v3)
        rr <- range(box)
        plot(rr, rr,type="n", axes = FALSE, ann = FALSE)
    }
    xx <- as.vector(rbind(v1[,1], v2[,1], v3[,1], NA))
    yy <- as.vector(rbind(v1[,2], v2[,2], v3[,2], NA))
    polygon(xx, yy, col=col.fill, border=col.mesh)
}

render.grid <- function(v1, v2, v3, box, fill, col.fill, col.mesh,
                        add, polynum, newpage) {
    if (! add) {
        if (newpage) grid::grid.newpage()
        rr <- range(box)
        grid::pushViewport(grid::viewport(w = 0.8, h = 0.8,
                                          xscale = rr, yscale = rr,
                                          name = "misc3dScene"))
        on.exit(grid::upViewport())
    }

    xx <- as.vector(rbind(v1[,1], v2[,1], v3[,1]))
    yy <- as.vector(rbind(v1[,2], v2[,2], v3[,2]))
    n.tri <- nrow(v1)
    idlen <- rep(3, n.tri)
    start <- 1
    end <- start + polynum - 1
    while (start <= n.tri) {
        end <- min(end, n.tri)
        j <- start : end
        j3 <- (3 * start - 2) : (3 * end)
        gp <- grid::gpar(fill = col.fill[j], col = col.mesh[j])
        grid::grid.polygon(x = xx[j3], y = yy[j3], default.units = "native",
                           gp = gp, id.lengths = idlen[j])
        start <- start + polynum
        end <- start + polynum
    }
}

makePerspMatrix <- function(d) {
    rbind(c(1, 0,  0, 0),
          c(0, 1,  0, 0),
          c(0, 0,  1, 0),
          c(0, 0, -1, 1 / d))
}


drawScene <- function(scene, light = c(0, 0, 1),
                      screen = list(z = 40, x = -60),
                      scale = TRUE,
                      R.mat = diag(4),
                      perspective = FALSE,
                      distance = if (perspective) 0.2 else 0, 
                      fill = TRUE,
                      xlim = NULL, ylim = NULL, zlim = NULL,
                      aspect = c(1, 1),
                      col.mesh = if (fill) NA else "black",
                      polynum = 100,
                      lighting = phongLighting,
                      add = FALSE,
                      engine = "standard",
                      col.bg = "transparent", depth = 0, newpage = TRUE) {
    scene <- colorScene(scene)
    sr <- sceneRanges(scene, xlim, ylim, zlim)
    if (add)
        rot.mat <- R.mat
    else
        rot.mat <- makeViewTransform(sr, scale, aspect, screen, R.mat)
    scene <- transformScene(scene, rot.mat)
    scene <- lightScene(scene, lighting, light)
    if (distance > 0) {
        scene <- addPerspective(scene, distance)
        rot.mat <- makePerspMatrix(distance) %*% rot.mat
    }
    box <- as.matrix(expand.grid(sr$xlim, sr$ylim, sr$zlim))
    box <- trans3dto3d(box, rot.mat)
    renderScene(scene, box, fill, col.mesh, add, engine, polynum,
                col.bg, depth, newpage)
    invisible(t(rot.mat))
}

loadRGL <- function() {
    if (! suppressWarnings(require(rgl,quietly=TRUE)))
        stop("rgl is not available")
}

phongLighting <- function(normals, view, light, color, color2, alpha,
                          material = "default") {
    if (length(light) == 4) {
        LI <- light[4]
        light <- light[1:3]
    }
    else LI <- 1
    if (is.character(material))
        material <- getMaterial(material)
    ambient <- material$ambient
    diffuse <- material$diffuse
    specular <- material$specular
    exponent <- material$exponent
    sr <- material$sr
    V <- view / sqrt(sum(view^2))
    L <- light / sqrt(sum(light^2))
    H <- (L + V) / sqrt(sum((L + V)^2))
    sgn <- as.vector(normals %*% V) > 0
    N <- ifelse(sgn,1, -1) * normals
    Is <- as.vector(specular * abs(N %*% H) ^ exponent)
    Id <-  as.vector(diffuse * pmax(N %*% L,0))
    rgbcol <- t(col2rgb(ifelse(sgn, color, color2)) / 255)
    Lrgbcol <- pmin(LI * ((ambient + Id + sr * Is) * rgbcol + (1 - sr) * Is),
                    1)
    Lrgbcol[is.na(Lrgbcol)] <- 0
    rgb(Lrgbcol[,1], Lrgbcol[,2], Lrgbcol[,3], alpha)
}
materials.database <- new.env(hash = TRUE)
registerMaterial <- function(name, ambient = 0.3, diffuse = 0.7, specular = 0.1, exponent = 10, sr = 0) {
    value <- list(ambient = ambient, diffuse = diffuse, specular = specular, exponent = exponent, sr = sr)
    assign(name, value, materials.database)
}
registerMaterial("shiny", ambient = 0.3, diffuse = 0.6, specular = 0.9, exponent = 20, sr = 0)
registerMaterial("dull", ambient = 0.3, diffuse = 0.8, specular = 0.0, exponent = 10, sr = 0)
registerMaterial("metal", ambient = 0.3, diffuse = 0.3, specular = 1.0, exponent = 25, sr = 0.5)
registerMaterial("default", ambient = 0.3, diffuse = 0.7, specular = 0.1, exponent = 10, sr = 0)
registerMaterial("metal", ambient = 0.45, diffuse = 0.45, specular = 1.5, exponent = 25, sr = 0.5)
registerMaterial("shiny", ambient = 0.36, diffuse = 0.72, specular = 1.08, exponent = 20, sr = 0)
getMaterial <- function(name) {
    if (exists(name, materials.database, inherits = FALSE))
        get(name, materials.database)
    else get("default", materials.database, inherits = FALSE)
}

perspLighting <- function(normals, view, light, color, color2, alpha,
                          material = "default") {
    if (length(light) == 4) {
        LI <- light[4]
        light <- light[1:3]
    }
    else LI <- 1
    if (is.character(material))
        material <- getMaterial(material)
    exponent <- material$exponent
    V <- view / sqrt(sum(view^2))
    L <- light / sqrt(sum(light^2))
    sgn <- as.vector(normals %*% V) > 0
    N <- ifelse(sgn,1, -1) * normals
    I <-  (pmax(1 + as.vector(N %*% L), 0) / 2) ^ (exponent * (3 / 40))
    Lrgbcol <- I * LI * t(col2rgb(ifelse(sgn, color, color2)) / 255)
    rgb(Lrgbcol[,1], Lrgbcol[,2], Lrgbcol[,3], alpha)
}

triangleNormalsPhong <- function(triangles) {
    N <- triangleNormals(triangles)
    ve <- t2ve(triangles)
    vt <- vertexTriangles(ve)
    VN <- vertexNormals(vt, N)
    interpolateVertexNormals(VN, ve$ib)
}

triangleNormalsPhongEX <- function(triangles, reps = 1) {
    N <- triangleNormals(triangles)
    ve <- t2ve(triangles)
    vt <- vertexTriangles(ve)
    VN <- vertexNormals(vt, N)
    vb <- ve$vb
    ib <- ve$ib
    n.tri <- nrow(N)
    while (reps > 0) {
	reps <- reps - 1
	n.ver <- nrow(VN)
	mt <- triangleMidTriangles(vb, ib, VN)
	vb <- cbind(vb, mt$vb)
	VN <- rbind(VN, mt$VN)
	mtib <- mt$ib + n.ver
	ib <- matrix(rbind(ib[1,], mtib[1,],mtib[3,],
			   mtib[1,], ib[2,],mtib[2,],
			   mtib[2,], ib[3,], mtib[3,],
			   mtib), nrow = 3)

	for (i in seq(along = triangles))
	    if (length(triangles[[i]]) == n.tri)
		triangles[[i]] <- rep(triangles[[i]], each = 4)
	   n.tri <- 4 * n.tri
    }
    triangles$N <- interpolateVertexNormals(VN, ib)
    triangles$v1 <- t(vb[,ib[1,]])
    triangles$v2 <- t(vb[,ib[2,]])
    triangles$v3 <- t(vb[,ib[3,]])
    triangles
}

triangleNormalsPhongEX <- function(triangles, reps = 1) {
    N <- triangleNormals(triangles)
    ve <- t2ve(triangles)
    vt <- vertexTriangles(ve)
    VN <- vertexNormals(vt, N)
    vb <- ve$vb
    ib <- ve$ib
    n.tri <- nrow(N)
    color <- rep(triangles$color, length = n.tri)
    color2 <- rep(triangles$color2, length = n.tri)
    col.mesh <- rep(triangles$col.mesh, length = n.tri)
    color2 <- ifelse(is.na(color2), color, color2)
    col.mesh <- ifelse(is.na(col.mesh), color, col.mesh)
    VC <- vertexColors(vt, color)
    VC2 <- vertexColors(vt, color2)
    VCm <- vertexColors(vt, col.mesh)

    while (reps > 0) {
        reps <- reps - 1
        n.ver <- nrow(VN)
        edges <- triangleEdges(vb, ib)
        VC <- rbind(VC, (VC[edges[1,],] + VC[edges[2,],]) / 2)
        VC2 <- rbind(VC2, (VC2[edges[1,],] + VC2[edges[2,],]) / 2)
        VCm <- rbind(VCm, (VCm[edges[1,],] + VCm[edges[2,],]) / 2)
        mt <- triangleMidTriangles(vb, ib, VN)
        vb <- cbind(vb, mt$vb)
        VN <- rbind(VN, mt$VN)
        mtib <- mt$ib + n.ver
        ib <- matrix(rbind(ib[1,], mtib[1,],mtib[3,],
                           mtib[1,], ib[2,],mtib[2,],
                           mtib[2,], ib[3,], mtib[3,],
                           mtib), nrow = 3)

        for (i in seq(along = triangles))
            if (length(triangles[[i]]) == n.tri)
                triangles[[i]] <- rep(triangles[[i]], each = 4)
        n.tri <- 4 * n.tri
    }   
    triangles$color <- interpolateVertexColors(VC, ib)
    triangles$color2 <- interpolateVertexColors(VC2, ib)
    triangles$color.mesh <- interpolateVertexColors(VCm, ib)
    triangles$N <- interpolateVertexNormals(VN, ib)
    triangles$v1 <- t(vb[,ib[1,]])
    triangles$v2 <- t(vb[,ib[2,]])
    triangles$v3 <- t(vb[,ib[3,]])
    triangles
}

lightTriangles <- function(triangles, lighting, light) {
    view <- c(0, 0, 1)
    normals <- triangleNormals(triangles)
    smooth <- if (is.null(triangles$smooth)) 0 else triangles$smooth
    if (smooth == 0)
        normals <- triangleNormals(triangles)
    else if (smooth == 1)
        normals <- triangleNormalsPhong(triangles)
    else {
        triangles <- triangleNormalsPhongEX(triangles, reps = smooth - 1)
        normals <- triangles$N
    }
    n.tri <- nrow(normals)
    color <- rep(triangles$color, length = n.tri)
    color2 <- rep(triangles$color2, length = n.tri)
    color2 <- ifelse(is.na(color2), color, color2)
    alpha <- rep(triangles$alpha, length = n.tri)
    mat <- triangles$material
    triangles$col.light <- lighting(normals, view, light, color, color2, alpha, mat)
    triangles
}

lightScene <- function(scene, lighting, light) {
    if (is.Triangles3D(scene))  lightTriangles(scene, lighting, light)
    else lapply(scene, lightTriangles, lighting, light)
}
