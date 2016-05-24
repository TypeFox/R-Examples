########## R function: curvHDRfilter ###########
 
# For obtaining a curvHDR filter from a set of
# multivariate continuous data. The current version
# supports univaraite, bivariate and trivariate data.

# Last changed: 07 MAR 2016

#' @importFrom rgl rgl.clear rgl.bg rgl.spheres
#' @importFrom misc3d contour3d
#' @importFrom graphics polygon plot hist
#' @importFrom stats IQR sd qnorm
#' @importFrom feature featureSignif
#' @importFrom grDevices contourLines
#' @importFrom misc3d contour3d updateTriangles drawScene.rgl
#' @importFrom geometry convhulln
#' @importFrom ptinpoly pip3d
#' @importFrom ks hpi Hpi.diag kde contourLevels
#' @importFrom hdrcde hdr
#' @importFrom KernSmooth dpik

curvHDRfilter <- function(x,HDRlevel=0.1,growthFac=NULL,signifLevel=0.05,bwFac=1,
                          gridsize=NULL,removeDebri=TRUE,minSampSize=NULL,
                          HpiGridSize=NULL,quiet=TRUE,graphChk=FALSE)
{

   # Make sure input data is a matrix with between 1 and 3 columns:

   x <- as.matrix(x)
   dmn <- ncol(x)
   xOrig <- x

   if (!any(ncol(x)==(1:3)))
      stop("x should be a vector or a matrix with between 1 and 3 columns.")
   
   # Set various defaults:

   if (is.null(growthFac))
      growthFac <- 2^dmn
  
   if (is.null(gridsize))
   {
     if (dmn==1) gridsize <- 801
     if (dmn==2) gridsize <- rep(151,2)
     if (dmn==3) gridsize <- rep(51,3)
   } 
   
   if ((dmn==3)&(is.null(HpiGridSize)))
      HpiGridSize <- rep(21,3)
  
   if (is.null(minSampSize)) minSampSize <- 50*(2^(dmn-1))

   # Remove debri if flag is set:

   if (removeDebri)
   {
      mins <- apply(x,2,min)   
      maxs <- apply(x,2,max)
      omitInds <- NULL
      for (j in 1:dmn)
         omitInds <- c(omitInds,(1:nrow(x))[(x[,j]==mins[j])|(x[,j]==maxs[j])])
      if (!is.null(omitInds))
         x <- as.matrix(x[-omitInds,])
   }

   # Perform robust standardisation of data:

   stDevs <- apply(x,2,sd)
   IQRvals <- apply(x,2,IQR)
   corrFac <- qnorm(3/4) - qnorm(1/4)
   sigHatsx <- apply(cbind(stDevs,IQRvals/corrFac),1,min)
   meanx <- apply(x,2,mean)
   for (j in 1:dmn)
      x[,j] <- (x[,j]-meanx[j])/sigHatsx[j]   
      
   # Obtain normal scale bandwidths:

   sampSizeFac <- (4/((dmn+6)*nrow(x)))^(1/(dmn+8))
   bwNS <- rep(sampSizeFac,dmn)

   # Set up plot if `graphChk' is TRUE:

   if (graphChk)
   {
      if (dmn==1)
      {
         histObj <- hist(x,freq=FALSE,breaks=100,col="orange",bty="l",
                         xlab="standardised data",main="graphical check in progress",
                         cex.main=1.7)
         histMax <- max(histObj$density)
         stripWidth <- histMax/40
         curPolLow <- 0 ; curPolUpp <- stripWidth
      }
      if (dmn==2)
        plot(x[,1],x[,2],col="orange",bty="l",xlab="standardised first variable",pch=".",
              ylab="standardised first variable",main="graphical check in progress",
              cex.main=1.7)

      if (dmn==3)
      {
           if (nrow(x)<=2000) xSub <- x
           if (nrow(x)> 2000) xSub <- x[sample(1:nrow(x),2000),]
           rgl.clear() ; rgl.bg(color="white")
           rgl.spheres(xSub[,1],xSub[,2],xSub[,3],col="orange",alpha=0.2,radius=0.025)
      }
        
      message("Hit Enter to continue.") ; ans <- readline()
   }

   if (!quiet)
   {
      ## cat("\n")
      message("Computing high curvature regions ...")
   }
   
   # Obtain high curvature regions:

   fSobj <- featureSignif(x,bw=bwFac*bwNS,gridsize=gridsize,signifLevel=signifLevel)
   
   if (dmn==1)
   {
      x.grid <- fSobj$fhat$x.grid[[1]]
      y.grid <- sign(as.numeric(fSobj$curv))
      curvPolys <- contour1D(0,x.grid,y.grid)
   }
   
   if (dmn==2)
       curvPolys <- contourLines(fSobj$fhat$x.grid[[1]],fSobj$fhat$x.grid[[2]],
                                 fSobj$curv,levels=0.5)
       
   if (dmn==3)
   {    
      contour3dObjCurv <- contour3d(fSobj$curv,0.5,fSobj$fhat$x.grid[[1]],
                                fSobj$fhat$x.grid[[2]],
                                fSobj$fhat$x.grid[[3]],engine="none")

      # Separate into contiguous components:

      curvPolys <- separateTriPolyh(contour3dObjCurv)
   }

   numCurvPolys <- length(curvPolys)

   # Show significant curvature polygons if `graphChk' is TRUE:

   if (graphChk)
   {
      for (j in 1:numCurvPolys)
      {
         if (dmn==1) 
         {
            curPoly <- list(x=c(curvPolys[[j]],rev(curvPolys[[j]])),
                            y=c(rep(curPolLow,2),rep(curPolUpp,2)))
            polygon(curPoly,col="darkmagenta")
            curPolLow <- curPolLow + stripWidth ; curPolUpp <- curPolUpp + stripWidth
         }
           
         if (dmn==2) 
            polygon(curvPolys[[j]],border="darkmagenta",lwd=2)
        
         if (dmn==3)
         {      
            curvPolys[[j]] <- updateTriangles(curvPolys[[j]],color="darkmagenta",alpha=0.15)
            drawScene.rgl(curvPolys[[j]],add=TRUE)   
         }          
      }
      message("Hit Enter to continue.") ; ans <- readline()
   }
   
   if (!quiet)
      message("Growing convexified high curvature regions ...")

   # Convexify polygons/polyhedra (if d>1).
   # Then grow them by specified growth factor
   # using  bisection searches:

   grownPolys <- list()
   for (j in 1:numCurvPolys)
   {
      currCurvPoly <- curvPolys[[j]]
           
      if (dmn==1)
      {
         currInitPoly <- currCurvPoly
         Vinit <- abs(currInitPoly[2]-currInitPoly[1])
      }

      if (dmn==2)
      {
         currCurvPolyMat <- cbind(currCurvPoly$x,currCurvPoly$y)
         hpts <- chull(currCurvPolyMat[!is.na(rowSums(currCurvPolyMat)),])	
	 hpts <- c(hpts,hpts[1]) 	
	 convHullPolyMat <- cbind(currCurvPolyMat[hpts,1],currCurvPolyMat[hpts,2])
         currChull <- list(x=convHullPolyMat[,1],y=convHullPolyMat[,2])
         currInitPoly <- currChull
         Vinit <- areaPolygon(currInitPoly)
      }
      
      if (dmn==3)
      {
         uniqVertices <- unique(rbind(currCurvPoly$v1,currCurvPoly$v2,
                                   currCurvPoly$v3))
         currChull <- convhulln(uniqVertices,options="FA")
         currChull$v1 <- uniqVertices[currChull$hull[,1],]
         currChull$v2 <- uniqVertices[currChull$hull[,2],]
         currChull$v3 <- uniqVertices[currChull$hull[,3],]
         currInitPoly <- currChull
         Vinit <- currInitPoly$vol
      }

      if (dmn==1) rUppFac <- 1
      if (dmn==2) rUppFac <- 1
      if (dmn==3) rUppFac <- 3/4
      rLow <- 0 ; rUpp <- (rUppFac*Vinit/pi)^(1/dmn)*(growthFac^(1/dmn)-1)

      # Make sure that root is captured:

      rootCaptured <- FALSE
      while (!rootCaptured)
      {
         if (dmn==1) Vupp <- intervalGrow(currInitPoly,rUpp)$volume
         if (dmn==2) Vupp <- polygonGrow(currInitPoly,rUpp)$volume
         if (dmn==3) Vupp <- polyhedronGrow(currInitPoly,rUpp)$volume
         if (Vupp>(growthFac*Vinit))
            rootCaptured <- TRUE
         else
            rUpp <- 1.2*rUpp
      }

      tolerance <- 0.0001
      relErr <- tolerance + 1
      while (relErr>tolerance)
      {
         rMid <- (rLow + rUpp)/2

         if (dmn==1)
         {
            newVolLow <- intervalGrow(currInitPoly,rLow)$volume
            newVolMid <- intervalGrow(currInitPoly,rMid)$volume
            newVolUpp <- intervalGrow(currInitPoly,rUpp)$volume
         }
         
         if (dmn==2)
         {
            newVolLow <- polygonGrow(currInitPoly,rLow)$volume
            newVolMid <- polygonGrow(currInitPoly,rMid)$volume
            newVolUpp <- polygonGrow(currInitPoly,rUpp)$volume
         }
         
         if (dmn==3)
         {
            newVolLow <- polyhedronGrow(currInitPoly,rLow)$volume
            newVolMid <- polyhedronGrow(currInitPoly,rMid)$volume
            newVolUpp <- polyhedronGrow(currInitPoly,rUpp)$volume
         }
           
         fLow <- newVolLow - growthFac*Vinit
         fUpp <- newVolUpp - growthFac*Vinit
         fMid <- newVolMid - growthFac*Vinit
         if (fMid<0) rLow <- rMid
         if (fMid>0) rUpp <- rMid
         relErr <- max(c((rUpp-rLow)/rMid,(fUpp-fLow)/fMid))
      }

      if (dmn==1)
         grownPolys[[j]] <- intervalGrow(currInitPoly,rMid)$grownInterval
      if (dmn==2)
         grownPolys[[j]] <- polygonGrow(currInitPoly,rMid)$grownPolygon
      if (dmn==3)
         grownPolys[[j]] <- polyhedronGrow(currInitPoly,rMid)$grownPolyHedron
   }

   # Add grown polyhedra if graphChk is TRUE:
      
   if (graphChk)
   {
      for (j in 1:length(grownPolys))
      {
         if (dmn==1)
         {
            curPoly <- list(x=c(grownPolys[[j]],rev(grownPolys[[j]])),
                            y=c(rep(curPolLow,2),rep(curPolUpp,2)))
            polygon(curPoly,col="forestgreen")
            curPolLow <- curPolLow + stripWidth ; curPolUpp <- curPolUpp + stripWidth
         }
         if (dmn==2)
            polygon(grownPolys[[j]],border="forestgreen",lwd=2)
         if (dmn==3)
         {
            grownPolys[[j]] <- updateTriangles(grownPolys[[j]],color="forestgreen",alpha=0.15)
            drawScene.rgl(grownPolys[[j]],add=TRUE)
         }      
      }  
      message("Hit Enter to continue.") ; ans <- readline()
   }

   if (!quiet)
         message("Obtaining highest density regions ...")

   # Obtain highest density regions based on data within each
   # grown polyhedron:

   curvHDRpolys <- list()
   jHDR <- 1 
   for (jGrow in 1:length(grownPolys))
   {
      if (dmn==1) interiorIndics <- ((x>=grownPolys[[jGrow]][1])&(x<=grownPolys[[jGrow]][2]))
      
      if (dmn==2)
      {  
         grownPolyMat <- cbind(grownPolys[[jGrow]]$x,grownPolys[[jGrow]]$y)
         interiorIndics <- as.logical(inpolygon(x,grownPolyMat))
      }
         
      if (dmn==3)
      {
         pip3dInfo <- t2ve(grownPolys[[jGrow]])
         currVertices <- t(pip3dInfo$vb)
         currFaces    <- t(pip3dInfo$ib)
         interiorIndics <- as.logical(pip3d(currVertices,currFaces,x)==1)
      }
        
      currSampSize <- sum(as.numeric(interiorIndics))

      if (currSampSize>=minSampSize)
      {
         kdeData  <- x[interiorIndics==TRUE,]

         # Obtain bandwidth(s):
         
         if (dmn==1) hVal <- hpi(kdeData)
         if (dmn==2) Hmat <- Hpi.diag(kdeData,binned=TRUE)
         if (dmn==3) Hmat <- Hpi.diag(kdeData,pilot="samse",binned=TRUE,
                                      bgridsize=HpiGridSize)
         
         # Obtain kernel density estimator:

         if (dmn>1)
            kDest <- kde(kdeData,H=Hmat,binned=TRUE)

         # Obtain HDR polygons:
         
         if (dmn==1)
         {
            hdr.obj <- hdr(kdeData,prob=100*(1-HDRlevel),h=dpik(kdeData))
            numReg <- length(hdr.obj$hdr)/2
            HDRobj <- NULL
            for (j in 1:numReg)
               HDRobj[[j]] <- hdr.obj$hdr[(2*j-1):(2*j)]
         }

         if (dmn>1)
            hdrLevel <- contourLevels(kDest,prob=HDRlevel)

         if (dmn==2)
            HDRobj <- contourLines(x=kDest$eval.points[[1]],y=kDest$eval.points[[2]],
                                   z=kDest$estimate,levels=hdrLevel)
         

         if (dmn==3)
         {  
            contour3dObjHDR <-  contour3d(kDest$estimate,level=hdrLevel, 
                                          kDest$eval.points[[1]],
                                          kDest$eval.points[[2]], 
                                          kDest$eval.points[[3]],  
                                          engine="none")
            HDRobj <- separateTriPolyh(contour3dObjHDR)
         }

         for (jInn in 1:length(HDRobj))
         {
            curvHDRpolys[[jHDR]] <- HDRobj[[jInn]]
            jHDR <- jHDR + 1
         }
      }
   }      

   # In the bivariate case, remove the "$level" component of
   # the polygons:

   if (dmn==2)
      for (j in 1:length(HDRobj))
         HDRobj[[j]] <- list(x=HDRobj[[j]]$x,y=HDRobj[[j]]$y)
     
   # Add HDR polyhedra if graphChk is TRUE:
   
   if (graphChk)
   {
      for (j in 1:length(curvHDRpolys))
      {
         if (dmn==1)
         {  
            curPoly <- list(x=c(curvHDRpolys[[j]],rev(curvHDRpolys[[j]])),
                            y=c(rep(curPolLow,2),rep(curPolUpp,2)))
            polygon(curPoly,col="blue")
            curPolLow <- curPolLow + stripWidth ; curPolUpp <- curPolUpp + stripWidth
         }
         if (dmn==2)
            polygon(curvHDRpolys[[j]],border="blue",lwd=2)

         if (dmn==3)
         {
            curvHDRpolys[[j]] <- updateTriangles(curvHDRpolys[[j]],color="blue",alpha=0.15)
            drawScene.rgl(curvHDRpolys[[j]],add=TRUE)
         }
      }  
      message("Hit Enter to continue.") ; ans <- readline()   
   }

   for (j in 1:length(curvHDRpolys))
   {
      # Transform to original scale:

      if (dmn==1)
         curvHDRpolys[[j]] <- sigHatsx*curvHDRpolys[[j]] + meanx
      
      if (dmn==2)
      {
         curvHDRpolys[[j]]$x <- sigHatsx[1]*curvHDRpolys[[j]]$x + meanx[1]
         curvHDRpolys[[j]]$y <- sigHatsx[2]*curvHDRpolys[[j]]$y + meanx[2]
      }
         
      if (dmn==3)
      {
         for (k in 1:3)
         {
            curvHDRpolys[[j]]$v1[,k] <- sigHatsx[k]*curvHDRpolys[[j]]$v1[,k] + meanx[k]
            curvHDRpolys[[j]]$v2[,k] <- sigHatsx[k]*curvHDRpolys[[j]]$v2[,k] + meanx[k]
            curvHDRpolys[[j]]$v3[,k] <- sigHatsx[k]*curvHDRpolys[[j]]$v3[,k] + meanx[k]
         }
      }

      # Adjust the colour and alpha values of output polyhedra (trivariate case only):

      if (dmn==3)
         curvHDRpolys[[j]] <- updateTriangles(curvHDRpolys[[j]],color="blue",alpha=0.3)
   }

   # Obtain logical vector named `insideFilter' that indicates which of the
   # input data points are inside at least one of the intervals/polygons/polyhedra:

   inFilterNum <- rep(0,nrow(xOrig))
   for (j in 1:length(curvHDRpolys))
   {
      if (dmn==1)
         inFilterNum <- inFilterNum + as.numeric((xOrig>=curvHDRpolys[[j]][1])&
                                                 (xOrig<=curvHDRpolys[[j]][2]))   
      if (dmn==2)
      {
         currPolyMat <- cbind(curvHDRpolys[[j]]$x,curvHDRpolys[[j]]$y)
         inFilterNum <- inFilterNum + as.numeric(inpolygon(xOrig,currPolyMat))
      }

      if (dmn==3)
      {
         pip3dInfo <- t2ve(curvHDRpolys[[j]])
         currVertices <- t(pip3dInfo$vb)
         currFaces    <- t(pip3dInfo$ib)
         inFilterNum <- inFilterNum + as.numeric(pip3d(currVertices,currFaces,xOrig)==1)
      }
   }
   
   insideFilter <- as.logical(inFilterNum>0)
   
   # Prepare return objects:
   
   outObj <- list(data=xOrig,insideFilter=insideFilter,polys=curvHDRpolys,
                  HDRlevel=HDRlevel)
   
   class(outObj) <- "curvHDRfilter"

   if (!quiet)
     message("Finished obtaining curvHDR filters.")
   
   return(outObj)
}

############ End of curvHDRfilter ###########
