########## R function: plot.curvHDRfilter ###########
 
# For visualisation of a curvHDRfilter object.

# Last changed: 16 NOV 2010

plot.curvHDRfilter <- function(x,removeDebri=TRUE,pch=NULL,cex=NULL,
                               bty=NULL,col=NULL,main=NULL,...)
{
   curvHDRfilterObj <- x

   # Extract relevant information from input object:

   data <- curvHDRfilterObj$data
   insideFilter <- curvHDRfilterObj$insideFilter
   polys <- curvHDRfilterObj$polys
   HDRlevel <- curvHDRfilterObj$HDRlevel
   dmn <- ncol(data)
   
   # Remove debri if flag is set:

   if (removeDebri)
   {
      mins <- apply(data,2,min)   
      maxs <- apply(data,2,max)
      omitInds <- NULL
      for (j in 1:dmn)
         omitInds <- c(omitInds,(1:nrow(data))[(data[,j]==mins[j])|(data[,j]==maxs[j])])
      if (!is.null(omitInds))
      {
         data <- as.matrix(data[-omitInds,])
         insideFilter <- insideFilter[-omitInds]
      }      
   }

   # Set plotting defaults:

   if (is.null(pch)) pch <- 1 ;   if (is.null(cex)) cex <- 0.1
   if (is.null(bty)) bty <- "l" 
   if (is.null(main)) main <- paste("curvHDR filter with HDR level=",HDRlevel)
   if (is.null(col))
   {
      if (dmn==1) col <- "orange"
      if (dmn>1)
      {
         col <- rep("orange",nrow(data))
         col[insideFilter] <- "red"
      }   
   }

   # Plot data and curvHDR gates:
   
   if (ncol(data)==1)    # Univariate data.
   {
      hist(data,freq=FALSE,breaks=100,bty=bty,
           col=col,main=main,cex.main=1.5,...)
      histObj <- hist(data,breaks=100,plot=FALSE)
      histMax <- max(histObj$density)
      stripWidth <- histMax/20
      addIntervs(polys,stripWidth)
   }
   
   if (ncol(data)==2)    # Bivariate data.
   {
      plot(data[,1],data[,2],pch=pch,cex=cex,bty=bty,
           col=col,main=main,cex.main=1.5,...)

      for (j in 1:length(polys))
         polygon(polys[[j]],border="blue",lwd=4)
   }

   if (ncol(data)==3)    # Trivariate data.
   {
      includePartial <- TRUE  # Note: this flag is to allow possible
                              # future fine-tuning that allows for
                              # polyhedra to be excluded if not fully
                              # inside plot limits.
      if (nrow(data)>1000)
      {
         subDataInds <- sample(1:nrow(data),1000)
         subData <- data[subDataInds,]
         subCol <- col[subDataInds]
      }
      if (nrow(data)<=1000)
      {
         subData <- data
         subCol <- col
      }
         
      pobj <- plot3D(subData[,1],subData[,2],subData[,3],
                     main=main,ptCol=subCol,ptAlpha=0.2,cex=0.5,...)

      avec <- c(attr(pobj,"a.x"),attr(pobj,"a.y"),attr(pobj,"a.z"))
      bvec <- c(attr(pobj,"b.x"),attr(pobj,"b.y"),attr(pobj,"b.z"))

      # Transform polyhedra to same space as plot3D():

      tpolys <- polys

      for (j in 1:length(polys))
      {
         for (k in 1:3)
         {
            tpolys[[j]]$v1[,k] <-  tranUnitInt(polys[[j]]$v1[,k],avec[k],bvec[k])
            tpolys[[j]]$v2[,k] <-  tranUnitInt(polys[[j]]$v2[,k],avec[k],bvec[k])
            tpolys[[j]]$v3[,k] <-  tranUnitInt(polys[[j]]$v3[,k],avec[k],bvec[k])
         }

         # Determine whether current polyhedron is
         # (a) fully within plot limits,
         # (b) partially within plot limits, or
         # (c) completely outside plot limits.
             
         intVec <- rep(NA,3)
         for (k in 1:3)
         {
            currVerts <- c(tpolys[[j]]$v1[,k],tpolys[[j]]$v2[,k],tpolys[[j]]$v3[,k])
            currMin <- min(currVerts) ; currMax <- max(currVerts)
            intVec[k] <- "partial"
            if (currMax<0) intVec[k] <- "outside"
            if (currMin>1) intVec[k] <- "outside"
            if ((currMin>=0)&(currMax<=1)) intVec[k] <- "inside"
         }

         type <- "partiallyInside"
         if (any(intVec=="outside")) type <- "fullyOutside"
         if (all(intVec=="inside")) type <- "fullyInside"

         # Continue only for polyhedra with some part within plot limits:

         if ((type=="fullyInside")|(type=="partiallyInside"))
         {
            # Plot those polyhedra fully within the plot limits:

            if (type=="fullyInside")              
               drawScene.rgl(tpolys[[j]],add=TRUE)         

            # Plot those polyhedra partially within the plot limits
            # if `includePartial' is true:
            
            if (type=="partiallyInside")
            {
               if (includePartial)
               {
                  # Truncate polyhedra with parts beyond the plot limits:
         
                  for (k in 1:3)
                  {
                     for (kk in 1:3)
                     {
                        currVerts <- tpolys[[j]][[kk]][,k]
                        currVerts[currVerts<0.02] <- 0.02
                        currVerts[currVerts>0.98] <- 0.98
                        tpolys[[j]][[kk]][,k] <- currVerts
                     }
                  }
                  drawScene.rgl(tpolys[[j]],add=TRUE) 
               }
            }
         }
      }
   }
   invisible()
}

############ End of plot.curvHDRfilter ###########
