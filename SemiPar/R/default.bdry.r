########### S function: default.bdry ############

# For formation of polygons around 
# data points 

# Last changed: 05 JAN 2005

default.bdry <- function(x,y,bdry=NULL)
{ 
   if (!is.null(bdry))
      bdry <- as.matrix(bdry)

   make.circle <- function(curr.point,edgx,edgy,pct=0.95)
   {
      circ.x.vert <- curr.point[1] + pct*edgx*c(1,1/2,-1/2,-1,-1,-1/2,1/2,1,1)
      circ.y.vert <- 
                curr.point[2] + pct*edgy*c(1/2,1,1,1/2,-1/2,-1,-1,-1/2,1/2)
      return(cbind(circ.x.vert,circ.y.vert))
   }

   poly.around.point <- function(poly,point) pointsInPoly(point,poly)


   if (missing(x)&missing(y)&!is.null(bdry))
   {
      x <- bdry[1,1]
      y <- bdry[1,2]
   }

   if (missing(x)&missing(y)&is.null(bdry))
      stop("Need to input either (x,y) or bdry")

   if (is.null(bdry))
   {         
#     get(getOption("device"))()
      par(mfrow=c(1,1))
      plot(x,y,pch=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")

      cat("\n\n")
      cat("        Using the left mouse button click\n")
      cat("        *clockwise* on the points you want\n")
      cat("        to use to form a polygon.\n\n")
      cat("        The polygon should be closed by\n")
      cat("        clicking inside the circle corres- \n")
      cat("        ponding to the starting vertex.\n")
   
      circ.fac <- 3

      finished <- F
      num.points <- 0
      while(!finished)
      {

         loc.out <- locator(1,type="p",pch=4)

         curr.point <- cbind(loc.out$x,loc.out$y)

         num.points <- num.points + 1

         if (num.points==1)
         {
            # Draw circle around starting point

            cex.val <- 0.025
            frame.ratio <- 185/135
            edgx <- cex.val*(max(x)-min(x))
            edgy <- frame.ratio*cex.val*(max(y)-min(y)) 

            circ.x.vert <- 
                curr.point[1] + 0.95*edgx*c(1,1/2,-1/2,-1,-1,-1/2,1/2,1,1)
            circ.y.vert <- 
                curr.point[2] + 0.95*edgy*c(1/2,1,1,1/2,-1/2,-1,-1,-1/2,1/2)
   
            polygon(circ.x.vert,circ.y.vert,density=0)
  
            points(curr.point,pch=18)
         }

         bdry <- rbind(bdry,curr.point)

         if (num.points>1)
            lines(c(curr.point[1],bdry[(num.points-1),1]),
                  c(curr.point[2],bdry[(num.points-1),2]))

         # Test to see if polygon is complete 

         pip.out <- pointsInPoly(curr.point,cbind(circ.x.vert,circ.y.vert))
   
         if ((pip.out==T)&(num.points>1)) finished <- T 
   
      }

      # Set last vertex of polygon to equal the first

      bdry[nrow(bdry),] <- bdry[1,]
   }

   # Plot data and boundary

   x.limits <- range(c(x,bdry[,1]))   
   y.limits <- range(c(y,bdry[,2]))

   plot(x,y,pch=1,xlim=x.limits,ylim=y.limits,xaxt="n",yaxt="n",
        bty="n",xlab="",ylab="")

   polygon(bdry,density=0)

   points(bdry,pch=4)
  
   return(bdry)
}

############## End default.bdry ###############
