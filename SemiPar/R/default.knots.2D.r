###### R-function: default.knots.2D ##########

# Computes default knots for a given
# x vector based on the `clara' algorithm
# in the R package  `cluster'.

# Last changed: 21 OCT 2005

default.knots.2D <- function(x1,x2,num.knots)
{

   # Set default value for num.knots

   if (missing(num.knots))
       num.knots <- max(10,min(50,round(length(x1)/4)))

   # Delete repeated values from x

   X <- cbind(x1,x2)

   dup.inds <- (1:nrow(X))[dup.matrix(X)==T]

   if (length(dup.inds) > 0)
      X <- X[-dup.inds,]

   # Obtain and output knots chosen using
   # coverage design principles

   knots <- cluster::clara(X,num.knots)$medoids

   # Display the knots.

   plot(x1,x2,xlab="",ylab="",bty="n",pch=1,xaxt="n",yaxt="n")
   points(knots[,1],knots[,2],col="red",lwd=3,cex=2)

   cat("\n\n")
   cat("   The knots chosen by default are displayed\n")
   cat("   in the accompanying graphics window (thick red circles),\n")
   cat("   along with the bivariate predictors (ordinary circles).\n")
   cat("   If the knots do not seem to provide adequate coverage\n")
   cat("   then you may want to try selecting them again using \n")
   cat("   the function default.knots.2D() with a different value of\n")
   cat("   the num.knots argument.\n\n")

   cat("   If you are satisfied with these knots then hit return\n")
   cat("   to continue.\n")

   ans <- readline()

   dev.off()

   # Ask user if knots are to be saved in a file.

   cat("\n\n Would you like to save these knots in a file? (y/n): ")

   ans <- readline()

   if (ans=="y")
   {
      cat("\n\n Enter filename for storage of knots: ")

      knots.filename <- readline()

      write(t(knots),file=knots.filename,ncolumns=2)
   }

   return(knots)
}

########## End of default.knots.2D ##########
