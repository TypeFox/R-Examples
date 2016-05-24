########## R-function: lines.spm ##########

# For adding a lines to a plot for a non-additive
# spm() fit.

# Last changed: 01 FEB 2005

lines.spm <- function(x,...)
{
   object <- x

   plot.params <- plotControl(...) 

   # Check whether or not model is a `bona fide' additive model,
   # since centring done only in this case.

   add.model <- add.model.flag(object$info)

   if (add.model==TRUE) stop("lines.spm() not applicable to additive models")

   plot.params <- set.plot.dflts(object,plot.params,1)

   drv <- plot.params$drv
   se <- plot.params$se
   shade <- plot.params$shade
   grid.size <- plot.params$grid.size
   xlim <- plot.params$xlim
   lty <- plot.params$lty
   lwd <- plot.params$lwd
   col <- plot.params$col
   se.lty <- plot.params$se.lty
   se.lwd <- plot.params$se.lwd
   se.col <- plot.params$se.col
   shade.col <- plot.params$shade.col

   # Extract coefficients
 
   X.coefs <- object$fit$coef$fixed
   Z.coefs <- object$fit$coef$random
   coefs <- c(X.coefs,Z.coefs)

   # Convert 'xlim' specifications to list if necessary.

   if ((!is.null(plot.params$xlim))&(!is.list(plot.params$xlim)))
   {
      if (length(plot.params$xlim)!=2) stop("illegal xlim")
      plot.params$xlim <- list(lower=rep(plot.params$xlim[1],1),
                               upper=rep(plot.params$xlim[2],1))
   }

   # Plot curve estimates

   lin.present <- !is.null(object$info$lin)
   pen.present <- !is.null(object$info$pen)

   if (lin.present|pen.present)
   {

      # Obtain plots for penalized components,
      # starting with the x-values.

      x.vals <- NULL
      if (lin.present)
         x.vals <- cbind(x.vals,object$info$lin$x)
      if (pen.present)
         x.vals <- cbind(x.vals,object$info$pen$x)

      if (lin.present)
         curve.type <- "lin"

      if (pen.present)
         curve.type <- "pen"

      num.curves <- length(curve.type)

      # Set default plotting parameters 

      plot.params <- set.plot.dflts(object,plot.params,num.curves)
    
      pen.num <- 1

      # Extract block indices

      block.inds <- object$aux$block.inds

      # Set up x grid 

      x.grid <- seq(plot.params$xlim$lower[1],
                    plot.params$xlim$upper[1],
                    length=plot.params$grid.size[1])
     
      if (curve.type=="lin") 
      {
         deg.val <- 1
         ncol.X.val <- 2
      }

      if (curve.type=="pen") 
      {
         deg.val <- object$info$pen$degree[pen.num]
           
         # Extract type of basis

         basis <- object$info$pen$basis

         if(basis=="trunc.poly")
            ncol.X.val <- 1 + deg.val
         if(basis=="tps")
            ncol.X.val <- 1 + (deg.val-1)/2
      }  

      X.grid <- rep(1,plot.params$grid.size[1])

      for (j in 2:ncol.X.val)
        X.grid <- cbind(X.grid,x.grid^(j-1))

      y.grid <- as.matrix(X.grid)%*%X.coefs

      if (curve.type=="pen")
      {
         # Set knot values

         knots <- object$info$pen$knots[[1]]

         Z.grid <- outer(x.grid,knots,"-")

         if (basis=="trunc.poly")
            Z.grid <- (Z.grid*(Z.grid>0))^deg.val
   
         if (basis=="tps")
         {
           Z.grid <- abs(Z.grid)^deg.val
           sqrt.Omega <-  object$info$trans.mat[[pen.num]]
           Z.grid <- t(solve(sqrt.Omega,t(Z.grid)))
         }                               
          
         y.grid <- y.grid + as.matrix(Z.grid)%*%Z.coefs
       }

       # Compute standard error bars if requested.

       if (se==TRUE)
       {  
          cov.mat <- object$aux$cov.mat

          curr.inds <- c(1,block.inds[[2]])

          C.grid <- X.grid
   
          if (curve.type=="pen")
             C.grid <- cbind(C.grid,Z.grid)
 
          cov.mat.curr <- as.matrix(cov.mat[curr.inds,curr.inds])

          se.grid <- sqrt(diag(C.grid%*%cov.mat.curr%*%t(C.grid)))

          lower <- y.grid - 2*se.grid
          upper <- y.grid + 2*se.grid

       
          if (shade==FALSE) 
          {
             lines(x.grid,lower,col=se.col,lty=se.lty,lwd=se.lwd)
             lines(x.grid,upper,col=se.col,lty=se.lty,lwd=se.lwd)
          }

          if (shade==TRUE)
          {
             # Make sure points are sorted according to their
             # abcissae         

             mat <- cbind(x.grid,lower,upper)
             mat <- mat[sort.list(mat[,1]),]
 
             x.grid <- mat[,1]
             lower <- mat[,2]
             upper <- mat[,3]
 
             polygon(c(x.grid,rev(x.grid)),c(lower,rev(upper)),
                     col=shade.col,border=FALSE)
  
         }
      }

      lines(x.grid,y.grid,lty=lty,lwd=lwd,col=col)
   }

   invisible()
}

########## End of R-function: lines.spm ##########
