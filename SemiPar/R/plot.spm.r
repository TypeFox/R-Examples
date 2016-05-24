########## R-function: plot.spm ##########

# For plotting results of spm()

# Last changed: 21 JAN 2005

plot.spm <- function(x,...)
{
   object <- x

   plot.params <- plotControl(...)  

   drv <- plot.params$drv
   se <- plot.params$se
   shade <- plot.params$shade

   # Work out if default `ylim' values need to be computed.

   default.ylim <- FALSE
   if (is.null(plot.params$ylim))
      default.ylim <- TRUE

   # Set flags for presence of components of 
   # various types

   lin.present <- !is.null(object$info$lin)
   pen.present <- !is.null(object$info$pen)
   krige.present <- !is.null(object$info$krige)
   random.present <- !is.null(object$info$random)
   const.only <- ((!lin.present)&(!pen.present)&
                  (!krige.present)&(!random.present))

   # Extract block indices

   block.inds <- object$aux$block.inds 

   # Extract coefficients

   if (pen.present|krige.present)
     coefs <- c(object$fit$coef$fixed,object$fit$coef$random)
   else
     coefs <- object$fit$coef$fixed

   # Obtain relevant information for standard errors

   if (se==TRUE)
      cov.mat <- object$aux$cov.mat

   # Remove random intercept coefficients and covariance
   # matrix component if present

   if (random.present)
   {
      random.inds <- block.inds[[length(block.inds)]]

      coefs <- coefs[-random.inds]

      if (se==TRUE)
         cov.mat <- cov.mat[-random.inds,-random.inds]
   }

   # Extract type of basis

   basis <- object$info$pen$basis

   # Set up first row for `underlying' design matrix of 
   # average values

   if (drv==0)
      ave.vals <- 1    # For intercept
   if (drv>0)
      ave.vals <- 0    # For intercept
  
   if (lin.present) # Add on entries for "lin" components
      ave.vals <- c(ave.vals,apply(object$info$lin$x,2,mean))

   if (pen.present) # Add on fixed effect entries for "pen" components
   {
      num.pen <- ncol(as.matrix(object$info$pen$x))
      for (ipen in 1:num.pen)
      {   
         deg.val <- object$info$pen$degree[ipen]
                
         if(basis=="trunc.poly")
            ncol.X.val <- deg.val
         if(basis=="tps")
            ncol.X.val <- (deg.val-1)/2
  
         ave.val <- mean(as.matrix(object$info$pen$x)[,ipen])
         ave.vals <- c(ave.vals,ave.val)
  
         if (ncol.X.val>1)
            for (ipow in 2:ncol.X.val)
               ave.vals <- c(ave.vals,ave.val^ipow)

      }
   }

   if (!pen.present)
      num.pen <- 0

   if (krige.present) # Add on fixed effect entries for "krige" component
   {
      x1.v <- object$info$krige$x[,1]
      x2.v <- object$info$krige$x[,2]

      m.krige <- (object$info$krige$degree/2) + 1

      X.kg <- NULL
      for (im in 2:m.krige)
      {
         X.kg <- cbind(X.kg,x1.v^(im-1))
         if (im>2)
            for (ipow in 2:(im-1))
               X.kg <- cbind(X.kg,(x1.v^(im-ipow))*(x2.v^(ipow-1)))
         X.kg <- cbind(X.kg,x2.v^(im-1))
      }
     
      ave.vals <- c(ave.vals,apply(X.kg,2,mean))
   }

   if (pen.present) # Add on random effect entries for "pen" components
   {
      num.pen <- ncol(as.matrix(object$info$pen$x))
      for (ipen in 1:num.pen)
      {
         deg.val <- object$info$pen$degree[ipen]
         knots <- object$info$pen$knots[[ipen]]

         ave.val <- mean(as.matrix(object$info$pen$x)[,ipen])
         if(basis=="trunc.poly")
         {
            ave.val.knots <- ave.val - knots
            ave.val.knots <- (ave.val.knots*(ave.val.knots>0))^deg.val
         }
         if (basis=="tps")
         {  
            ave.val.knots <- abs(ave.val - knots)^deg.val
            sqrt.Omega <-  object$info$trans.mat[[ipen]]
            ave.val.knots <- t(solve(sqrt.Omega,ave.val.knots)) 
         }
         ave.vals <- c(ave.vals,ave.val.knots)
      }
   }

   if (krige.present) # Add on random effect entries for "krige" component
   {
       ave.val.x1 <- mean(object$info$krige$x[,1])
       ave.val.x2 <- mean(object$info$krige$x[,2])

       knots <- object$info$krige$knots
       knots.1 <- knots[,1]
       knots.2 <- knots[,2]

       num.knots <- length(knots.1)

       ave.val.knots <- NULL

       dists.ave <- sqrt((ave.val.x1-knots.1)^2
                                +(ave.val.x2-knots.2)^2)         
       ave.val.knots <- tps.cov(dists.ave,m=m.krige,d=2)

       ave.val.knots <- solve(object$info$trans.mat[[length(
                        object$info$trans.mat)]],as.matrix(ave.val.knots))
       ave.vals <- c(ave.vals,ave.val.knots) 
   }
 
   # Generate the plots

   if (const.only)
   {   
      mean.est <- coefs[1]
      x.grid <- c(0.25,0.75)
      se.val <- sqrt(diag(object$aux$cov.mat))[1]
      lower <- mean.est - 2*se.val
      upper <- mean.est + 2*se.val

      plot(0,0,type="n",xlim=c(0,1),ylim=c(mean.est-2.5*se.val,
           mean.est+2.5*se.val),xaxt="n",xlab="",ylab="",bty="l")
        
      if (se==TRUE)
      {
         if (shade==FALSE) 
         {
            lines(x.grid,rep(lower,2),lty=2,lwd=2,col="black")
                  
            lines(x.grid,rep(upper,2),lty=2,lwd=2,col="black")
         }

         if (shade==TRUE) 
            polygon(c(x.grid,rev(x.grid)),c(rep(lower,2),rep(upper,2)),
                    border=FALSE,col="grey70")
      }

      lines(x.grid,rep(mean.est,2),lwd=2)
   }   

   # Set plotting defaults

   if (lin.present|pen.present)
   {  

      # Determine number of and types of curve estimates to be 
      # plotted

      curve.type <- NULL
      if (lin.present)
          curve.type <- c(curve.type,rep("lin",
                          ncol(as.matrix(object$info$lin$x))))

      if (pen.present)
          curve.type <- c(curve.type,rep("pen",
                          ncol(as.matrix(object$info$pen$x))))

      num.curves <- length(curve.type)
   }

   if (!(lin.present|pen.present)) 
      num.curves <- 0

   plot.params <- set.plot.dflts(object,plot.params,num.curves)

   # Convert 'xlim' and 'ylim' specifications to lists if necessary.

   if ((!is.null(plot.params$xlim))&(!is.list(plot.params$xlim)))
   {
      if (length(plot.params$xlim)!=2) stop("illegal xlim")
      plot.params$xlim <- list(lower=rep(plot.params$xlim[1],num.curves),
                               upper=rep(plot.params$xlim[2],num.curves))
   }

   if ((!is.null(plot.params$ylim))&(!default.ylim)
      &(!is.list(plot.params$ylim)))
   {
      if (length(plot.params$ylim)!=2) stop("illegal ylim")
      plot.params$ylim <- list(lower=rep(plot.params$ylim[1],num.curves),
                               upper=rep(plot.params$ylim[2],num.curves))
   }

   # Do the plots

   if (lin.present|pen.present)
   {   
      # Obtain plots for penalised components,
      # starting with the x-values.

      x.vals <- NULL

      if (lin.present)
         x.vals <- cbind(x.vals,object$info$lin$x)
      if (pen.present)
         x.vals <- cbind(x.vals,object$info$pen$x)

      # Set counter for fixed coefficients

      fc.stt.pos <- 2
      pen.num <- 1
      for (j in 1:num.curves)
      {
          # Determine grid size and form `underyling'
          # average matrix

          grid.size <- plot.params$grid.size[j]
  
          if(drv==0)
             C.grid <-  matrix(rep(ave.vals,grid.size),grid.size,
                                   length(ave.vals),byrow=TRUE)

          if (drv>0)
             C.grid <- matrix(0,grid.size,length(ave.vals))

          # Set up grid for jth component

          x.grid <- seq(plot.params$xlim$lower[j],
                        plot.params$xlim$upper[j],
                        length=plot.params$grid.size[j])

          if (curve.type[j]=="lin") 
             deg.val <- 1
          if (curve.type[j]=="pen") 
             deg.val <- object$info$pen$degree[pen.num]

          # Set up insertion for C.grid matrix for current component,

          X.g.inst <- NULL

          if (curve.type[j]=="lin") 
          {
             if (drv==0)
                X.g.inst <- cbind(X.g.inst,x.grid)

             if (drv==1)
                X.g.inst <- cbind(X.g.inst,rep(1,grid.size[j]))

             if (drv>1)
                X.g.inst <- cbind(X.g.inst,rep(0,grid.size[j]))
          }

          if (curve.type[j]=="pen")
          {            
             if(basis=="trunc.poly")
                ncol.X.val <- deg.val
             if(basis=="tps")
                ncol.X.val <- (deg.val-1)/2

             if (drv==0)
             {
                for (pow in 1:ncol.X.val)
                { 
                   new.col.data <- x.vals[,j]^pow
                   new.col <- x.grid^pow
                   X.g.inst <- cbind(X.g.inst,new.col)
                } 
             }

             if (drv>0)
             {
                for (pow in 1:ncol.X.val)
                { 
                   new.col.data <- x.vals[,j]^pow 
                   pow.drv <- pow - drv
                   if (pow.drv>=0)
                      new.col <- prod(pow:(pow.drv+1))*(x.grid^pow.drv)
                   else
                      new.col <- rep(0,length(x.grid))
                   X.g.inst <- cbind(X.g.inst,new.col)
                }
             }

         }

         C.g.inst <- X.g.inst

         # Set up the column indices for the insertion

         if (curve.type[j]=="lin")
         {
            inst.col.inds <- fc.stt.pos
            fc.stt.pos <- fc.stt.pos + 1
         }

         if (curve.type[j]=="pen")
         {
            if(basis=="trunc.poly")
               ncol.X.val <- deg.val
            if(basis=="tps")
               ncol.X.val <- (deg.val-1)/2

            inst.col.inds <- fc.stt.pos:(fc.stt.pos+ncol.X.val-1)
            fc.stt.pos <- fc.stt.pos + ncol.X.val
  
            # Set knot values

            knots <- object$info$pen$knots[[pen.num]]
            num.knots <- length(knots)

            Z.g.inst <- outer(x.grid,knots,"-")

            if (basis=="trunc.poly")
            {
               if (drv==0)
                  Z.g.inst <- (Z.g.inst*(Z.g.inst>0))^deg.val
               else
               {
                  if (deg.val >= drv)
                  {
                     mfac <- prod((deg.val-drv+1):deg.val)
                     Z.g.inst <- mfac*(Z.g.inst*(Z.g.inst>0))^(deg.val-drv)
                  }
                  else
                     Z.g.inst <- matrix(0,length(x.grid),length(knots))
               }             
            }

            if(basis=="tps")
            {  
               if(drv==0)
               {  
                  Z.g.inst <- abs(Z.g.inst)^deg.val
                  sqrt.Omega <-  object$info$trans.mat[[pen.num]]
                  Z.g.inst <- t(solve(sqrt.Omega,t(Z.g.inst)))
               }
               else
               {
                  if (deg.val>=drv)
                  {
                     mfac <- prod((deg.val-drv+1):deg.val)
                     Z.g.inst <- mfac*(Z.g.inst^(deg.val-drv-1)*abs(Z.g.inst)) 
                     sqrt.Omega <-  object$info$trans.mat[[pen.num]]
                     Z.g.inst <- t(solve(sqrt.Omega,t(Z.g.inst)))
                  }
                  else
                     Z.g.inst <- matrix(0,length(x.grid),length(knots))
                }
             }

             C.g.inst <- cbind(C.g.inst,Z.g.inst)
            
             # Update the column indices for the insertion
                 
             inst.col.inds <- c(inst.col.inds,
                                   block.inds[[j+1]][-(1:ncol.X.val)])
  
             pen.num <- pen.num + 1
          }

          # Perform the insertion and obtain the y.grid

          C.grid[,inst.col.inds] <- C.g.inst

          y.grid <- as.vector(C.grid%*%coefs)

          # Set default value for vertical range

          if (default.ylim)   
          {
             plot.params$ylim$lower[j] <- min(y.grid)
             plot.params$ylim$upper[j] <- max(y.grid)
          }
 
          # Compute standard error bars if requested.

          if (se==TRUE)
          {
             se.grid <- sqrt(diag(C.grid%*%cov.mat%*%t(C.grid)))
       
             lower <- y.grid - 2*se.grid
             upper <- y.grid + 2*se.grid

             if (default.ylim)
             {
                plot.params$ylim$lower[j] <- min(lower)
                plot.params$ylim$upper[j] <- max(upper)
             }

             ylim.val <- range(c(lower,upper))
          }

          if (plot.params$plot.it[j]==TRUE)
          {
             plot(x.grid,y.grid,type="n",bty=plot.params$bty[j],
                  main=plot.params$main[j],xlab=plot.params$xlab[j],
                  ylab=plot.params$ylab[j],
                  xlim=c(plot.params$xlim$lower[j],
                       plot.params$xlim$upper[j]),
                  ylim=c(plot.params$ylim$lower[j],
                       plot.params$ylim$upper[j]))

             if (se==TRUE)
             {
                if (shade==FALSE) 
                {
                   lines(x.grid,lower,lty=plot.params$se.lty[j],
                                      lwd=plot.params$se.lwd[j],
                                      col=plot.params$se.col[j])
                   lines(x.grid,upper,lty=plot.params$se.lty[j],
                                      lwd=plot.params$se.lwd[j],
                                      col=plot.params$se.col[j])
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
    		
                   # The following code ensures that shaded plots
                   # sent to PostScript files have a nice shade of grey.
       	
                   polygon(c(x.grid,rev(x.grid)),c(lower,rev(upper)),
                           col=plot.params$shade.col[j],border=FALSE)  

                }              
            }

            if((drv>=1)&(plot.params$zero.line==TRUE))
               abline(h=0,err=-1)

            lines(x.grid,y.grid,lty=plot.params$lty[j],lwd=plot.params$lwd[j],
                  col=plot.params$col[j])

            if (plot.params$jitter.rug==TRUE)
               rug.x <- jitter(x.vals[,j])

            if (plot.params$jitter.rug==FALSE)
               rug.x <- x.vals[,j]
            
            rug(rug.x,quiet=1,col=plot.params$rug.col[j])
         }
      }
   }

   if (krige.present)
   {  
      if (plot.params$plot.image==TRUE)
      { 
         # Set number of curves to zero if only krige present

         if ((!lin.present)&(!pen.present))
            num.curves <- 0

         # Determine image.grid size and form `underyling'
         # average matrix

         pixel.info <- get.pixel.info(plot.params)

         on.inds <- pixel.info$on.inds
         off.inds <- pixel.info$off.inds

         image.grid.size <- plot.params$image.grid.size

         num.pix <- image.grid.size[1]*image.grid.size[2] 
         # total number of pixels

         num.on.pix <- length(on.inds)    # total number of switched on pixels

         C.grid <-  matrix(rep(ave.vals,num.on.pix),num.on.pix,
                           length(ave.vals),byrow=TRUE)

         # Set vectorised version of z matrix for image()

         z <- rep(0,num.pix)

         # Obtain matrix of (x1,x2) pairs corresponding to
         # the switched on pictures

         X.g.inst <- as.matrix(pixel.info$newdata)

         x1.g.inst <- X.g.inst[,1] ; x2.g.inst <-  X.g.inst[,2]

         m.krige <- (object$info$krige$degree/2) + 1

         if (m.krige>2)
            for (im in 3:m.krige)
            {
               X.g.inst <- cbind(X.g.inst,x1.g.inst^(im-1))
               if (im>2)
                  for (ipow in 2:(im-1))
                     X.g.inst <- cbind(X.g.inst,(x1.g.inst^(im-ipow))
                                       *(x2.g.inst^(ipow-1)))
               X.g.inst <- cbind(X.g.inst,x2.g.inst^(im-1))
            }

         # Obtain mesh-wise Z matrix

         knots <- object$info$krige$knots
         knots.1 <- knots[,1]
         knots.2 <- knots[,2]

         num.knots <- length(knots.1)

         dists.g <- outer(x1.g.inst,knots.1,"-")^2
         dists.g <- sqrt(dists.g + outer(x2.g.inst,knots.2,"-")^2)
         prelim.Z.g.inst <- tps.cov(dists.g,m=m.krige,d=2)
     
         Z.g.inst <-  t(solve(object$info$trans.mat[[num.pen+1]],
                        t(prelim.Z.g.inst)))
     
         C.g.inst <- cbind(X.g.inst,Z.g.inst)

         # Set up the column indices for the insertion
 
         inst.col.inds <- block.inds[[num.curves+2]]

         # Perform the insertion and obtain the values of z
         # for switched on pixels.

         C.grid[,inst.col.inds] <- C.g.inst

         z[on.inds] <- as.vector(C.grid%*%coefs)

         # Set values of z for switched off pixels to NA  

         z[off.inds] <- NA
     
         # Convert to a matrix and obtain image plot 

         z <- matrix(z,image.grid.size[1],image.grid.size[2])

         imageFull(z,plot.params) 
     }
  }
  invisible()
 }

########## end of S-function: plot.spm ##########

