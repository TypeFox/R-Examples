########## S-function: spmDesign ##########

# For creating the design matrices for
# a semiparametric model.

# Last changed: 22 SEP 2005

spmDesign <- function(spm.info)
{  
   # Extract required information 

   degree.pen <- spm.info$pen$degree
   degree.krige <- spm.info$krige$degree
   knots.pen <- spm.info$pen$knots
   basis.pen <- spm.info$pen$basis
   x.pen <- spm.info$pen$x

   # Initialise matrices and counters

   X <- NULL
   Z <- NULL    
   block.inds <- NULL  
   re.block.inds <- NULL
   trans.mat <- NULL
   stt.pos <- 1
   re.stt.pos <- 1

   # Obtain X-matrix components corresponding 
   # to the "block.inds", starting with the
   # intercept.

   block.inds <- 1
   stt.pos <- stt.pos + 1

   intercept.only <- (is.null(spm.info$lin)&is.null(spm.info$pen)
                      &is.null(spm.info$krige))

   if (!is.null(spm.info$lin))
   {
      for (j in 1:ncol(as.matrix(spm.info$lin$x)))
         block.inds <- c(block.inds,list(j+stt.pos-1))
    
      stt.pos <- stt.pos + ncol(as.matrix(spm.info$lin$x))
   }

   if (!is.null(spm.info$pen))
   { 

      ncol.X.pen <- degree.pen
      if(basis.pen=="tps")
      {
         ncol.X.pen <- (ncol.X.pen-1)/2            
         if (any(floor(ncol.X.pen)!=ncol.X.pen)) 
            stop("Only odd degree thin plate splines supported.")
      }

      for (j in 1:ncol(as.matrix(spm.info$pen$x)))
      {
         block.inds <- c(block.inds,list(stt.pos:(stt.pos+ncol.X.pen[j]-1)))
         stt.pos <- stt.pos + ncol.X.pen[j]
      }

   }

   if (!is.null(spm.info$krige))
   {
      ncol.X.krige <- (degree.krige+2)*(degree.krige+4)/8 - 1

      block.inds <- c(block.inds,list(stt.pos:(stt.pos+ncol.X.krige-1)))
      stt.pos <- stt.pos + ncol.X.krige
   }

   # Now obtain the X and Z matrices, and update
   # "block.inds", "re.block.inds" and "trans.mat" sequentially.

   comp.num <- 2
   if (!is.null(spm.info$lin))
   {
      # Update the X matrix

      X <- cbind(X,spm.info$lin$x)

      comp.num <- comp.num + ncol(X)
   } 

   if (!is.null(spm.info$pen))
   {  

      for (j in 1:ncol(as.matrix(x.pen)))
      {   
         # Update the X matrix.
          
         for (ideg in 1:ncol.X.pen[j])
            X <- cbind(X,as.matrix(x.pen)[,j]^ideg)

         # Update the Z matrix.

         if (basis.pen=="tps")
         {  
            new.cols <- outer(as.matrix(x.pen)[,j],knots.pen[[j]],"-")  
            new.cols <- abs(new.cols)^degree.pen[j]

            sqrt.Omega <- matrix.sqrt(abs(outer(knots.pen[[j]],
                            knots.pen[[j]],"-"))^degree.pen[j])

            new.cols <- t(solve(sqrt.Omega,t(new.cols)))  
            trans.mat <- c(trans.mat,list(sqrt.Omega))
         }

         if (basis.pen=="trunc.poly")
         {
            new.cols <- outer(as.matrix(x.pen)[,j],knots.pen[[j]],"-")
            new.cols <- (new.cols*(new.cols>0))^degree.pen[j]
            trans.mat <- c(trans.mat,list(NULL))
         }

         Z <- cbind(Z,new.cols)

         re.block.inds <- c(re.block.inds,list(re.stt.pos:ncol(Z)))

         end.pos <- ncol(Z) + stt.pos - re.stt.pos
         block.inds[[comp.num]] <- c(block.inds[[comp.num]],stt.pos:end.pos)

         re.stt.pos <- ncol(Z) + 1
         stt.pos <-  end.pos + 1
         comp.num <- comp.num + 1
      }
   } 

   if (!is.null(spm.info$krige))
   {
      # Now add on the kriging contribution, starting with
      # obtaining the "Omega matrix" of intra-knot covariances,
      # but first work out range parameter.

      knots.krige <- spm.info$krige$knots

      # Update the X matrix

      x1.v <- spm.info$krige$x[,1]
      x2.v <- spm.info$krige$x[,2]

      m.krige <- (degree.krige/2) + 1

      for (im in 2:m.krige)
      {
         X <- cbind(X,x1.v^(im-1))
         if (im>2)
            for (ipow in 2:(im-1))
               X <- cbind(X,(x1.v^(im-ipow))*(x2.v^(ipow-1)))
         X <- cbind(X,x2.v^(im-1))
      }

      # Update the Z matrix

      num.knots.krige <- nrow(knots.krige)

      dist.mat <- matrix(0,num.knots.krige,num.knots.krige)
      dist.mat[lower.tri(dist.mat)] <- dist(as.matrix(knots.krige))
      dist.mat <- dist.mat + t(dist.mat)

      Omega <- tps.cov(dist.mat,m=m.krige,d=2)

      # Obtain preliminary Z matrix of knot to data covariances,
      # for the kriging component.

      x.knot.diffs.1 <- outer(spm.info$krige$x[,1],knots.krige[,1],"-")
      x.knot.diffs.2 <- outer(spm.info$krige$x[,2],knots.krige[,2],"-")
      x.knot.dists <- sqrt(x.knot.diffs.1^2+x.knot.diffs.2^2)

      prelim.Z <- tps.cov(x.knot.dists,m=m.krige,d=2)

      # Transform to mixed model canonical form 
      # (covariance matrix of random effects is multiple
      # of identity)

      sqrt.Omega <- matrix.sqrt(Omega)

      trans.mat <- c(trans.mat,list(sqrt.Omega))

      # Combine to form final Z matrix

      Z <- cbind(Z,t(solve(sqrt.Omega,t(prelim.Z))))

      re.block.inds <- c(re.block.inds,list(re.stt.pos:ncol(Z)))        

      end.pos <- ncol(Z) + stt.pos - re.stt.pos
      block.inds[[comp.num]] <- c(block.inds[[comp.num]],stt.pos:end.pos)

      re.stt.pos <- ncol(Z) + 1
      stt.pos <- end.pos + 1
      comp.num <- comp.num + 1
   }

   # Add on column of 1's for the intercept term

   if (intercept.only)
      X <- as.matrix(rep(1,length(spm.info$y)))

   if (!intercept.only)
      X <- cbind(rep(1,nrow(X)),X)

   # Save spline component of Z separately.

   Z.spline <- Z

   # If random term present then add on Kronecker-type
   # intercept structure to Z matrix

   if (!is.null(spm.info$random))
   {
      group.inds <- spm.info$random$group.inds

      num.groups <- length(group.inds)

      Z.add <- matrix(0,length(unlist(group.inds)),num.groups)
   
      for (j in 1:num.groups)
         Z.add[group.inds[[j]],j] <- 1

      # Append to end of Z matrix

      Z <- cbind(Z,Z.add)

      re.block.inds <- c(re.block.inds,list(re.stt.pos:ncol(Z)))  

      end.pos <- ncol(Z) + stt.pos - re.stt.pos

      block.inds <- c(block.inds,list(stt.pos:end.pos))

      trans.mat <- c(trans.mat,list(NULL))
   }

   if (is.null(spm.info$krige)) sqrt.Omega.krige <- NULL

   # Return spm() fit object

   return(list(spm.info=spm.info,X=X,Z=Z,Z.spline=Z.spline,
               trans.mat=trans.mat,block.inds=block.inds,
               re.block.inds=re.block.inds))   
                
}

########## End of spmDesign ##########







