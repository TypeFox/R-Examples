########## S-function: spmDesign ##########

# For creating the design matrices for
# a semiparametric model.

spmDesignOS <- function(spm.info)
{  
   # Extract required information 

   degree.pen <- spm.info$pen$degree



   degree.krige <- spm.info$krige$degree
   knots.pen <- spm.info$pen$knots
   basis.pen <- spm.info$pen$basis
   x.pen <- spm.info$pen$x

  if (!is.null(basis.pen))
    for (j in 1:length(degree.pen))
      if (basis.pen=="os" & length(degree.pen[[j]])==1){
        if ((degree.pen[[j]]+1)/2%%1 != 0) {warning("If only degree of B-spline basis is given, degree must be chosen such that q=(p+1)/2 is an integer. Set to its default p=3, q=2."); degree.pen[[j]] =c(3,2)}
        else degree.pen[[j]] <- c(degree.pen[[j]], degree.pen[[j]])
      }

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
      if(basis.pen=="os") ncol.X.pen <- lapply(ncol.X.pen,function(x) x[2]-1)
      if(basis.pen=="tps")
      {
         ncol.X.pen <-lapply(ncol.X.pen,function(x) (x-1)/2) #(ncol.X.pen-1)/2
         if (any(sapply(ncol.X.pen,floor)!=ncol.X.pen))
            stop("Only odd degree thin plate splines supported.")
      }

      for (j in 1:ncol(as.matrix(spm.info$pen$x)))
      {
         block.inds <- c(block.inds,list(stt.pos:(stt.pos+ncol.X.pen[[j]]-1)))
         stt.pos <- stt.pos + ncol.X.pen[[j]]
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

      X <- cbind(X,as.matrix(spm.info$lin$x))

      comp.num <- comp.num + ncol(X)
   } 

   if (!is.null(spm.info$pen))
   {  

      for (j in 1:ncol(as.matrix(x.pen)))
      {   

if (basis.pen=="tps"|basis.pen=="trunc.poly"){
         # Update the X matrix.
          
        Xnew=numeric()
         for (ideg in 1:ncol.X.pen[[j]])
            Xnew <- cbind(Xnew,as.matrix(x.pen)[,j]^ideg)
#            X <- cbind(X,as.matrix(x.pen)[,j]^ideg)
############
#centering
n=nrow(Xnew)
colsum= (colSums(Xnew))
if (length(colsum)==1) Xnew=matrix(t(apply(Xnew,1,function(x) x-colsum/n)))
else Xnew=t(apply(Xnew,1,function(x) x-colsum/n))
#Xnew=(diag(rep(1,n))-1/n)%*%Xnew
############
X=cbind(X,Xnew)

         # Update the Z matrix.

         if (basis.pen=="tps")
         {  
            new.cols <- outer(as.matrix(x.pen)[,j],knots.pen[[j]],"-")  
            new.cols <- abs(new.cols)^degree.pen[[j]]

            sqrt.Omega <- matrix.sqrt(abs(outer(knots.pen[[j]],
                            knots.pen[[j]],"-"))^degree.pen[[j]])

            new.cols <- t(solve(sqrt.Omega,t(new.cols)))  
            trans.mat <- c(trans.mat,list(sqrt.Omega))
         }

         if (basis.pen=="trunc.poly")
         {
            new.cols <- outer(as.matrix(x.pen)[,j],knots.pen[[j]],"-")
            new.cols <- (new.cols*(new.cols>0))^degree.pen[[j]]
            trans.mat <- c(trans.mat,list(NULL))
         }

          ############
          #centering
            colsum= colSums(new.cols)
            new.cols=t(apply(new.cols,1,function(x) x-colsum/n))
          #############

          Z <- cbind(Z,new.cols)
        }
        else {
          if (basis.pen!="os") stop("Currently only trunc.poly, tps, os basis functions allowed.")

          m=  spm.info$pen$degree[[j]]
          nk=length(knots.pen[[j]])

          xx= as.matrix(x.pen)[,j]
          data.xx = data.frame(x = xx)
          names(data.xx) <- spm.info$pen$name[j]

          knotsx= seq(min(xx),max(xx),length=nk+2)[-c(1,nk+2)]
          names(knots) <- NULL

          smooth=list(m=m,bs.dim=nk,term=spm.info$pen$name[j],p.order=m,knots=knotsx)
          class(smooth)= "ospline.smooth"


          CZ.temp <- Predict.matrix.lme(smooth,data.xx)
          Z <- cbind(Z,CZ.temp$Z)
          if (m[1]!=1) X <- cbind(X,CZ.temp$C)

          trans.mat <- c(trans.mat,list(NULL))
          knots.pen[[j]]=CZ.temp$knots
        }


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

      # Transform mixed model to canonical form 
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


   if (length(spm.info$y)==1) spm.info$y=rep(spm.info$y,length(spm.info$pen$x))
   if (intercept.only)
      X <- as.matrix(rep(1,length(spm.info$y)))

   if (!intercept.only)
      X <- cbind(rep(1,length(spm.info$y)),X)

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

  knots.pen -> spm.info$pen$knots

   return(list(spm.info=spm.info,X=X,Z=Z,Z.spline=Z.spline,
               trans.mat=trans.mat,block.inds=block.inds,
               re.block.inds=re.block.inds))   
                
}

########## End of spmDesign ##########







