########## R function: predict.spm ##########

# For obtaining predictions and corresponding
# standard errors from an spm() fit object.

# Last changed: 02 FEB 2005

predict.spm <- function(object,newdata,se=FALSE,...)
{
   # Set flags for presence of components of various types

   lin.present <- !is.null(object$info$lin)
   pen.present <- !is.null(object$info$pen)
   krige.present <- !is.null(object$info$krige)
   random.present <- !is.null(object$info$random)

   # Extract block indices

   block.inds <- object$aux$block.inds 

   # Extract coefficients

   if (pen.present|krige.present)
      coefs <- c(object$fit$coef$fixed,object$fit$coef$random)
   if(!(pen.present|krige.present))
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
   
   # Build newdata.lin, newdata.pen and 
   # newdata.krige
   
   num.vars <- ncol(newdata)
   
   newdata.lin <- NULL
   if (lin.present)
      for (curr.name in object$info$lin$name)
      {
         col.num <- (1:num.vars)[names(newdata)==curr.name]
         newdata.lin <- cbind(newdata.lin,newdata[,col.num])
      }

   newdata.pen <- NULL
   if (pen.present)
      for (curr.name in object$info$pen$name)
      { 
         col.num <- (1:num.vars)[names(newdata)==curr.name]
         newdata.pen <- cbind(newdata.pen,newdata[,col.num])
      }
 
   newdata.krige <- NULL
   if (krige.present)
   {
      for (curr.name in object$info$krige$name)
      {
         col.num <- (1:num.vars)[names(newdata)==curr.name]
         newdata.krige <- cbind(newdata.krige,newdata[,col.num])
      }
   }
      
   # Build up newdata design matrix, starting
   # with the intercept.
   
   C.newdata <- as.matrix(rep(1,nrow(newdata)))
   
   if (lin.present) # Add on columns for "lin" components
      C.newdata <- cbind(C.newdata,newdata.lin)

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
     
         new.col <- newdata.pen[,ipen]
         C.newdata <- cbind(C.newdata,new.col)
     
         if (ncol.X.val>1)
            for (ipow in 2:ncol.X.val)
               C.newdata <- cbind(C.newdata,new.col^ipow)
      }
   }

   if (!pen.present)
      num.pen <- 0
   
   if (krige.present) # Add on fixed effect entries for "krige" component
   {
      x1.v <- newdata.krige[,1]
      x2.v <- newdata.krige[,2]
   
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
   
      C.newdata <- cbind(C.newdata,X.kg)
   }
   
   if (pen.present) # Add on random effect entries for "pen" components
   {
      num.pen <- ncol(as.matrix(newdata.pen))
      for (ipen in 1:num.pen)
      {
         deg.val <- object$info$pen$degree[ipen]
         knots <- object$info$pen$knots[[ipen]]
   
         new.col <- as.matrix(newdata.pen)[,ipen]
         if(basis=="trunc.poly")
         {
            new.cols.spline <- outer(new.col,knots,"-")
            new.cols.spline <- (new.cols.spline*(new.cols.spline>0))^deg.val
         }
         if (basis=="tps")
         {  
            new.cols.spline <- abs(outer(new.col,knots,"-"))^deg.val
            sqrt.Omega <-  object$info$trans.mat[[ipen]]
            new.cols.spline <- t(solve(sqrt.Omega,t(new.cols.spline))) 
         }
         C.newdata <- cbind(C.newdata,new.cols.spline)
      }
   }
   
   if (krige.present) # Add on random effect entries for "krige" component
   {
       new.col.x1 <- newdata.krige[,1]
       new.col.x2 <- newdata.krige[,2]
   
       knots <- object$info$krige$knots
       knots.1 <- knots[,1]
       knots.2 <- knots[,2]
   
       num.knots <- length(knots.1)
   
       new.cols.spline <- NULL
   
       dists.mat <- sqrt(outer(new.col.x1,knots.1,"-")^2
                                +outer(new.col.x2,knots.2,"-")^2)         

       new.cols.spline <- tps.cov(dists.mat,m=m.krige,d=2)
   
       sqrt.Omega <-object$info$trans.mat[[length(object$info$trans.mat)]]

       new.cols.spline <- t(solve(sqrt.Omega,t(new.cols.spline)))

       C.newdata <- cbind(C.newdata,new.cols.spline) 
   }
 
   fit <- as.vector(C.newdata%*%coefs)
   
   if (se==TRUE)
   {
      se.vec <- sqrt(diag(C.newdata%*%cov.mat%*%t(C.newdata)))
      return(list(fit=fit,se=se.vec))  
   }

   if (se==FALSE)
      return(fit)
}

########## End of predict.spm ##########