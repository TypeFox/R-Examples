########## R function: predict.spm ##########

# For obtaining predictions and corresponding
# standard errors from an asp2() fit object.


predict.asp <- function(object,newdata,se=FALSE,...)
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
      for (curr.name in object$info$lin$name.orig) {
        
        if (is.factor(newdata[[curr.name]])){
          temp=model.matrix(as.formula(paste("~",curr.name)),data=newdata,constrasts=object$info$lin$contrasts[[curr.name]])
          temp=temp[,-1,drop=F]
          temp=as.matrix(temp)
          dimnames(temp)[[2]]=NULL
          newdata.lin <- cbind(newdata.lin,temp) 
        }
        else {
          col.num <- (1:num.vars)[names(newdata)==curr.name]
          newdata.lin <- cbind(newdata.lin,newdata[,col.num])
        }
        
        
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

   if (pen.present){ # Add on fixed effect entries for "pen" components

      if (basis=="tps"|basis=="trunc.poly"){
        num.pen <- ncol(as.matrix(object$info$pen$x))
        for (ipen in 1:num.pen){
           deg.val <- object$info$pen$degree[[ipen]]

           if(basis=="trunc.poly")
              ncol.X.val <- deg.val
           if(basis=="tps")
              ncol.X.val <- (deg.val-1)/2

           new.col <- newdata.pen[,ipen]
           Cnew=as.matrix(new.col)
#           C.newdata <- cbind(C.newdata,new.col)

           if (ncol.X.val>1)
              for (ipow in 2:ncol.X.val)
                 Cnew <- cbind(Cnew,new.col^ipow)

          ############
          #centering
          original= as.matrix(object$info$pen$x)[,ipen,drop=F]
            if (ncol.X.val>1)
              for (ipow in 2:ncol.X.val)
                 original <- cbind(original,object$info$pen$x[,ipen]^ipow)
          n=nrow(original)
          colsum= (colSums(original))
          if (length(colsum)==1) Cnew=matrix(t(apply(Cnew,1,function(x) x-colsum/n)))
          else Cnew=t(apply(Cnew,1,function(x) x-colsum/n))
          #Xnew=(diag(rep(1,n))-1/n)%*%Xnew
          ############

          C.newdata <- cbind(C.newdata,Cnew)
        }
      }else if (basis=="os"){
        num.pen <- ncol(as.matrix(object$info$pen$x))
        for (ipen in 1:num.pen){
          m<- object$info$pen$degree[[ipen]]

          nk=length(object$info$pen$knots[[ipen]])

          xx= as.matrix(newdata.pen)[,ipen]
          data.xx = data.frame(x = xx)
          names(data.xx) <- object$info$pen$name[ipen]

          knotsx= seq(min(xx),max(xx),length=nk+2)[-c(1,nk+2)]
          names(knots) <- NULL

          smooth=list(m=m,bs.dim=nk,knots= object$info$pen$knots[[ipen]],term=object$info$pen$name[ipen],p.order=m,knots=knotsx)
          class(smooth)= "ospline.smooth"

          CZ.temp <- Predict.matrix.lme(smooth,data.xx,center=F)
          
          # centering with respect to curve at data points
            data1 <- data.frame(x= (object$info$pen$x)[,ipen,drop=F])
            names(data1) <- object$info$pen$name[ipen]
            CZj1= Predict.matrix.lme(smooth,data1,center=F)
            if (class(smooth)=="ospline.smooth") Cj1= CZj1$C[,drop=F]  else Cj1= CZj1$C[,-1,drop=F]
            CZ.temp$C=t(apply(CZ.temp$C,1,function(x) x-colSums(Cj1)/nrow(Cj1)))[,,drop=T]

          if (m[1]!=1) C.newdata <- cbind(C.newdata,CZ.temp$C)
        }
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
   
   if (pen.present){ # Add on random effect entries for "pen" components
      if (basis=="tps"|basis=="trunc.poly"){
        num.pen <- ncol(as.matrix(newdata.pen))
        for (ipen in 1:num.pen){
           deg.val <- object$info$pen$degree[[ipen]]
           knots <- object$info$pen$knots[[ipen]]

           new.col <- as.matrix(newdata.pen)[,ipen]
           if(basis=="trunc.poly"){
              new.cols.spline <- outer(new.col,knots,"-")
              new.cols.spline <- (new.cols.spline*(new.cols.spline>0))^deg.val
           }
           if (basis=="tps"){
              new.cols.spline <- abs(outer(new.col,knots,"-"))^deg.val
              sqrt.Omega <-  object$info$trans.mat[[ipen]]
              new.cols.spline <- t(solve(sqrt.Omega,t(new.cols.spline)))
           }
#          C.newdata <- cbind(C.newdata,new.cols.spline)
            Cnew=new.cols.spline

          ############
          #centering
           new.col <- object$info$pen$x[,ipen]
            if(basis=="trunc.poly"){
              new.cols.spline <- outer(new.col,knots,"-")
              new.cols.spline <- (new.cols.spline*(new.cols.spline>0))^deg.val
           }
           if (basis=="tps"){
              new.cols.spline <- abs(outer(new.col,knots,"-"))^deg.val
              sqrt.Omega <-  object$info$trans.mat[[ipen]]
              new.cols.spline <- t(solve(sqrt.Omega,t(new.cols.spline)))
           }


          n=nrow(new.cols.spline)
          colsum= (colSums(new.cols.spline))
          if (length(colsum)==1) Cnew=matrix(t(apply(Cnew,1,function(x) x-colsum/n)))
          else Cnew=t(apply(Cnew,1,function(x) x-colsum/n))
          #Xnew=(diag(rep(1,n))-1/n)%*%Xnew
          ############

          C.newdata <- cbind(C.newdata,Cnew)
        }
      }else if (basis=="os"){
        num.pen <- ncol(as.matrix(object$info$pen$x))
        for (ipen in 1:num.pen){
          m<- object$info$pen$degree[[ipen]]

          nk=length(object$info$pen$knots[[ipen]])

          xx= as.matrix(newdata.pen)[,ipen]
          data.xx = data.frame(x = xx)
          names(data.xx) <- object$info$pen$name[ipen]

          knotsx= seq(min(xx),max(xx),length=nk+2)[-c(1,nk+2)]
          names(knots) <- NULL

          smooth=list(m=m,bs.dim=nk,knots= object$info$pen$knots[[ipen]],term=object$info$pen$name[ipen],p.order=m,knots=knotsx)
          class(smooth)= "ospline.smooth"

          CZ.temp <- Predict.matrix.lme(smooth,data.xx,center=F)

        # centering with respect to curve at data points
            data1 <- data.frame(x= (object$info$pen$x)[,ipen,drop=F])
            names(data1) <- object$info$pen$name[ipen]
            CZj1= Predict.matrix.lme(smooth,data1,center=F)
            Zj1= CZj1$Z
            CZ.temp$Z=t(apply(CZ.temp$Z,1,function(x) x-colSums(Zj1)/nrow(Zj1)))



          C.newdata <- cbind(C.newdata,CZ.temp$Z)

        #  C.newdata <- cbind(C.newdata,Cnew)
        }
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