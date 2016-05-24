########## R-function: spm ##########

# For fitting semiparametric regression models.

# Last changed: 21 NOV 2005 by MPW


spm <- function(form,random=NULL,group=NULL,family="gaussian",
               spar.method="REML",omit.missing=NULL)
{
  X.Declan <- NULL
  rm(X.Declan)
  Z.Jaida <- NULL
  rm(Z.Jaida)
  Z.spline.Jaida <- NULL
  rm(Z.spline.Jaida)
  
  e <- new.env()
##   require("nlme")

   # Read in primary information

   spm.info <- spmFormRead(form,omit.missing)

   # Read in random effect information if present

   
   random.info <- NULL
   if (!is.null(random))
   {
      random.info <- random.read(random,group)
      spm.info <- c(spm.info,list(random=random.info))
      group <- as.numeric(factor(group))
   }

   # Put in safeguard against zero smoothing parameters.

   if (!is.null(unlist(spm.info$pen$spar)))
   {
      if (any(unlist(spm.info$pen$spar)==0))
         stop("zero smoothing parameters not supported in current version.")
   }

   # Create required matrices

   design.info <- spmDesign(spm.info)

   X <- design.info$X
   Z <- design.info$Z
   Z.spline <- design.info$Z.spline
   y <- spm.info$y

   # Add on transformation matrix information

   trans.mat <- design.info$trans.mat
   spm.info <- c(spm.info,list(trans.mat=trans.mat))

   block.inds <- design.info$block.inds

   re.block.inds <- design.info$re.block.inds
   col.ones <- rep(1,nrow(X))

   # Determine if automatic smoothing parameter selection to be done

   auto.spar.select <- FALSE

   if ((is.null(spm.info$pen))&(is.null(spm.info$krige)))
      auto.spar.select <- TRUE

   if (!is.null(spm.info$pen))
   {
      auto.spar <- 0

      for (j in 1:length(spm.info$pen$name))
         auto.spar <- (auto.spar +
         (is.null(spm.info$pen$spar[[j]])&(spm.info$pen$adf[[j]]=="miss")))

       if (auto.spar > 0)
          auto.spar.select <- TRUE
   }

   if (!is.null(spm.info$krige))
   {
      if ((is.null(spm.info$krige$spar))
           &(spm.info$krige$adf[[1]]=="miss"))
        auto.spar.select <- TRUE
   }

   # If not automatic smoothing parameter selection do
   # "df-to-spar" conversion if required, and create
   # G-matrix

   if (auto.spar.select==FALSE)
   {
      # If "adf" specified then convert to spar values

      if (!is.null(spm.info$lin))
         num.lin <- ncol(as.matrix(spm.info$lin$x))

      if (is.null(spm.info$lin))
         num.lin <- 0

      compon.num <- 1
      if (!is.null(spm.info$pen))
      {
         basis.type <- spm.info$pen$basis

         for (j in 1:ncol(as.matrix(spm.info$pen$x)))
         {
            if (!is.null(spm.info$pen$adf[[j]])&
                 (spm.info$pen$adf[[j]]!="miss"))
            {

               deg.val <- spm.info$pen$degree[j]

               # Extract X matrix for current component

               stt.ind <- block.inds[[compon.num+num.lin+1]][1]

               ncol.Xj <- length(block.inds[[compon.num+num.lin+1]])
               ncol.Xj <- ncol.Xj - length(re.block.inds[[compon.num]])

               Xj <- X[,stt.ind:(stt.ind+ncol.Xj-1)]

               # Add column of ones for intercept

               Xj <- cbind(col.ones,Xj)

               # Extract Z matrix for current component

               Zj <- Z[,re.block.inds[[compon.num]]]

               adf.val <- spm.info$pen$adf[[j]]

               if (family=="gaussian")
                  spar.val <- df.to.spar(adf.val+1,Xj,Zj)

               if (!family=="gaussian")
                  spar.val <- glm.df.to.spar(adf.val+1,y,Xj,Zj,family)

               if (basis.type=="trunc.poly")
                   spm.info$pen$spar[[j]] <- exp(log(spar.val)/(2*deg.val))
               else
                   spm.info$pen$spar[[j]] <- exp(log(spar.val)/deg.val)

               compon.num <- compon.num + 1
            }
         }
      }
      if (!is.null(spm.info$krige))
      {
         if ((!is.null(spm.info$krige$adf))&(spm.info$krige$adf[[1]]!="miss"))
         {
            deg.val <- spm.info$krige$degree

            # Extract X matrix for current component

            stt.ind <- block.inds[[compon.num+num.lin+1]][1]
            ncol.Xj <- length(block.inds[[compon.num+num.lin+1]])
            ncol.Xj <- ncol.Xj - length(re.block.inds[[compon.num]])

            Xj <- X[,stt.ind:(stt.ind+ncol.Xj-1)]

            # Add column of ones for intercept

            Xj <- cbind(col.ones,Xj)

            Zj <- Z[,re.block.inds[[compon.num]]]
            adf.val <- spm.info$krige$adf[[1]]

            if (family=="gaussian")
               spar.val <- df.to.spar(adf.val+1,Xj,Zj)

            if (!family=="gaussian")
               if (is.null(spm.info$off.set))
                  spar.val <- glm.df.to.spar(adf.val+1,y,Xj,Zj,family)
               else
                  spar.val <- glm.df.to.spar(adf.val+1,y-spm.info$off.set,
                                             Xj,Zj,family)

            spm.info$krige$spar <- exp(log(spar.val)/deg.val)
         }
      }

      # Create G matrix

      diag.G <- NULL
      if (!is.null(spm.info$pen))
         for (j in 1:ncol(as.matrix(spm.info$pen$x)))
         {
            deg.val <- spm.info$pen$degree[j]
            spar.val <- spm.info$pen$spar[[j]]
            num <- length(spm.info$pen$knots[[j]])
            if (basis.type=="trunc.poly")
               diag.G <- c(diag.G,
                           rep(1/(exp((2*deg.val)*log(spar.val))),num))
            else
               diag.G <- c(diag.G,
                           rep(1/(exp((deg.val)*log(spar.val))),num))
         }

      if (!is.null(spm.info$krige))
       {
         spar.val <- spm.info$krige$spar
         num.knots <- nrow(spm.info$krige$knots)
         diag.G <- c(diag.G,rep((1/spar.val^2),num.knots))
      }
   }

   # Set up inputs for mixed model functions.

   dummy.group.vec <- col.ones
   fdummy.group.vec <- factor(dummy.group.vec)

   # Assign unusual names to formula variables
   # to overcome the `scoping' problem.

   
   
   assign("dummy.group.vec.Handan",dummy.group.vec,envir = as.environment(e))
   assign("fdummy.group.vec.Handan",fdummy.group.vec,envir = as.environment(e))
   assign("group.Handan",group,envir = as.environment(e))
   assign("X.Declan",X,envir = as.environment(e))
 
 
 X.Declan <<- X

  dummy.group.vec.Handan <- get("dummy.group.vec.Handan" , envir = as.environment(e))
	
   if (!is.null(Z))
   {
     X.Declan <<- X
     Z.Jaida <<- Z
     Z.spline.Jaida <<- Z.spline
        
      if (is.null(random)) # Calls to lme() and glmmPQL() are of 1 type.
      {
        X.Declan <<- X
        Z.spline.Jaida <<- Z.spline
        
        Z.Jaida <<- Z
        
               
        
         data.fr <- groupedData(y~-1 + X.Declan | dummy.group.vec.Handan,
                    data=data.frame(y,X.Declan,Z.Jaida,
                    dummy.group.vec.Handan))

         Z.block  <-  list()
         for (i in 1:length(re.block.inds))
            Z.block[[i]] <- as.formula(paste("~Z.Jaida[,c(",paste(
                                 re.block.inds[[i]],collapse=","),")]-1"))

         if (length(re.block.inds)==1)
         {
            if (family=="gaussian")
               lme.fit <- nlme::lme(y~-1+X.Declan,random=pdIdent(~-1+Z.Jaida),
                              data=data.fr,method=spar.method)
            if (family!="gaussian")
            {
               lme.fit <- MASS::glmmPQL(y~-1+X.Declan,
               random=list(dummy.group.vec.Handan=pdIdent(~-1+Z.Jaida)),
                           data=data.fr,family=family)
            }
         }

         if (length(re.block.inds)>1)
         {
            if (family=="gaussian")
               lme.fit <- nlme::lme(y~-1+X.Declan,
                           random=list(dummy.group.vec.Handan=
                           pdBlocked(Z.block,
                           pdClass=rep("pdIdent",length(Z.block)))),
                           data=data.fr,method=spar.method)

            if (family!="gaussian")
            {
               lme.fit <- MASS::glmmPQL(y~-1+X.Declan,
                           random=list(dummy.group.vec.Handan=
                           pdBlocked(Z.block,pdClass=
                           rep("pdIdent",length(Z.block)))),
                           data=data.fr,family=family)
            }
         }
      }

      if (!is.null(random)) # Calls to lme() and glmmPQL() are of other type.
      {
         
         
         Z.Jaida <<- Z
         
     
         
         Z.spline.Jaida <<- Z.spline
         
        
         
         data.fr <- groupedData(y~-1+X.Declan|dummy.group.vec.Handan,
                    data=data.frame(y,X.Declan,Z.Jaida,
                    dummy.group.vec.Handan))

         Z.block  <-  list()
         for (i in 1:length(re.block.inds))
            Z.block[[i]] <- as.formula(paste("~Z.Jaida[,c(",paste(
                                 re.block.inds[[i]],collapse=","),")]-1"))

         if (length(re.block.inds)==1)  # Must be parametric random
                                        # effects model.
         {
            if (family=="gaussian")
               stop("implement later; use pigs to test")

            if (family!="gaussian")
            {
               stop("implement later; use Poisson simul to test")
            }
         }

         if (length(re.block.inds)>1)
         {
           
            if (family=="gaussian")
            {

               Z.block <- list(list(fdummy.group.vec.Handan
                               =pdIdent(~Z.spline.Jaida-1)),
                               list(group.Handan=pdIdent(~1)))

               Z.block <- unlist(Z.block,recursive=FALSE)
               data.fr <- groupedData(y~-1+X.Declan|fdummy.group.vec.Handan,
                          data=data.frame(y,X.Declan,Z.spline.Jaida,group))


               lme.fit <- nlme::lme(y~-1+X.Declan,data=data.fr,random=Z.block,
                              method=spar.method)
            }

            if (family!="gaussian")
            {
            
               lme.fit <- MASS::glmmPQL(y~-1+X.Declan,
                           random=list(dummy.group.vec.Handan=
                           pdBlocked(Z.block,pdClass=
                           rep("pdIdent",length(Z.block)))),
                           data=data.fr,family=family)
            }
         }
      }

      # Add sigma to object

      lme.fit <- c(lme.fit,list(sigma=summary(lme.fit)$sigma))
   }

   if (is.null(Z))
   {
     
     
      data.fr <- cbind(y,X.Declan,dummy.group.vec.Handan)
      G <- NULL

      if (family=="gaussian")
      {
     

         lm.fit <- lm(y~-1+X.Declan)
         lme.fit <- list(coef=list(fixed=lm.fit$coef),
                         sigma=summary(lm.fit)$sigma)
      }

      if (family!="gaussian")
      {
      
         if (!is.null(spm.info$off.set))
         {
      
            if (!is.null(X))
               glm.fit <- glm(y~-1+X.Declan,
                              offset=spm.info$off.set,family=family)
            if (is.null(X))
               glm.fit <- glm(y~1,offset=spm.info$off.set,family=family)
         }

         if (is.null(spm.info$off.set))
         {
 
            if (!is.null(X))
               glm.fit <- glm(y~-1+X.Declan,family=family)
            if (is.null(X))
               glm.fit <- glm(y~1,family=family)
         }

         # Build lme.fit object

         lme.fit <- list()
         lme.fit$coef$fixed <-  glm.fit$coef
         lme.fit$coef$random <- NULL
         lme.fit$loglik <- NULL
      }
   }

   # Coerce the random coefficients into an array
   # and store the G matrix.

   if (!is.null(Z))
   {
      lme.fit$coef$random <- unlist(lme.fit$coef$random)
      sig.u.hat <- lme.fit$sigma*exp(unlist(lme.fit$modelStruct))

      diag.sqrt.G <- NULL
      for (ib in 1:length(re.block.inds))
         diag.sqrt.G <- c(diag.sqrt.G,
                            rep(sig.u.hat[ib],length(re.block.inds[[ib]])))

      G <- diag(diag.sqrt.G^2)
   }

   # Store the smoothing parameters and
   # the REML-based estimate of the residual variance

   resid.var <- lme.fit$sigma^2

   if (auto.spar.select==FALSE)
   {
      if (family=="gaussian")
      {
         # Obtain  required QR decomposition

         G <- resid.var*diag(diag.G)

         qr.out <- lmeFitQr(y,X,Z,G,resid.var=resid.var)

         coef.ests <- qr.out$coefficients[1:(ncol(X)+ncol(Z))]

         # Build lme.fit object

         lme.fit <- list()

         lme.fit$coef$fixed <- coef.ests[1:ncol(X)]
         lme.fit$coef$random <- coef.ests[(1+ncol(X)):length(coef.ests)]
      }

      if ((family!="gaussian")&(!is.null(Z)))
      {

         G <- diag(diag.G)

         C.mat <- cbind(X,Z)
         ridge.vec <- c(rep(0,ncol(X)),1/diag.G)

         if (!is.null(spm.info$off.set))
            ridge.reg.fit <- irls.ridge(C.mat,y,off.var=spm.info$off.set,
                                        ridge.vec=ridge.vec,
                                        max.it=50,acc=0.000001,family=family)
         else
            ridge.reg.fit <- irls.ridge(C.mat,y,ridge.vec=ridge.vec,
                                   max.it=50,acc=0.000001,family=family)

         # Build lme.fit object

         lme.fit <- list()

         lme.fit$coef$fixed <-  ridge.reg.fit$coef[1:ncol(X)]
         lme.fit$coef$random <- ridge.reg.fit$coef[(1+ncol(X)):ncol(C.mat)]

         lme.fit$loglik <- NULL
      }
   }

   if (auto.spar.select==TRUE)  # Store the smoothing parameters
   {
      if ((!is.null(spm.info$pen))|(!is.null(spm.info$krige)))
      {
         sigu2.hat <- rep(0,length(re.block.inds))
         if(!is.null(Z))
            for (ib in 1:length(re.block.inds))
               sigu2.hat[ib] <- diag(G)[re.block.inds[[ib]][1]]

         if (is.null(spm.info$krige))  # pen only
         {
            basis.type <- spm.info$pen$basis
            for (ip in 1:ncol(as.matrix(spm.info$pen$x)))
            {
               deg.val <- spm.info$pen$degree[ip]
               if (basis.type=="trunc.poly")
                  spm.info$pen$spar[[ip]] <-
                     exp(log(resid.var/sigu2.hat[ip])/(2*deg.val))
               else
                  spm.info$pen$spar[[ip]] <-
                     exp(log(resid.var/sigu2.hat[ip])/deg.val)
            }
         }

         if (is.null(spm.info$pen))   # krige only
         {
            deg.val <- spm.info$krige$degree
            spm.info$krige$spar <- exp(log(resid.var/sigu2.hat)/deg.val)
         }

         if ((!is.null(spm.info$pen))&(!is.null(spm.info$krige))) # both
         {
            basis.type <- spm.info$pen$basis
            for (ip in 1:ncol(as.matrix(spm.info$pen$x)))
            {
               deg.val <- spm.info$pen$degree[ip]
               var.rats <- (resid.var/sigu2.hat[ip])
               if (basis.type=="trunc.poly")
                  spm.info$pen$spar[[ip]] <- var.rats^(1/(2*deg.val))
               else
                  spm.info$pen$spar[[ip]] <- var.rats^(1/(deg.val))
            }

            num.pen <-  ncol(as.matrix(spm.info$pen$x))
            var.rat <- (resid.var/sigu2.hat[num.pen+1])
            deg.val <- spm.info$krige$degree
            spm.info$krige$spar <- exp(log(var.rat)/deg.val)
         }
      }
   }

   # Obtain auxiliary information

   if (family=="gaussian")
      aux.info <- lmeAux(X,Z,G,resid.var,block.inds)

   if (family!="gaussian")
   {
      if (is.null(Z))
      {
         R <- glm.fit$R
         rinv <- backsolve(R,diag(ncol(X)))
         cov.mat <- rinv%*%t(rinv)
         df.fit <- ncol(X)
         df.res <- length(y) - df.fit
         df <- rep(1,ncol(X))
         random.var <- NULL

         aux.info <- list(cov.mat=cov.mat,df=df,block.inds=
                          block.inds,random.var=random.var,
                          df.fit=df.fit,df.res=df.res)
      }

      if (!is.null(Z))
      {
         if (auto.spar.select==TRUE)
         {
            C.mat <- cbind(X,Z)
            diag.G <- diag(G)
            ridge.vec <- c(rep(0,ncol(X)),1/diag.G)

            if (!is.null(spm.info$off.set))
               ridge.reg.fit <- irls.ridge(C.mat,y,off.var=spm.info$off.set,
                                     ridge.vec=ridge.vec,
                                     max.it=50,acc=0.000001,family=family)
            else
               ridge.reg.fit <- irls.ridge(C.mat,y,ridge.vec=ridge.vec,
                                     max.it=50,acc=0.000001,family=family)
         }

         aux.info <- glmeAux(X,Z,G,block.inds,ridge.reg.fit,family)
      }
   }

   # Determine the range of fitted values for each
   # each component of the model

   if(!is.null(Z))
   {
      coef.ests <- c(lme.fit$coef$fixed,lme.fit$coef$random)
      C.mat <- cbind(X,Z)
   }
   if (is.null(Z))
   {
      coef.ests <- lme.fit$coef$fixed
      C.mat <- X
   }

   mins <- NULL
   maxs <- NULL

   for (j in 1:length(block.inds))
   {
      fitted.j <- as.matrix(C.mat[,block.inds[[j]]])%*%
                    as.matrix(coef.ests[block.inds[[j]]])

      mins[j] <- min(fitted.j)
      maxs[j] <- max(fitted.j)
   }

   aux.info <- c(aux.info,list(mins=mins,maxs=maxs))

   # Determine fitted values and residuals

   if (family=="gaussian")
   {
      fitted <- as.vector(C.mat%*%coef.ests)
      resids <- y - fitted

      lme.fit$fitted <- fitted
      lme.fit$residuals <- resids
   }

   if (family!="gaussian")
   {
      eta.hat <- C.mat%*%coef.ests

      mu.hat <- inv.link(eta.hat,family)
      fitted <- mu.hat

      if (family=="binomial")
         resids <-  binomial()$dev.resids(y,mu.hat,rep(1,length(y)))

      if (family=="poisson")
         resids <- poisson()$dev.resids(y,mu.hat,rep(1,length(y)))

      Dev <- sum(resids^2)

      # Obtain WLS approximation to deviance

      wt <- dinv.link(eta.hat,family)
      Dev.wls <- sum((y-inv.link(eta.hat,family))^2/wt)

      lme.fit <- c(lme.fit,list(fitted=fitted,residuals=resids,
                   deviance=Dev,deviance.wls=Dev.wls))
   }

   names(lme.fit$coef$fixed) <- NULL
   names(lme.fit$coef$random) <- NULL

   spm.fit.obj <- list(fit=lme.fit,info=spm.info,aux=aux.info)

   class(spm.fit.obj) <- "spm"

   return(spm.fit.obj)
}

########## End of spm ##########
