fabatch <-
function(x, y, batch, nbf=NULL, minerr=1e-06, probcrossbatch=TRUE, maxiter=100, maxnbf=12) {   
    
  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 
  if(!(is.factor(y) & all(levels(y)==(1:2))))
    stop("'y' has to be of class 'factor' with levels '1' and '2'.")
  if(max(table(batch)) >= ncol(x))
    stop("The sample size within each batch has to be smaller than the number of variables.")
  if((length(batch)!=length(y)) | (length(y)!=nrow(x)))
    stop("At least one pair out of 'x', 'y' and 'batch' is incompatible.")
	
   ##require("glmnet")

   batches = levels(batch) 
   nbatches = length(batches)
   nvars <- ncol(x)

   biggestbatch <- names(table(batch))[which.max(table(batch))]

   ylevels = levels(y)


   # If a variable is constant in a batch, infinite values
   # or missing values will arise, when scaling by batch.
   # To avoid this, the following is done:
   # Variables, for which this phenomen occurs, are left
   # out during scaling and factor estimation and only
   # adjusted with respect to their mean.
	 
   # Calculate batch-specific standard deviations:
   sdb = as.list(rep(0,nbatches))
   for (i in 1:nbatches) {
      sdb[[i]] = apply(x[batch==batches[i],],2,sd)
   }
   sdb0=sdb

   badvariableslist <- lapply(sdb, function(x) which(x==0))
   badvariables <- sort(unique(unlist(badvariableslist)))
   goodvariables <- setdiff(1:ncol(x), badvariables)

   if(length(badvariables) > 0) {
     xsafe <- x
     x <- x[,goodvariables]
     sdb <- lapply(sdb, function(x) x[goodvariables])
     sdb0=sdb
   }

   # overall variable means:
   
   meanoverall <- colMeans(x)
   
   # First center and scale the variables per batch to apply the Lasso
   # for obtaining the probability estimates:
    
   # Remove batch specific means:
   if(nbatches > 1) {
     adjustmentmod = lm(x~batch)
     design = model.matrix(~batch)
     adjustmentcoef = coef(adjustmentmod)
     xb = x-design%*%adjustmentcoef
     adjustmentcoef0 = adjustmentcoef
   } else
     xb <- scale(x, scale=FALSE)  

   # Pooled variance:

   pooledsds <- apply(xb, 2, sd)   


   # 'scaledxb' is X centered and scaled per batch:    
   scaledxb = xb
   for (i in 1:nbatches) {
         scaledxb[batch==batches[i],] = scale(xb[batch==batches[i],],center=rep(0,ncol(x)),scale=sdb[[i]])
   }

   # Now lasso is trained on the batch-centered and batch-scaled X and the fitted probability
   # of each y to belong to class 1 or 2 calculated - 'proba' (There is probably a problem with overfitting:
   # When the probabilities are much more shifted towards 0 and 1 in the training data than in the 
   # test data, the training and test data might not be that comparable anymore -> Solution might be
   # to use cross validation to estimate the probabilities here):   
   cvmod = glmnet::cv.glmnet(x=as.matrix(scaledxb),y=y,family="binomial",type.measure="deviance", alpha = 0)
   lambda.min = cvmod$lambda.min
   mod = glmnet::glmnet(x=as.matrix(scaledxb),y=y,family="binomial",lambda=lambda.min, alpha = 0)


   if(probcrossbatch) {

   if(nbatches > 1) {
     proba <- rep(NA, nrow(x))
     for(i in 1:nbatches) {
       modtemp = glmnet::glmnet(x=as.matrix(scaledxb[batch!=batches[i],]),y=y[batch!=batches[i]],family="binomial",lambda=lambda.min, alpha = 0)
       b0temp = modtemp$a0
       btemp = modtemp$b  
	   proba[batch==batches[i]] = as.vector(1/(1+exp(-b0temp-as.matrix(scaledxb[batch==batches[i],])%*%btemp))) # mean(ifelse(p1x>0.5,2,1)!=ytest)
	 }   
   } else {
     warning("Number of Training batches equal to one. Using ordinary cross-validation for preliminary class probability estimation.")
     cvmod = glmnet::cv.glmnet(x=as.matrix(scaledxb),y=y,family="binomial",type.measure="deviance", keep=TRUE, lambda=c(lambda.min, lambda.min+(1/1000)*lambda.min), alpha = 0)
     proba = cvmod$fit.preval[,1]
   }   

  }
  else {
    cvmod = glmnet::cv.glmnet(x=as.matrix(scaledxb),y=y,family="binomial",type.measure="deviance", keep=TRUE, lambda=c(lambda.min, lambda.min+(1/1000)*lambda.min), alpha = 0)
    proba = cvmod$fit.preval[,1]
  }      

   # 'b0' is intercept from Lasso and 'b' coefficient vector for variables:
   b0 = mod$a0
   b = mod$b


   # Calculate the class-specific mean vectors 'm1' and 'm2':
   m1 = apply(scaledxb[y==ylevels[1],],2,mean)
   m2 = apply(scaledxb[y==ylevels[2],],2,mean) 


   nbfvec <- rep(NA, nbatches)

   # Determine number of factors if not given:    
   if (!is.null(nbf)) {
     nbfvec <- rep(nbf, times=nbatches)
     nbfinput <- NULL
   }
   else
     nbfinput <- nbf

   # Calculate the factors on 'cdta' (batch-centered-scaled and class-removed):


   falist <- vector("list", nbatches)
   
   criterionall <- list()

   scaledxbfa <- scaledxb

   for (i in 1:nbatches) {
   
       if(is.null(nbf)) {
           ##require("mnormt")
           maxnbf2 <- min(c(maxnbf, floor(sum(batch==batches[i])/2)))
           nbfobj = nbfactors(scale(sweep(scale(scaledxb[batch==batches[i],],center=m1,scale=FALSE), 1, 1-proba[batch==batches[i]], "*") + 
              sweep(scale(scaledxb[batch==batches[i],],center=m2,scale=FALSE), 1, proba[batch==batches[i]], "*")), maxnbfactors=maxnbf2, minerr=minerr, maxiter=maxiter)
           nbfvec[i] <- nbfobj$optimalnbfactors
       }

       if(is.na(nbfvec[i])) {
         warning("There occured an issue in the factor estimation. Number of factors set to zero.")
         nbfvec[i] <- 0
       }   
    
	 if(nbfvec[i] > 0) {
        # Calculate the factors on 'cdta' (batch-centered-scaled and class-removed):
        fa = emfahighdim(sweep(scale(scaledxb[batch==batches[i],],center=m1,scale=FALSE), 1, 1-proba[batch==batches[i]], "*") + 
          sweep(scale(scaledxb[batch==batches[i],],center=m2,scale=FALSE), 1, proba[batch==batches[i]], "*"),nbf=nbfvec[i],minerr=minerr, maxiter=maxiter)

        # Remove the factor influences:
        scaledxbfa[batch==batches[i],] = scaledxb[batch==batches[i],] - fa$Factors%*%t(fa$B)

        falist[[i]] <- fa
		
     } else {
       scaledxbfa[batch==batches[i],] <- scaledxb[batch==batches[i],]
       fa <- NULL
     }

   } 
   
   means2batch <- sd2batch <- matrix(nrow=length(levels(batch)), ncol=ncol(scaledxbfa))

   # scale again:    
   for (i in 1:nbatches) {
         means2batch[i,] <- colMeans(scaledxbfa[batch==batches[i],])
         sd2batch[i,] <- apply(scaledxbfa[batch==batches[i],], 2, sd)
         scaledxbfa[batch==batches[i],] = scale(scaledxbfa[batch==batches[i],], center=means2batch[i,], scale=sd2batch[i,])
   }
   
   xfa <- sweep(sweep(scaledxbfa, 2, pooledsds, "*"), 2, meanoverall, "+")
   
   if(length(badvariables)>0) {
     xbadvar <- xsafe[,badvariables, drop=FALSE]
	 
	 sdbbad = as.list(rep(0,nbatches))
     for (i in 1:nbatches) {
       sdbbad[[i]] = apply(xbadvar[batch==batches[i],],2,sd)
     }
     whichzerosdmat <- t(sapply(sdbbad, function(x) x==0))

     sdbbad0to1 <- lapply(sdbbad, function(x) {
       x[x==0] <- 1
       x
     })
	 
     pooledsdsbad <- mapply(function(x, y) {
      if(!all(y)) {
        return(sd(x[batch %in% which(!y)]))
      }
      else
        return(1)
     }, as.data.frame(xbadvar), as.data.frame(whichzerosdmat))
	 	 
     mubadvar <- colMeans(xbadvar)
     xbadvaradj <- matrix(nrow=nrow(xsafe), ncol=length(badvariables), data=mubadvar, byrow=TRUE)
     for(i in 1:nbatches)
       xbadvaradj[batch==i,] <- xbadvaradj[batch==i,,drop=FALSE] + sweep(scale(xbadvar[batch==i,,drop=FALSE], center=TRUE, scale=sdbbad0to1[[i]]), 2, pooledsdsbad, "*")
     xfanew <- matrix(nrow=nrow(xsafe), ncol=ncol(xsafe))
     xfanew[,goodvariables] <- xfa
     xfanew[,badvariables] <- xbadvaradj
     xfa <- xfanew
   }

     meanvectorwhole <- rep(NA, nvars)
	 meanvectorwhole[goodvariables] <- meanoverall
	 if(length(badvariables)>0)
	   meanvectorwhole[badvariables] <- mubadvar

	 sdvectorwhole <- rep(NA, nvars)
	 sdvectorwhole[goodvariables] <- pooledsds
	 if(length(badvariables)>0)
	   sdvectorwhole[badvariables] <- 1
	   
   params <- list(xadj=xfa, m1=m1,m2=m2,b0=b0,b=b, pooledsds=sdvectorwhole, meanoverall=meanvectorwhole, minerr=minerr,
     nbfinput=nbfinput, badvariables=badvariables, nbatches=nbatches, batch=batch, nbfvec=nbfvec)
	 
  class(params) <- "fabatch"
	 
  return(params)
	 
}
