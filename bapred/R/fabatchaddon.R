fabatchaddon <-
function(params, x, batch) {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 
  if(max(table(batch)) >= ncol(x))
    stop("The sample size within each batch has to be smaller than the number of variables.")
	
   if(class(params) != "fabatch")
     stop("Input parameter 'params' has to be of class 'fabatch'.")
	 	 
   if(ncol(params$xadj) != ncol(x))
     stop("Number of variables in test data matrix different to that of training data matrix.")	 

   batches = levels(batch)
   nbatches = length(batches)
   nvars <- ncol(x)

   sdb = as.list(rep(0,nbatches))
   for (i in 1:nbatches) {
      sdb[[i]] = apply(x[batch==batches[i],],2,sd)
   }
   sdb0=sdb
      
   badvariableslisthere <- lapply(sdb, function(x) which(x==0))
   badvariableshere <- sort(unique(unlist(badvariableslisthere)))
   goodvariableshere <- setdiff(1:nvars, badvariableshere)
    	  
   if(nbatches > 1) {
     # Remove batch specific means:
     adjustmentmod = lm(x~batch)
     design = model.matrix(~batch)
     adjustmentcoef = coef(adjustmentmod)
     xb = x-design%*%adjustmentcoef
     adjustmentcoef0 = adjustmentcoef
   }
   else
     xb = scale(x, scale=FALSE)
	 
   pooledsds <- apply(xb, 2, sd)   
   
   # 'scaledxb' is X centered and scaled per batch:    
   scaledxb = xb
   for (i in 1:nbatches) {
     scaledxb[batch==batches[i],setdiff(1:nvars, badvariableslisthere[[i]])] = 
	   scale(xb[batch==batches[i],setdiff(1:nvars, badvariableslisthere[[i]])],center=rep(0,nvars),scale=sdb[[i]][setdiff(1:nvars, badvariableslisthere[[i]])])
   }   
   
   # Predict probabilities using the model estimated on the training data:
   p1x = as.vector(1/(1+exp(-params$b0-as.matrix(scaledxb)[,setdiff(1:nvars, params$badvariables)]%*%params$b))) # mean(ifelse(p1x>0.5,2,1)!=ytest)

  
  herealsogood <- lapply(badvariableslisthere, 
    function(x) which(setdiff(1:nvars, params$badvariables) %in% setdiff(1:nvars, x)))
  goodinboth <- lapply(badvariableslisthere,
    function(x) intersect(setdiff(1:nvars, params$badvariables), setdiff(1:nvars, x)))

	
   scaledxbfa <- scaledxb
   
   for (i in 1:nbatches) {
   
       # Determine number of factors if not given:    
       if (is.null(params$nbfinput)){
         ##require("mnormt")
           nbf = nbfactors(scale(sweep(scale(scaledxb[batch==batches[i],goodinboth[[i]]],center=params$m1[herealsogood[[i]]],scale=FALSE), 1, 1-p1x[batch==batches[i]], "*") + 
              sweep(scale(scaledxb[batch==batches[i],goodinboth[[i]]],center=params$m2[herealsogood[[i]]],scale=FALSE), 1, p1x[batch==batches[i]], "*")), maxnbfactors=min(c(floor(sum(batch==batches[i])/2), 12)), minerr=params$minerr)$optimalnbfactors
       }
       else
          nbf <- params$nbfinput  
	   
       if(is.na(nbf)) {
         warning("There occured an issue in the factor estimation. Number of factors set to zero.")
         nbf <- 0
       }   
    
	  if(nbf > 0) {
        # Calculate the factors on 'cdta' (batch-centered-scaled and class-removed):
        fa = emfahighdim(sweep(scale(scaledxb[batch==batches[i],goodinboth[[i]]],center=params$m1[herealsogood[[i]]],scale=FALSE), 1, 1-p1x[batch==batches[i]], "*") + 
          sweep(scale(scaledxb[batch==batches[i],goodinboth[[i]]],center=params$m2[herealsogood[[i]]],scale=FALSE), 1, p1x[batch==batches[i]], "*"),nbf=nbf,minerr=params$minerr)

        # Remove the factor influences:
        scaledxbfa[batch==batches[i],goodinboth[[i]]] = scaledxb[batch==batches[i],goodinboth[[i]]] - fa$Factors%*%t(fa$B)
      }
      else {
        scaledxbfa[batch==batches[i],goodinboth[[i]]] <- scaledxb[batch==batches[i],goodinboth[[i]]]
        fa <- NULL
      }

   } 
   
   means2batch <- sd2batch <- matrix(nrow=length(levels(batch)), ncol=ncol(scaledxbfa))

   # scale again:    
   for (i in 1:nbatches) {
         means2batch[i,] <- colMeans(scaledxbfa[batch==batches[i],])
         sd2batch[i,] <- apply(scaledxbfa[batch==batches[i],], 2, sd)
		 sd2batch[i,][sd2batch[i,]==0] <- 1
         scaledxbfa[batch==batches[i],] = scale(scaledxbfa[batch==batches[i],], center=means2batch[i,], scale=sd2batch[i,])
   }
   
   xadj <- sweep(sweep(scaledxbfa, 2, params$pooledsds, "*"), 2, params$meanoverall, "+")

   if(length(params$badvariables) > 0) {
     for(i in 1:nbatches) {
	   if(length(setdiff(params$badvariables, badvariableslisthere[[i]])) > 0) {
    	 heregood <- setdiff(params$badvariables, badvariableslisthere[[i]])
	     meanssave <- colMeans(xadj[batch==batches[i],heregood, drop=FALSE])
		 xadj[batch==batches[i],heregood, drop=FALSE] <- scale(xadj[batch==batches[i],heregood, drop=FALSE], center=TRUE, scale=TRUE)
         xadj[batch==batches[i],heregood, drop=FALSE] <- sweep(sweep(xadj[batch==batches[i],heregood, drop=FALSE], 2, pooledsds[heregood], "*"), 2, meanssave, "+")  		 
       }
     }	 
   }
   
   return(xadj)
   
}
