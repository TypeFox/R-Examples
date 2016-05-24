gnsc.train <-
function(x, col.struc=NULL, row.struc=NULL, standardize = T, nlambda = NULL, lambda.max = 10, lambda = NULL, verbose = TRUE){
	  gcinfo(FALSE)
	  d = nrow(x);
	  n = ncol(x);
	  fit = list()
	  if(is.null(col.struc)) col.struc = 1:n
    if(is.null(row.struc)) row.struc = 1:d
    if(length(col.struc)!=n)  stop("col.struc must have the same length as the column of the input data.")
    if(length(row.struc)!=d)  stop("row.struc must have the same length as the row of the input data.")    
    x=x[order(row.struc),]; row.struc=row.struc[order(row.struc)]    

	  if(!is.null(lambda)) nlambda = length(lambda)
	  if(is.null(lambda)){
		   if(is.null(nlambda)) nlambda = 10
		   lambda.min = 0.1*lambda.max
		   lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
		   rm(lambda.max,lambda.min)
		   gc()
    }
    
    ## summary statistics for row and column class
    row.class = unique(row.struc); K = length(row.class) 
    col.class = unique(col.struc); M = length(col.class)
    if(K < 2){ cat("no. of variables should be higher than 1.\n"); return()} 
    if(M < 2){ cat("no. of sample classes should be higher than 1.\n"); return()}    
    ## standardize the data
    if(standardize)
       x = (x-apply(x,1,mean))/apply(x,1,sd)
    ## calculate the hatmu
  	mu = matrix(0, nrow=d[1],ncol=M)
	  for(ic in 1:M)
		    mu[,ic] = rowMeans(matrix(x[, col.struc == col.class[ic]],ncol=sum(col.struc == col.class[ic])))
	  all.mean = 0
    if(!standardize) all.mean=apply(x,1,mean)
    ## calculate the omega term
    omega = sqrt(sapply(1:K,function(k) sum(row.struc == row.class[k]))%*%t(sapply(1:M,function(m) 1/sum(col.struc == col.class[m]))))

    ## perform the shrunkening
	  tilde.mu = array(0, dim=c(d,M,nlambda))
	  relate.mat = array(0, dim=c(K,M,nlambda))
	  Thresh.mat = array(0, dim=c(K,M,nlambda))
	  
	  if(verbose)  cat("Conducting Group Nearest Shrunken Centroids...\n")
  	clust.se.id =  gnsc.restruc(col.struc)$path.se.id
  	reorder.path = gnsc.restruc(row.struc); path.se.id = reorder.path$path.se.id; path.n = reorder.path$path.n
    inv.cov.mat = gnsc.icov(x, pathway=row.struc, path.n=path.n, path.se.id=path.se.id)
    
    for(ci in 1:M){ 
		   cimat=x[,clust.se.id[ci,1]:clust.se.id[ci,2]]
		   hat_mu.all = (mu[,ci] - all.mean)		
		   for(k in 1:K){
			     tmp.id = (path.se.id[k,1]:path.se.id[k,2])
			     hat_mu = hat_mu.all[tmp.id]
			     tmp = sqrt(hat_mu%*%inv.cov.mat[[k]]%*%hat_mu)
			     for(sim in 1:nlambda){
               Thresh.mat[k,ci,sim] = 1-lambda[sim]*omega[k,ci]/tmp
               thresh = max(Thresh.mat[k,ci,sim], 0)
               if(thresh != 0 ) relate.mat[k,ci,sim] = 1
               tilde.mu[tmp.id,ci,sim] = thresh * hat_mu
           }
	     }
		   if(verbose)  cat("Conducting Group Nearest Shrunken Centroids...", floor(100*ci/M),"%", "\r")
   }  
   
   nonzero = NULL
   for(sim in 1:nlambda)
      nonzero = c(nonzero, sum(relate.mat[,,sim])) 
   
	 colnames(tilde.mu) = col.class
	 rownames(tilde.mu) = row.struc
	 colnames(relate.mat) = col.class
	 rownames(relate.mat) = row.class
	 colnames(Thresh.mat) = col.class
	 rownames(Thresh.mat) = row.class         
   
   yhat = matrix(0,n,nlambda); errors=rep(0,nlambda)
   for(sim in 1:nlambda){
       out = gnsc.predict(x, x.class=col.struc, all.mean, tilde.mu[,,sim], row.struc, col.struc, inv.cov.mat, path.se.id)
       yhat[,sim] = out$yhat
       errors[sim] = out$error
       if(verbose) cat(floor(100*sim/nlambda),"%")
   }
   colnames(yhat) = lambda
   
   fit$tilde.mu = tilde.mu
   fit$relate.mat = relate.mat
   fit$Thresh.mat = Thresh.mat
   fit$lambda = lambda
   fit$nlambda = nlambda
   fit$yhat = yhat
   fit$errors = errors 
   fit$icov = inv.cov.mat
   fit$nonzero = nonzero
   fit$all.mean = all.mean   
   
   fit$row.struc = row.struc
   fit$col.struc = col.struc
   fit$path.se.id = path.se.id
   fit$standardize = standardize
   rm(tilde.mu,relate.mat,Thresh.mat,inv.cov.mat,nonzero)
      
   if(verbose) cat("done\n")
   class(fit) = "gnsc"
   return(fit)
}
