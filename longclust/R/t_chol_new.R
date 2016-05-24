longclustEM <- function(x, Gmin, Gmax, class = NULL, linearMeans = FALSE, modelSubset = NULL, initWithKMeans = FALSE, criteria = "BIC", equalDF = FALSE, gaussian = FALSE,userseed=1004) {
	if(criteria != "BIC" && criteria != "ICL")
	  stop("criteria should be \"BIC\" or \"ICL\" ");
  is.int <- function(no) {
	  abs(no - round(no)) < 1e-15
  }
	if(!is.int(Gmin) || !is.int(Gmax) || Gmin > Gmax || Gmin < 2)
	  stop("Gmin and Gmax should be integers, with 2 <= Gmin <= Gmax");
  
	if(!is.logical(linearMeans))
	  stop("linearMeans should be TRUE or FALSE");
	if(!is.logical(equalDF))
	  stop("equalDF should be TRUE or FALSE");
  if(!is.logical(gaussian))
    stop("gaussian should be TRUE or FALSE");
	if(!is.logical(initWithKMeans))
	  stop("initWithKMeans should be TRUE or FALSE");

  df_flag <- 0;
	if(equalDF)
	 df_flag <- 1;

  gaussian_flag <-0; 
  if(gaussian)
   gaussian_flag <- 1;

	models_all <- c("VVI", "VVA", "EEI", "EEA", "VEI", "VEA", "EVA", "EVI");
	models_selected <- c(0,0,0,0,0,0,0,0);
  if(is.null(modelSubset))
	  modelSubset = models_all;
  for(mod in modelSubset) {
	  found = FALSE;
  	for(i in 1:8) {
	    if(mod == models_all[i]) {
			  found = TRUE;
				models_selected[i] = 1;
				break;
			}
	  }
		if(found == FALSE)
		  stop("modelSubset should be subset of {\"VVI\", \"VVA\", \"EEI\", \"EEA\", \"VEI\", \"VEA\", \"EVA\", \"EVI\"} ");
	}

	as.matrix(x);
  N <- dim(x)[1];
  p <- dim(x)[2];
#  if(!is.loaded("tcholmeans_entry") || !is.loaded("tchol_entry"))
#    dyn.load('t_chol.so'); 
  set.seed(userseed);

  use_kmeans <- 0;
	kmeansz <- NULL;
	if(initWithKMeans) {
	  if(!is.null(class)) {
		  print("Warning: Ignoring the option of initializing with k-means, since class is non-empty.");
		} else {
		  use_kmeans <- 1; 
      kmeansz <- matrix(0, N,Gmax+1-Gmin);
      for(G in Gmin:Gmax) 
  	    kmeansz[,G+1-Gmin] <- kmeans(x, G, nstart=10)$cluster;
		}
  }

  if(!is.null(class)) {
    if(any(!is.int(class)) || min(class) < 0) 
		  stop("The vector class may contain non-negative integers only.");
		if(length(class) != N) 
		  stop("The vector class must have length equal to the number of samples.");
		if(max(class) > Gmin)
		  stop("The parameter Gmin cannot be less than max(class).");
  } else{
	  class <- rep(0, N);
	}

  if(linearMeans == TRUE) {
  	time1 <- matrix(0,p,2);
    time1[,1] <- rep(1,p);
    time1[,2] <- 1:p;
		
    res <- .C("tcholmeans_entry", as.double(x), as.integer(class), as.double(time1), as.integer(N), as.integer(p), as.integer(Gmin), as.integer(Gmax), as.integer(df_flag), as.integer(gaussian_flag), as.integer(models_selected), as.integer(use_kmeans), as.double(kmeansz), Gbest = as.integer(0), Gbest_icl = as.integer(0), zbest = double(N*(Gmax+1)),zbest_icl = double(N*(Gmax+1)), nubest = double(Gmax+1), nubest_icl = double(Gmax+1), mubest = double((Gmax+1)*p), mubest_icl = double((Gmax+1)*p), Tbest = double((Gmax+1)*p*p), Tbest_icl =double((Gmax+1)*p*p), Dbest = double((Gmax+1)*p*p), Dbest_icl = double((Gmax+1)*p*p), llres = double(8*(Gmax+1)), bicres = double(8*(Gmax+1)), iclres = double(8*(Gmax+1)), package="lonclust");
	} else {
    res <- .C("tchol_entry", as.double(x), as.integer(class), as.integer(N), as.integer(p), as.integer(Gmin), as.integer(Gmax), as.integer(df_flag), as.integer(gaussian_flag), as.integer(models_selected), as.integer(use_kmeans), as.double(kmeansz), Gbest = as.integer(0), Gbest_icl = as.integer(0), zbest = double(N*(Gmax+1)), zbest_icl = double(N*(Gmax+1)), nubest = double(Gmax+1), nubest_icl = double(Gmax+1), mubest = double((Gmax+1)*p), mubest_icl = double((Gmax+1)*p), Tbest = double((Gmax+1)*p*p), Tbest_icl = double((Gmax+1)*p*p), Dbest = double((Gmax+1)*p*p), Dbest_icl = double((Gmax+1)*p*p), llres = double(8*(Gmax+1)), bicres = double(8*(Gmax+1)), iclres = double(8*(Gmax+1)), package="longclust");
	}
  
  if(criteria == "BIC") {
    Gbest <- res[['Gbest']];
    nubest <- res[['nubest']][1:Gbest];
	  mubest <- matrix(res[['mubest']][1:(Gbest*p)], Gbest, p);
	  temp <- matrix(res[['zbest']], N, Gmax+1);
	}	else {
	  Gbest <- res[['Gbest_icl']];
    nubest <- res[['nubest_icl']][1:Gbest];
	  mubest <- matrix(res[['mubest_icl']][1:(Gbest*p)], Gbest, p);
	  temp <- matrix(res[['zbest_icl']], N, Gmax+1);
	}
	zbest <- temp[, 1:Gbest];

  Tbest <- list();
	Dbest <- list();
	prev = 1;
	if(Gbest >= 2) {
	for(g in 1:Gbest){
	  if(criteria == "BIC") {
  	  Tbest[[g]] <- matrix(res[['Tbest']][prev:(g*p*p)], p, p);
	  	Dbest[[g]] <- matrix(res[['Dbest']][prev:(g*p*p)], p, p);
	  } else {
  	  Tbest[[g]] <- matrix(res[['Tbest_icl']][prev:(g*p*p)], p, p);
	  	Dbest[[g]] <- matrix(res[['Dbest_icl']][prev:(g*p*p)], p, p);
		}
		prev <- g*p*p+1;
	}
	}

  rows <- c(Gmin:Gmax);
	models <- c("VVI", "VVA", "EEI", "EEA", "VEI", "VEA", "EVA", "EVI");
  model_indices <- vector();
  for(i in 1:8) {
	  if(models_selected[i] == 1)
		  model_indices <- append(model_indices, i);
	}

	temp <- matrix(res[['llres']], 8, Gmax+1);
	llres <- t(temp[, (Gmin+1):(Gmax+1)]);

	for(i in 1:(Gmax+1-Gmin)) {
    for(j in 1:8) {
		  if(llres[i,j] < -(2^27))
			  llres[i,j] <- -Inf; 
			if(llres[i,j] > (2^27)) 
			  llres[i,j] <- Inf; 
		} 
	}
	colnames(llres) <- models;
	rownames(llres) <- rows;
	llres <- llres[,model_indices];

	temp <- matrix(res[['bicres']], 8, Gmax+1);
	bicres <- t(temp[, (Gmin+1):(Gmax+1)]);
	for(i in 1:(Gmax+1-Gmin)) {
    for(j in 1:8) {
		  if(bicres[i,j] < -(2^27))
			  bicres[i,j] <- -Inf;
			if(bicres[i,j] > (2^27))
			  bicres[i,j] <- Inf;
		}
	}
  colnames(bicres) <- models;
	rownames(bicres) <- rows;
	bicres <- bicres[,model_indices];


	temp <- matrix(res[['iclres']], 8, Gmax+1);
	iclres <- t(temp[, (Gmin+1):(Gmax+1)]);
	for(i in 1:(Gmax+1-Gmin)) {
    for(j in 1:8) {
		  if(iclres[i,j] < -(2^27))
			  iclres[i,j] <- -Inf;
			if(iclres[i,j] > (2^27))
			  iclres[i,j] <- Inf;
		}
	}
  colnames(iclres) <- models;
	rownames(iclres) <- rows;
	iclres <- iclres[,model_indices];

	result <- list(Gbest=Gbest, zbest=zbest, nubest=nubest, mubest = mubest, Tbest=Tbest, Dbest = Dbest, llres=llres, bicres=bicres, iclres=iclres, gaussian=gaussian);
  class(result) <- "longclust"
  result

}

summary.longclust<- function(object, ...) {
  if(object$Gbest >= 2) {
    cat("Number of components:", object$Gbest, "\n");
    cat("z:\n");
    print(object$zbest);
	  cat("\n");
	  for(g in 1:object$Gbest) {
      cat("Cluster: ", g, "\n");
			if(object$gaussian == FALSE)
	      cat("v: ", object$nubest[g], "\n");
	    cat("mean:", object$mubest[g,], "\n");
  	  cat("T:\n");
	    print(object$Tbest[[g]]);
	    cat("D:\n");
  	  print(object$Dbest[[g]]);
	    cat("\n");
    }
	}
	cat("bicres:\n");
  printCoefmat(object$bicres);
  cat("iclres:\n");
  printCoefmat(object$iclres);
}

print.longclust<- function(x, ...) {
  if(x$Gbest >= 2) {
    cat("Number of components:", x$Gbest, "\n");
    cat("z:\n");
    print(x$zbest);
  	cat("\n");
  	for(g in 1:x$Gbest) {
      cat("Cluster: ", g, "\n");
			if(x$gaussian == FALSE)
	      cat("v: ", x$nubest[g], "\n");
	    cat("mean:", x$mubest[g,], "\n\n");
    }
	} else{
	  print("None of the model-component combinations have a finite score.");
	}
}

plot.longclust<- function(x, data, ...) {
  if(x$Gbest < 2)
	  return();
	as.matrix(data);
	if(dim(x$zbest)[1] != dim(data)[1])
	  stop("The data provided is incompatible with the longxt object");
	N <- dim(data)[1];
  mp <- rep(0, N);
  mump <- rep(1: (x$Gbest));
	for(i in 1:N) {
	  max <- 0;
		bestj <- 1;
	  for(j in 1:(x$Gbest)){
		  if(max < x$zbest[i,j]) {
			  max = x$zbest[i,j];
				bestj <- j;
			}
		}
		mp[i] <- bestj;
	}
  par(mfrow=c(1,1),ask=FALSE)
  matplot(t(data), type="l", lty=2, col=mp, xlab = "Time", ylab = "Values", main = "Components");
  matplot(t(x$mu), type="l", lty=1, lwd = 4, col=mump, add=TRUE);
  par(mfrow=c(1,2),ask=TRUE)

	p <- 2;
	q <- 1;
	flag <- 1;
	while( (p*q) < (x$Gbest) ) {
	  if(flag == 1) {
		  q <- q + 1;
		} else {
		  p <- p + 1;
		}
		flag = 1 - flag;
	}
	par(mfrow = c(p,q));
	for(i in 1:(x$Gbest)) {
	  ind <- which(mp == i);
    if(length(ind) == 0) {
      cat("Cluster ", i, " is empty.");
      next;
    }
		str <- sprintf("Component %d", i);  
		matplot(t(data)[,ind], type="l", lty=1, col=mp[ind], xlab = "Time", ylab = "Values", main=str);
    matplot(x$mu[i,], type="l", lty=1, lwd=4, col=mp[ind][1], add=TRUE);
	}
}
