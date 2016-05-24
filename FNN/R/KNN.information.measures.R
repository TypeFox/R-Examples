################################################################################
# KNN Information Theory Measures                                              #
# File:   KNN.information.measures.R                                           #
# Author: Shengqiao Li                                                         #
# Date:   December 12, 2008                                                    #
#         2010-3-25 add entropy and crossentropy                               #
################################################################################
entropy<- function (X, k=10, algorithm=c("kd_tree", "brute"))
{
  algorithm<- match.arg(algorithm);

  #check data
  if (!is.numeric(X))  stop("Data non-numeric")
  if (any(is.na(X)))   stop("Data include NAs")
  if (!is.matrix(X))   X <- matrix(X)
  if (storage.mode(X) == "integer")  storage.mode(X) <- "double"

  n <- nrow(X)
  p <- ncol(X)
  if (k >= n) stop("k must less than the sample size!")

  Cname<- switch(algorithm,
              kd_tree = "KNN_MLD_kd",
              brute = "KNN_MLD_brute"
  );

  MLD<- .C(Cname, t(X), as.integer(k), p, n, MLD = double(k), DUP = FALSE)$MLD
  # mean of log dist
  
  #digamma(n) is more popular than log(n)  
  #H <- log(n) - digamma(1:k) + p/2*log(pi) - lgamma(p/2+1) + p*MLD 
  H <- digamma(n) - digamma(1:k) + p/2*log(pi) - lgamma(p/2+1) + p*MLD 
  
  return(H)
}

crossentropy <- function(X, Y, k=10, algorithm=c("kd_tree", "cover_tree", "brute"))
{
  algorithm<- match.arg(algorithm);
	if (!(is.numeric(X)&& is.numeric(Y))) stop("Data non-numeric");
	if (any(is.na(X), is.na(Y)))  stop("Data include NAs");

	if (!is.matrix(X))  X<- matrix(X);
	if (!is.matrix(Y))  Y<- matrix(Y);


	n <- nrow(X); m<- nrow(Y);
	p <- ncol(X);

	if(k>=n) stop("k must less than the sample size!");

	dnn<- knnx.dist(Y, X, k=k, algorithm=algorithm);
  MLD<- colMeans(log(dnn));
	#H<- p*MLD + p/2*log(pi) - lgamma(p/2+1) + log(m) - digamma(1:k)
  H <- digamma(m) - digamma(1:k) + p/2*log(pi) - lgamma(p/2+1) + p*MLD 

	return (H);
}

KL.divergence<- function(X, Y, k=10, algorithm=c("kd_tree", "cover_tree", "brute"))
{
  #Kullback-Leibler Distance
  algorithm<- match.arg(algorithm);

  if (!is.matrix(X))  X<- matrix(X);
  if (!is.matrix(Y))  Y<- matrix(Y);

  n<- nrow(X); p<- ncol(X);
  m<- nrow(Y);

  log(m/n) + p*(colMeans(log(knnx.dist(Y, X, k=k, algorithm)))- colMeans(log(knn.dist(X, k=k, algorithm))));

}

KL.dist<- function(X, Y, k=10,  algorithm=c("kd_tree", "cover_tree", "brute"))
{
  #Symmetric Kullback-Leibler divergence. i.e. Kullback-Leibler distance
  algorithm<- match.arg(algorithm);

  KL.divergence(X, Y, k, algorithm) + KL.divergence(Y, X, k, algorithm)

}

KLx.divergence<- function (X, Y, k = 10, algorithm="kd_tree")
{
  #Kullback-Leibler divergence
  algorithm<- match.arg(algorithm);
  if(storage.mode(X)=="integer") storage.mode(X)<- "double";
  if(storage.mode(Y)=="integer") storage.mode(Y)<- "double";

  if(!is.matrix(X)) X<- as.matrix(X);
  if(!is.matrix(Y)) Y<- as.matrix(Y);

  n<- nrow(X); m<- nrow(Y);
  d <- ncol(X); p<- ncol(Y);

  if(d!=p) stop("Number of columns must be same!.");
  if(k>=n) warning("k should be less than sample size!");

  .C("KL_divergence", t(X), t(Y), as.integer(k), d, n, m, KL = double(k), DUP=FALSE)$KL;

}

KLx.dist<- function (X, Y, k = 10, algorithm="kd_tree")
{
  #Symmetric Kullback-Leibler divergence. i.e. Kullback-Leibler distance
  algorithm<- match.arg(algorithm);
  if(storage.mode(X)=="integer") storage.mode(X)<- "double";
  if(storage.mode(Y)=="integer") storage.mode(Y)<- "double";

  if(!is.matrix(X)) X<- as.matrix(X);
  if(!is.matrix(Y)) Y<- as.matrix(Y);

  n<- nrow(X); m<- nrow(Y);
  d<- ncol(X); p<- ncol(Y);

  if(d!=p) stop("Number of columns must be same!.");
  if(k>=n) warning("k should be less than sample size!");

  .C("KL_dist", t(X), t(Y), as.integer(k), d, n, m, KLD = double(k), DUP=FALSE)$KLD;

}


mutinfo.R<- function(X, Y, k = 10, direct=TRUE)
{
  #mutual information of X and Y 
  #KSG method by Kraskov, Stogbauer and Grassberger
  #slow. for test only
  
  if(is.vector(X)+is.vector(Y) < 2) stop("X and Y must be vectors.");
  
  n<- length(X);

  nx<- ny<- integer(n);
  di<- double(n);
  
  for (i in 1:n){
    dx<- abs(X[i]- X);
    dy<- abs(Y[i]- Y);    
    r<- ifelse (dx > dy, dx, dy);

    di[i]<- max(r[rank(r)<= (k+1)]); #exclude self-match

    nx[i]<- sum(dx < di[i]); #include self-match
    ny[i]<- sum(dy < di[i]); #include self-match
  }      
  
  mi<- digamma(n) + digamma(k) - mean(digamma(nx) + digamma(ny));

  return(mi);
}


mutinfo<- function(X, Y, k = 10, direct=TRUE)
{
  #mutual information of X and Y
  
  #check data
  if (!is.numeric(X))  stop("Data non-numeric")
  if (any(is.na(X)))   stop("Data include NAs")  
  if (!is.matrix(X))   X <- matrix(X)
  if (storage.mode(X) == "integer")  storage.mode(X) <- "double"
    
  if (!is.numeric(Y))  stop("Data non-numeric")
  if (any(is.na(Y)))   stop("Data include NAs")
  if (!is.matrix(Y))   Y <- matrix(Y)
  if (storage.mode(Y) == "integer")  storage.mode(Y) <- "double"

  
  n <- nrow(X)

  if (k >= n) stop("k must less than the sample size!")

  p1<- ncol(X);
  p2<- ncol(Y);

  
  if (direct){
    #KSG method by Kraskov, Stogbauer and Grassberger
    
    #res<- .C("mutinfo", rbind(X, Y), as.integer(k), n, nx = integer(n), ny = integer(n), DUP=FALSE)    
    res<- .C("mdmutinfo", t(X), t(Y), as.integer(p1), as.integer(p2), as.integer(k), n, nx = integer(n), ny = integer(n), DUP=FALSE)    
    nx<- res$nx;
    ny<- res$ny; 

    mi<- digamma(n) + digamma(k) - mean(digamma(nx) + digamma(ny));

  }
  else {
    algorithm = "kd_tree";
    Hx<- entropy(X, k, algorithm);
    Hy<- entropy(Y, k, algorithm)
    H<-  entropy (cbind(X, Y), k, algorithm);
    mi<- Hx + Hy - H;
  }  
  
  return(mi);
}
