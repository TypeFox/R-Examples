NbClust <-function(data = NULL, diss=NULL, distance ="euclidean", min.nc=2, max.nc=15, method =NULL, index = "all", alphaBeale = 0.1)
{
    
    x<-0
    min_nc <- min.nc
    max_nc <- max.nc
    
    if(is.null(method))    
      stop("method is NULL")
    method <- pmatch(method, c("ward.D2", "single", "complete", "average", 
                               "mcquitty", "median", "centroid", "kmeans","ward.D"))
    
        
    indice <- pmatch(index, c("kl","ch","hartigan","ccc","scott","marriot","trcovw","tracew","friedman",
                              "rubin","cindex","db","silhouette","duda","pseudot2","beale","ratkowsky","ball",
                              "ptbiserial","gap", "frey", "mcclain",  "gamma", "gplus", "tau", "dunn", 
                              "hubert", "sdindex", "dindex", "sdbw", "all","alllong"))
    if (is.na(indice))
      stop("invalid clustering index")
    
    if (indice == -1)
      stop("ambiguous index")
    
    if ((indice == 3)|| (indice == 5)|| (indice == 6)|| (indice == 7)|| (indice == 8)|| (indice == 9)|| (indice == 10)|| (indice == 11)|| (indice == 18)|| (indice == 27)|| (indice == 29)|| (indice == 31)|| (indice == 32))
    { 
      if((max.nc-min.nc)<2)
        stop("The difference between the minimum and the maximum number of clusters must be at least equal to 2")
    }
    
    
    
    if(is.null(data))
    {
       
      if(method==8)
      {
        stop("\n","method = kmeans, data matrix is needed")
      }
      else
      {  
        if ((indice == 1 )|| (indice == 2)|| (indice == 3 )|| (indice == 4 )|| (indice == 5)|| (indice == 6)|| (indice == 7)|| (indice == 8)|| (indice == 9)|| (indice == 10)|| (indice == 12)|| (indice == 14)|| (indice == 15)|| (indice == 16)|| (indice == 17)|| (indice == 18)
          || (indice == 19)|| (indice == 20)|| (indice == 23)|| (indice == 24)|| (indice == 25) || (indice == 27)|| (indice == 28)
          || (indice == 29) || (indice == 30)|| (indice == 31) || (indice == 32))
          stop("\n","Data matrix is needed. Only frey, mcclain, cindex, sihouette and dunn can be computed.", "\n")
       
        if(is.null(diss))
          stop("data matrix and dissimilarity matrix are both null")
        else  
          cat("\n","Only frey, mcclain, cindex, sihouette and dunn can be computed. To compute the other indices, data matrix is needed","\n") 
      }
    }
    
    else
    {
      jeu1 <- as.matrix(data)
      numberObsBefore <- dim(jeu1)[1]
      jeu <- na.omit(jeu1) # returns the object with incomplete cases removed 
      nn <- numberObsAfter <- dim(jeu)[1]
      pp <- dim(jeu)[2]    
      TT <- t(jeu)%*%jeu   
      sizeEigenTT <- length(eigen(TT)$value)
      eigenValues <- eigen(TT/(nn-1))$value
      
      # Only for indices using vv : CCC, Scott, marriot, tracecovw, tracew, friedman, rubin
      
      if (any(indice == 4) || (indice == 5) || (indice == 6) || (indice == 7) || (indice == 8) || (indice == 9) || (indice == 10) || (indice == 31) || (indice == 32))
      {
        for (i in 1:sizeEigenTT) 
        {
          if (eigenValues[i] < 0) {
            #cat(paste("There are only", numberObsAfter,"nonmissing observations out of a possible", numberObsBefore ,"observations."))
            stop("The TSS matrix is indefinite. There must be too many missing values. The index cannot be calculated.")
          } 
        }
        s1 <- sqrt(eigenValues)
        ss <- rep(1,sizeEigenTT)
        for (i in 1:sizeEigenTT) 
        {
          if (s1[i]!=0) 
            ss[i]=s1[i]
        }
        vv <- prod(ss)  
      } 
            
    }
    
     
    
  
      
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
  #                                                                                                                      #
  #                                              Distances                                                               #
  #                                                                                                                      #
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
  
if(is.null(distance))
  distanceM<-7
if(!is.null(distance))
  distanceM <- pmatch(distance, c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
  
if (is.na(distanceM)) 
{
  stop("invalid distance")
} 
    
if(is.null(diss))
{  
      
    if (distanceM == 1) 
    {
    		md <- dist(jeu, method="euclidean")	# "dist" function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
	  }
    if (distanceM == 2) 
      {
    		md <- dist(jeu, method="maximum")	
	    }
    if (distanceM == 3) 
      {
    		md <- dist(jeu, method="manhattan")	
	    }
    if (distanceM == 4) 
      {
    		md <- dist(jeu, method="canberra")	
	    }
    if (distanceM == 5) 
      {
    		md <- dist(jeu, method="binary")	
	    }
    if (distanceM == 6) 
     {
    		md <- dist(jeu, method="minkowski")	
	   }

   if (distanceM == 7) 
    {		  
     stop("dissimilarity matrix and distance are both NULL")		
    } 
}

if(!is.null(diss))
{
  if((distanceM==1)||(distanceM==2)|| (distanceM==3)|| (distanceM==4)|| (distanceM==5)|| (distanceM==6))
    stop("dissimilarity matrix and distance are both not null")
  else
    md <- diss
}
   

  
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
    #                                                                                                                      #
    #                                              Methods                                                                 #
    #                                                                                                                      #
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
  
   
    res <- array(0, c(max_nc-min_nc+1,30))
    x_axis <- min_nc:max_nc
    resCritical <- array(0, c(max_nc-min_nc+1,4))
    rownames(resCritical) <- min_nc:max_nc
    colnames(resCritical) <- c("CritValue_Duda", "CritValue_PseudoT2", "Fvalue_Beale", "CritValue_Gap")
    rownames(res) <- min_nc:max_nc
    colnames(res) <- c("KL","CH","Hartigan","CCC","Scott","Marriot", "TrCovW", "TraceW","Friedman","Rubin","Cindex","DB",
                       "Silhouette", "Duda", "Pseudot2", "Beale", "Ratkowsky", "Ball", "Ptbiserial", "Gap", "Frey", "McClain","Gamma", "Gplus", "Tau", "Dunn", 
                       "Hubert", "SDindex", "Dindex", "SDbw")   
    
    if (is.na(method))
	     stop("invalid clustering method")
    if (method == -1)
	     stop("ambiguous method")
    if (method == 1) 
    {
        hc<-hclust(md,method = "ward.D2")      
    }
    if (method == 2) 
    {
        hc<-hclust(md,method = "single")		
	  }
    if (method == 3)
     {
        hc<-hclust(md,method = "complete")		
     }
	   
    if (method == 4) 
    {
        hc<-hclust(md,method = "average")
    }
	  
    if (method == 5) 
    {
        hc<-hclust(md,method = "mcquitty")
		
   	}
    if (method == 6) 
    {
        hc<-hclust(md,method = "median")
			
	  }
    if (method == 7) 
    {
        hc<-hclust(md,method = "centroid")
		 
	  }
    if (method == 9) 
    {
      hc<-hclust(md,method = "ward.D")
  
    }

   
  

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                                                                      #
#                                              Indices                                                                 #
#                                                                                                                      #
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
    


    
##############################
#                            #
#        SD and SDbw         #
#                            #
##############################    
    

centers<-function(cl,x)
{
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    centers <- matrix(nrow = k, ncol = ncol(x))
    {
        for (i in 1:k) 
        {
            for (j in 1:ncol(x)) 
            {
                centers[i, j] <- mean(x[cl == i, j])
            }
        }
    }
    return(centers)
}    

Average.scattering <- function (cl, x)
{
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    centers.matrix <- centers(cl,x)
    
    cluster.size <- numeric(0)  
    variance.clusters <- matrix(0, ncol = ncol(x), nrow = k)
    var <- matrix(0, ncol = ncol(x), nrow = k)
    
    for (u in 1:k) 
      cluster.size[u] <- sum(cl == u)

    for (u in 1:k) 
    {  
   		for (j in 1:ncol(x)) 
      { 
			   for(i in 1:n) 
         {     				   
           if(cl[i]==u)                   
              variance.clusters[u,j]<- variance.clusters[u,j]+(x[i, j]-centers.matrix[u,j])^2 
         }
	    }            
    }

    for (u in 1:k) 
    {    
       for (j in 1:ncol(x)) 
          variance.clusters[u,j]= variance.clusters[u,j]/ cluster.size[u]   
    }
    
     
     variance.matrix <- numeric(0)
     for(j in 1:ncol(x)) 
        variance.matrix[j]=var(x[,j])*(n-1)/n

     
      Somme.variance.clusters<-0
      for (u in 1:k) 
         Somme.variance.clusters<-Somme.variance.clusters+sqrt((variance.clusters[u,]%*%(variance.clusters[u,])))
         

      # Standard deviation
      stdev<-(1/k)*sqrt(Somme.variance.clusters)
      
      #Average scattering for clusters  
      scat<- (1/k)* (Somme.variance.clusters /sqrt(variance.matrix %*% variance.matrix))
      
      scat <- list(stdev=stdev, centers=centers.matrix, variance.intraclusters= variance.clusters, scatt=scat)
      return(scat)
}

density.clusters<-function(cl, x)
{
   x <- as.matrix(x)
   k <- max(cl)
   n <- length(cl)
         
   distance <- matrix(0, ncol = 1, nrow = n)
   density <-  matrix(0, ncol = 1, nrow = k)
   centers.matrix<-centers(cl,x)
   stdev<-Average.scattering(cl,x)$stdev 
   for(i in 1:n) 
   {        
       u=1
       while(cl[i] != u )
          u<-u+1
       for (j in 1:ncol(x))   
       {               
           distance[i]<- distance[i]+(x[i,j]-centers.matrix[u,j])^2 
       }     
       distance[i]<-sqrt(distance[i])            
       if (distance[i] <= stdev)
          density[u]= density[u]+1                      
   }  
    dens<-list(distance=distance, density=density)    
	  return(dens)          
 
}


density.bw<-function(cl, x)
{
   x <- as.matrix(x)
   k <- max(cl)
   n <- length(cl)   
   centers.matrix<-centers(cl,x)
   stdev<-Average.scattering(cl,x)$stdev 
   density.bw<- matrix(0, ncol = k, nrow = k)
   u<-1
   
   for(u in 1:k)
   {
     for(v in 1:k)
     {
       if(v!=u)
       {  
          distance<- matrix(0, ncol = 1, nrow = n)
          moy<-(centers.matrix[u,]+centers.matrix[v,])/2
          for(i in 1:n)
          {
            if((cl[i]==u)||(cl[i]==v))
            {
              for (j in 1:ncol(x))   
              {               
                 distance[i]<- distance[i]+(x[i,j]-moy[j])^2 
              }   
              distance[i]<- sqrt(distance[i])
              if(distance[i]<= stdev)
              {
                density.bw[u,v]<-density.bw[u,v]+1                  
              }  
            }           
          }
        }       
       }
      }
     density.clust<-density.clusters(cl,x)$density 
     S<-0
     for(u in 1:k)
       for(v in 1:k)
       {  
         if(max(density.clust[u], density.clust[v])!=0)
            S=S+ (density.bw[u,v]/max(density.clust[u], density.clust[v]))
       }   
     density.bw<-S/(k*(k-1))
     return(density.bw) 
  
 }      
           
            

Dis <- function (cl, x)
{   # Dis : Total separation between clusters

    x <- as.matrix(x)
    k <- max(cl)
    centers.matrix <- centers(cl,x)
    Distance.centers<-dist(centers.matrix)
    Dmin<-min(Distance.centers)
    Dmax<-max(Distance.centers)
    Distance.centers<-as.matrix(Distance.centers)
    s2<-0
    for (u in 1:k)
    {
       s1=0
       for(j in 1:ncol(Distance.centers))
       {s1<-s1 + Distance.centers[u,j]       
       }
       s2<-s2+1/s1      
    }  
    Dis<-(Dmax/Dmin)*s2  
    return(Dis)
}  
 

    
##################################
#                                #  
#         Hubert index           #
#                                #
##################################
    
    
    
Index.Hubert<-function(x, cl)
{
      
      k <- max(cl)
      n<-dim(x)[1]
      y <- matrix(0, ncol = dim(x)[2], nrow = n)
      P<- as.matrix(md)
      meanP<-mean(P)
      variance.matrix <- numeric(0)
      md <- dist(x, method="euclidean")
      for(j in 1:n) 
        variance.matrix[j]=var(P[,j])*(n-1)/n
      varP<-sqrt(variance.matrix %*% variance.matrix)
      
      centers.clusters<-centers(cl,x)
      for(i in 1:n)
      {
        for(u in 1:k)
        {
          if(cl[i]==u)
            y[i,]<-centers.clusters[u,]
        }   
      }  
      
      Q<- as.matrix(dist(y, method="euclidean"))
      meanQ<-mean(Q)
      for(j in 1:n) 
        variance.matrix[j]=var(Q[,j])*(n-1)/n
      varQ<-sqrt(variance.matrix %*% variance.matrix)
      
      M<-n*(n-1)/2
      S<-0
      n1<-n-1
      
      for(i in 1:n1)
      { 
        j<-i+1
        while(j<=n)
        {
          S<-S+(P[i,j]-meanP)*(Q[i,j]-meanQ)
          j<-j+1
        }
        
      } 
      gamma<-S/(M*varP*varQ)
      
      return(gamma)
}   



##################################
#                                #  
#   Gamma, Gplus and Tau         #
#                                #
##################################    


Index.sPlussMoins <- function (cl1,md)
{
    cn1 <- max(cl1)
    n1 <- length(cl1)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
    di <- list()
    for (u in 1:cn1) {
        cluster.size[u] <- sum(cl1 == u)
        du <- as.dist(dmat[cl1 == u, cl1 == u])
        within.dist1 <- c(within.dist1, du)
        average.distance[u] <- mean(du)
        median.distance[u] <- median(du)
        bv <- numeric(0)
        for (v in 1:cn1) {
            if (v != u) {
                suv <- dmat[cl1 == u, cl1 == v]
                bv <- c(bv, suv)
                if (u < v) {
                  separation.matrix[u, v] <- separation.matrix[v,u] <- min(suv)
                  between.dist1 <- c(between.dist1, suv)
                }
            }
        }
    }

    nwithin1 <- length(within.dist1)
    nbetween1 <- length(between.dist1)
    meanwithin1 <- mean(within.dist1)
    meanbetween1 <- mean(between.dist1)
    
    s.plus <- s.moins <- 0 
    #s.moins<-sum(rank(c(within.dist1,between.dist1),ties="first")[1:nwithin1]-rank(within.dist1,ties="first"))
    #s.plus  <-sum(rank(c(-within.dist1,-between.dist1),ties="first")[1:nwithin1]-rank(-within.dist1,ties="first"))
    for (k in 1: nwithin1)
    {
      s.plus <- s.plus+(colSums(outer(between.dist1,within.dist1[k], ">")))
      s.moins <- s.moins+(colSums(outer(between.dist1,within.dist1[k], "<")))
    }    
    
    Index.Gamma <- (s.plus-s.moins)/(s.plus+s.moins)
    Index.Gplus <- (2*s.moins)/(n1*(n1-1))
    t.tau  <- (nwithin1*nbetween1)-(s.plus+s.moins)
    Index.Tau <- (s.plus-s.moins)/(((n1*(n1-1)/2-t.tau)*(n1*(n1-1)/2))^(1/2))

    results <- list(gamma=Index.Gamma, gplus=Index.Gplus, tau=Index.Tau)
    return(results)
}



    
##################################
#                                #  
#      Frey and McClain          #
#                                #
################################## 
    
    
    

Index.15and28  <- function (cl1,cl2,md)
{
    cn1 <- max(cl1)
    n1 <- length(cl1)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
    di <- list()
    for (u in 1:cn1) 
      {
        cluster.size[u] <- sum(cl1 == u)
        du <- as.dist(dmat[cl1 == u, cl1 == u])
        within.dist1 <- c(within.dist1, du)
        #average.distance[u] <- mean(du)
        #median.distance[u] <- median(du)
        #bv <- numeric(0)
        for (v in 1:cn1) {
            if (v != u) {
                suv <- dmat[cl1 == u, cl1 == v]
                #bv <- c(bv, suv)
                if (u < v) {
                  separation.matrix[u, v] <- separation.matrix[v,u] <- min(suv)
                  between.dist1 <- c(between.dist1, suv)
                }
            }
        }
    }
    cn2 <- max(cl2)
    n2 <- length(cl2)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist2 <- between.dist2 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn2, nrow = cn2)
    di <- list()
    for (w in 1:cn2) {
        cluster.size[w] <- sum(cl2 == w)
        dw <- as.dist(dmat[cl2 == w, cl2 == w])
        within.dist2 <- c(within.dist2, dw)
        #average.distance[w] <- mean(dw)
        #median.distance[w] <- median(dw)
        bx <- numeric(0)
        for (x in 1:cn2) {
            if (x != w) {
                swx <- dmat[cl2 == w, cl2 == x]
                bx <- c(bx, swx)
                if (w < x) {
                  separation.matrix[w, x] <- separation.matrix[x,w] <- min(swx)
                  between.dist2 <- c(between.dist2, swx)
                }
            }
        }
    }
    nwithin1 <- length(within.dist1)
    nbetween1 <- length(between.dist1)
    meanwithin1 <- mean(within.dist1)
    meanbetween1 <- mean(between.dist1)
    meanwithin2 <- mean(within.dist2)
    meanbetween2 <- mean(between.dist2)
    Index.15 <- (meanbetween2-meanbetween1)/(meanwithin2-meanwithin1)
    Index.28 <- (meanwithin1/nwithin1)/(meanbetween1/nbetween1)

    results <- list(frey=Index.15,mcclain=Index.28)
    return(results)
}

    
##################################
#                                #  
#      Point-biserial            #
#                                #
##################################  
    
  
    
Indice.ptbiserial <- function (x,md,cl1)
{
	nn <- dim(x)[1]
	pp <- dim(x)[2]

	md2 <- as.matrix(md)
	m01 <- array(NA, c(nn,nn))
	nbr <- (nn*(nn-1))/2
	pb <- array(0,c(nbr,2))
	
	m3 <- 1
	for (m1 in 2:nn)
	{
	     m12 <- m1-1
	   for (m2 in 1:m12)
	   {
		if (cl1[m1]==cl1[m2]) m01[m1,m2]<-0
		if (cl1[m1]!=cl1[m2]) m01[m1,m2]<-1
		pb[m3,1] <- m01[m1,m2]
		pb[m3,2] <- md2[m1,m2]
		m3 <- m3+1
	   }
	}

	y <- pb[,1]
	x <- pb[,2] 

	biserial.cor <- function (x, y, use = c("all.obs", "complete.obs"), level = 1) 
	{
	    if (!is.numeric(x)) 
	        stop("'x' must be a numeric variable.\n")
	    y <- as.factor(y)
	    if (length(levs <- levels(y)) > 2) 
	        stop("'y' must be a dichotomous variable.\n")
	    if (length(x) != length(y)) 
	        stop("'x' and 'y' do not have the same length")
	    use <- match.arg(use)
	    if (use == "complete.obs") {
	        cc.ind <- complete.cases(x, y)
	        x <- x[cc.ind]
	        y <- y[cc.ind]
	    }
	    ind <- y == levs[level]
	    diff.mu <- mean(x[ind]) - mean(x[!ind])
	    prob <- mean(ind)
	    diff.mu * sqrt(prob * (1 - prob))/sd(x)
	}

    ptbiserial <- biserial.cor(x=pb[,2], y=pb[,1], level = 2)
    return(ptbiserial)
}

    
##########################################
#                                        #
#       Duda, pseudot2 and beale         #
#                                        #
##########################################
    

Indices.WKWL <- function (x,cl1=cl1,cl2=cl2)
{
   dim2 <- dim(x)[2]
   wss <- function(x) 
    {
	      x <- as.matrix(x)
        n <- length(x)
        centers <- matrix(nrow = 1, ncol = ncol(x))

        if (ncol(x) == 1) 
          {	centers[1, ] <- mean(x) 	}
        if (is.null(dim(x))) 
          {
		      bb <- matrix(x,byrow=FALSE,nrow=1,ncol=ncol(x))
        	centers[1, ] <- apply(bb, 2, mean)
    		  }
    		else 
          {
                centers[1, ] <- apply(x, 2, mean)
		      }

        x.2 <- sweep(x,2,centers[1,],"-")
        withins <- sum(x.2^2)
        wss <- sum(withins)
        return(wss)
    }

     ncg1 <- 1
     ncg1max <- max(cl1)
     while((sum(cl1==ncg1)==sum(cl2==ncg1)) && ncg1 <=ncg1max) 
     { ncg1 <- ncg1+1 }
     g1 <- ncg1

     ncg2 <- max(cl2)
     nc2g2 <- ncg2-1
     while((sum(cl1==nc2g2)==sum(cl2==ncg2)) && nc2g2 >=1) 
     { 
	     ncg2 <- ncg2-1 
	     nc2g2 <- nc2g2-1
     }
     g2 <- ncg2

     NK <- sum(cl2==g1)
     WK.x <- x[cl2==g1,]
     WK <- wss(x=WK.x)

     NL <- sum(cl2==g2)
     WL.x <- x[cl2==g2,]
     WL <- wss(x=WL.x)

     NM <- sum(cl1==g1)
     WM.x <- x[cl1==g1,]
     WM <- wss(x=WM.x)

     duda <- (WK+WL)/WM

     BKL <- WM-WK-WL
     pseudot2 <- BKL/((WK+WL)/(NK+NL-2))

     beale <- (BKL/(WK+WL))/(((NM-1)/(NM-2))*(2^(2/dim2)-1))

    results <- list(duda=duda,pseudot2=pseudot2,NM=NM,NK=NK,NL=NL,beale=beale)
    return(results)
}


########################################################################
#                                                                      #
#       ccc, scott, marriot, trcovw, tracew, friedman and rubin        #
#                                                                      #
########################################################################    
    
   

Indices.WBT <- function(x,cl,P,s,vv) 
{
  n <- dim(x)[1]
  pp <- dim(x)[2]
  qq <- max(cl)
  z <- matrix(0,ncol=qq,nrow=n)
  clX <- as.matrix(cl)

  for (i in 1:n)
    for (j in 1:qq)
    {
	    z[i,j]==0
	    if (clX[i,1]==j) 
	    {z[i,j]=1}
    }

  xbar <- solve(t(z)%*%z)%*%t(z)%*%x
  B <- t(xbar)%*%t(z)%*%z%*%xbar
  W <- P-B
  marriot <- (qq^2)*det(W)
  trcovw <- sum(diag(cov(W)))
  tracew <- sum(diag(W))
  if(det(W)!=0)
     scott <- n*log(det(P)/det(W))
  else {cat("Error: division by zero!")}
  friedman <- sum(diag(solve(W)*B))
  rubin <- sum(diag(P))/sum(diag(W))
  

  R2 <- 1-sum(diag(W))/sum(diag(P))
  v1 <- 1
  u <- rep(0,pp)
  c <- (vv/(qq))^(1/pp)
  u <- s/c
  k1 <- sum((u>=1)==TRUE)
  p1 <- min(k1,qq-1)
  if (all(p1>0,p1<pp))
  {
    for (i in 1:p1)
    v1 <- v1*s[i]
    c <- (v1/(qq))^(1/p1)
    u <- s/c
    b1 <- sum(1/(n+u[1:p1]))
    b2 <- sum(u[p1+1:pp]^2/(n+u[p1+1:pp]),na.rm=TRUE)
    E_R2 <- 1-((b1+b2)/sum(u^2))*((n-qq)^2/n)*(1+4/n)
    ccc <- log((1-E_R2)/(1-R2))*(sqrt(n*p1/2)/((0.001+E_R2)^1.2))
  }else 
  {
    b1 <- sum(1/(n+u))
    E_R2 <- 1-(b1/sum(u^2))*((n-qq)^2/n)*(1+4/n)
    ccc <- log((1-E_R2)/(1-R2))*(sqrt(n*pp/2)/((0.001+E_R2)^1.2))
  }
 results <- list(ccc=ccc,scott=scott,marriot=marriot,trcovw=trcovw,tracew=tracew,friedman=friedman,rubin=rubin)
 return(results)
}

   

    ########################################################################
    #                                                                      #
    #                   Kl, Ch, Hartigan, Ratkowsky and Ball               #
    #                                                                      #
    ########################################################################     
    
    


Indices.Traces <- function(data, d, clall, index = "all") 
{
	x <- data
	cl0 <- clall[,1]
	cl1 <- clall[,2]
	cl2 <- clall[,3]
	clall <- clall
	nb.cl0 <- table(cl0)
	nb.cl1 <- table(cl1)
	nb.cl2 <- table(cl2)
	nb1.cl0 <- sum(nb.cl0==1)
	nb1.cl1 <- sum(nb.cl1==1)
	nb1.cl2 <- sum(nb.cl2==1)

  gss <- function(x, cl, d) 
    {
        n <- length(cl)
        k <- max(cl)
        centers <- matrix(nrow = k, ncol = ncol(x))
        for (i in 1:k) 
          {

            	  if (ncol(x) == 1)
                 {
                  	centers[i, ] <- mean(x[cl == i, ])
            	   }
                if (is.null(dim(x[cl == i, ])))
                  {
		                bb <- matrix(x[cl == i, ],byrow=FALSE,nrow=1,ncol=ncol(x))
        	          centers[i, ] <- apply(bb, 2, mean)
    	          	}
    		        else 
                  {
                    centers[i, ] <- apply(x[cl == i, ], 2, mean)
		              }

        }
    	allmean <- apply(x, 2, mean)
    	dmean <- sweep(x, 2, allmean, "-")
    	allmeandist <- sum(dmean^2)
      withins <- rep(0, k)
      x.2 <- (x - centers[cl, ])^2
      for (i in 1:k) {
            withins[i] <- sum(x.2[cl == i, ])
        }
      wgss <- sum(withins)
    	bgss <- allmeandist - wgss

    results <- list(wgss=wgss,bgss=bgss, centers=centers)
    return(results)
    }

    vargss <- function(x, clsize, varwithins) 
    {
        nvar <- dim(x)[2]
        n <- sum(clsize)
        k <- length(clsize)
        varallmean <- rep(0, nvar)
        varallmeandist <- rep(0, nvar)
        varwgss <- rep(0, nvar)
        for (l in 1:nvar) varallmean[l] <- mean(x[, l])
        vardmean <- sweep(x, 2, varallmean, "-")
        for (l in 1:nvar) {
            varallmeandist[l] <- sum((vardmean[, l])^2)
            varwgss[l] <- sum(varwithins[, l])
        }
        varbgss <- varallmeandist - varwgss
        vartss <- varbgss + varwgss
        zvargss <- list(vartss = vartss, varbgss = varbgss)
        return(zvargss)
    }
    varwithinss <- function(x, centers, cluster) {
        nrow <- dim(centers)[1]
        nvar <- dim(x)[2]
        varwithins <- matrix(0, nrow, nvar)
        x <- (x - centers[cluster, ])^2
        for (l in 1:nvar) {
            for (k in 1:nrow) {
                varwithins[k, l] <- sum(x[cluster == k, l])
            }
        }
        return(varwithins)
    }



    indice.kl <- function (x, clall, d = NULL, centrotypes = "centroids"){
	if (nb1.cl1 > 0){
		KL <- NA
		}
    	    if (sum(c("centroids", "medoids") == centrotypes) == 0) 
       	     stop("Wrong centrotypes argument")
    	    if ("medoids" == centrotypes && is.null(d)) 
             stop("For argument centrotypes = 'medoids' d cannot be null")
    	   if (!is.null(d)) {
            if (!is.matrix(d)) {
               d <- as.matrix(d)
            }
            row.names(d) <- row.names(x)
    	    }
    	    #if (is.null(dim(x))) {
            #	    dim(x) <- c(length(x), 1)
    	    #}
    	    m <- ncol(x)
    	    g <- k <- max(clall[, 2])
    	    KL <- abs((g - 1)^(2/m) * gss(x, clall[, 1], d)$wgss - 
        	g^(2/m) * gss(x, clall[, 2], d)$wgss)/abs((g)^(2/m) * 
       	 	gss(x, clall[, 2], d)$wgss - (g + 1)^(2/m) * 
       		gss(x, clall[, 3], d)$wgss)
	    return(KL)
    }



    indice.ch <- function (x, cl, d = NULL, centrotypes = "centroids"){
	if (nb1.cl1 > 0){
		CH <- NA
		}
    	    if (sum(c("centroids", "medoids") == centrotypes) == 0) 
       	     stop("Wrong centrotypes argument")
    	    if ("medoids" == centrotypes && is.null(d)) 
             stop("For argument centrotypes = 'medoids' d cannot be null")
    	   if (!is.null(d)) {
            if (!is.matrix(d)) {
               d <- as.matrix(d)
            }
            row.names(d) <- row.names(x)
    	    }
    	    #if (is.null(dim(x))) {
            #	    dim(x) <- c(length(x), 1)
    	    #}
    	    n <- length(cl)
    	    k <- max(cl)
	    CH <- (gss(x, cl, d)$bgss/(k-1))/
		(gss(x, cl, d)$wgss/(n-k))
	    return(CH)
    }
    

# hartigan

    indice.hart <- function(x, clall, d = NULL, centrotypes = "centroids"){
    	if (sum(c("centroids", "medoids") == centrotypes) == 0) 
    	    stop("Wrong centrotypes argument")
    	if ("medoids" == centrotypes && is.null(d)) 
    	    stop("For argument centrotypes = 'medoids' d cannot be null")
    	if (!is.null(d)) {
    	    if (!is.matrix(d)) {
    	        d <- as.matrix(d)
    	    }
    	    row.names(d) <- row.names(x)
    	}
    	#if (is.null(dim(x))) {
    	#    dim(x) <- c(length(x), 1)
    	#}
    	n <- nrow(x)
    	g <- max(clall[, 1])
    	HART <- (gss(x, clall[, 2], d)$wgss/gss(x, clall[, 3],d)$wgss - 1) * (n - g - 1)
	return(HART)
    }
    
   
    

    indice.ball <- function(x, cl, d = NULL, centrotypes = "centroids"){
  	wgssB <- gss(x, cl, d)$wgss
  	qq <- max(cl)
	  ball <- wgssB/qq
	  return(ball)
    }




    indice.ratkowsky <- function(x, cl, d, centrotypes = "centroids"){
	  qq <- max(cl)
	  clsize <- table(cl)
	  centers <- gss(x, cl, d)$centers
    varwithins <- varwithinss(x, centers, cl)
    zvargss <- vargss(x, clsize, varwithins)
    ratio <- mean(sqrt(zvargss$varbgss/zvargss$vartss))
    ratkowsky <- ratio/sqrt(qq)
    return(ratkowsky)
    }

    indice <- pmatch(index, c("kl", "ch", "hart", "ratkowsky", "ball", "all"))
    if (is.na(indice)) 
        stop("invalid clustering index")
    if (indice == -1) 
        stop("ambiguous index")
    vecallindex <- numeric(5)
    if (any(indice == 1) || (indice == 6)) 
        vecallindex[1] <- indice.kl(x,clall,d)
    if (any(indice == 2) || (indice == 6)) 
        vecallindex[2] <- indice.ch(x,cl=clall[,2],d)
    if (any(indice == 3) || (indice == 6)) 
        vecallindex[3] <- indice.hart(x,clall,d)
    if (any(indice == 4) || (indice == 6)) 
        vecallindex[4] <- indice.ratkowsky(x,cl=cl1, d)
    if (any(indice == 5) || (indice == 6)) 
        vecallindex[5] <- indice.ball(x,cl=cl1,d)
    names(vecallindex) <- c("kl", "ch", "hart", "ratkowsky", "ball")
    if (indice < 6) 
        vecallindex <- vecallindex[indice]
    return(vecallindex)

}


    
    ########################################################################
    #                                                                      #
    #                              C-index                                 #
    #                                                                      #
    ######################################################################## 
     
    
    
Indice.cindex <- function (d, cl) 
{
    d <- data.matrix(d)
    DU <- 0
    r <- 0
    v_max <- array(1, max(cl))
    v_min <- array(1, max(cl))
    for (i in 1:max(cl)) {
        n <- sum(cl == i)
        if (n > 1) {
            t <- d[cl == i, cl == i]
            DU = DU + sum(t)/2
            v_max[i] = max(t)
            if (sum(t == 0) == n) 
                v_min[i] <- min(t[t != 0])
            else v_min[i] <- 0
            r <- r + n * (n - 1)/2
        }
    }
    Dmin = min(v_min)
    Dmax = max(v_max)
    if (Dmin == Dmax) 
        result <- NA
    else result <- (DU - r * Dmin)/(Dmax * r - Dmin * r)
    result
}


    
    ########################################################################
    #                                                                      #
    #                                 DB                                   #
    #                                                                      #
    ########################################################################     
    
    
Indice.DB <- function (x, cl, d = NULL, centrotypes = "centroids", p = 2, q = 2) 
{
    if (sum(c("centroids") == centrotypes) == 0) 
        stop("Wrong centrotypes argument")
    if (!is.null(d)) {
        if (!is.matrix(d)) {
            d <- as.matrix(d)
        }
        row.names(d) <- row.names(x)
    }
    if (is.null(dim(x))) {
        dim(x) <- c(length(x), 1)
    }
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    dAm <- d
    centers <- matrix(nrow = k, ncol = ncol(x))
    if (centrotypes == "centroids") {
        for (i in 1:k) {
            for (j in 1:ncol(x)) {
                centers[i, j] <- mean(x[cl == i, j])
            }
        }
    }
    else {
        stop("wrong centrotypes argument")
    }
    S <- rep(0, k)
    for (i in 1:k) {
        ind <- (cl == i)
        if (sum(ind) > 1) {
            centerI <- centers[i, ]
            centerI <- rep(centerI, sum(ind))
            centerI <- matrix(centerI, nrow = sum(ind), ncol = ncol(x), 
                byrow = TRUE)
            S[i] <- mean(sqrt(apply((x[ind, ] - centerI)^2, 1, 
                sum))^q)^(1/q)
        }
        else S[i] <- 0
    }
    M <- as.matrix(dist(centers, p = p))
    R <- array(Inf, c(k, k))
    r = rep(0, k)
    for (i in 1:k) {
        for (j in 1:k) {
            R[i, j] = (S[i] + S[j])/M[i, j]
        }
        r[i] = max(R[i, ][is.finite(R[i, ])])
    }
    DB = mean(r[is.finite(r)])
    resul <- list(DB = DB, r = r, R = R, d = M, S = S, centers = centers)
    resul
}


    
    
    ########################################################################
    #                                                                      #
    #                             Silhouette                               #
    #                                                                      #
    ########################################################################     
    
    

Indice.S <- function (d, cl) 
{
    d <- as.matrix(d)
    Si <- 0
    for (k in 1:max(cl)) {
        if ((sum(cl == k)) <= 1) 
            Sil <- 1
        else {
            Sil <- 0
            for (i in 1:length(cl)) {
                if (cl[i] == k) {
                  ai <- sum(d[i, cl == k])/(sum(cl == k) - 1)
                  dips <- NULL
                  for (j in 1:max(cl)) if (cl[i] != j) 
                    if (sum(cl == j) != 1) 
                      dips <- cbind(dips, c((sum(d[i, cl == j]))/(sum(cl == 
                        j))))
                    else dips <- cbind(dips, c((sum(d[i, cl == 
                      j]))))
                  bi <- min(dips)
                  Sil <- Sil + (bi - ai)/max(c(ai, bi))
                }
            }
        }
        Si <- Si + Sil
    }
    Si/length(cl)
}


    
    ########################################################################
    #                                                                      #
    #                                  Gap                                 #
    #                                                                      #
    ######################################################################## 
    
Indice.Gap <- function (x, clall, reference.distribution = "unif", B = 10, 
    method = "ward.D2", d = NULL, centrotypes = "centroids") 
{
    GAP <- function(X, cl, referenceDistribution, B, method, d, centrotypes) 
      {
        set.seed(1)
        simgap <- function(Xvec) 
          {
            ma <- max(Xvec)
            mi <- min(Xvec)
            set.seed(1)
            Xout <- runif(length(Xvec), min = mi, max = ma)
            return(Xout)
          }
          pcsim <- function(X, d, centrotypes) 
          {
            if (centrotypes == "centroids") 
            {
                Xmm <- apply(X, 2, mean)
            }
            
            for (k in (1:dim(X)[2])) 
            {
                X[, k] <- X[, k] - Xmm[k]
            }
            ss <- svd(X)
            Xs <- X %*% ss$v
            Xnew <- apply(Xs, 2, simgap)
            Xt <- Xnew %*% t(ss$v)
            for (k in (1:dim(X)[2])) {
                Xt[, k] <- Xt[, k] + Xmm[k]
            }
            return(Xt)
        }
        if (is.null(dim(x))) 
        {
            dim(x) <- c(length(x), 1)
        }
        ClassNr <- max(cl)
        Wk0 <- 0
        WkB <- matrix(0, 1, B)
        for (bb in (1:B)) {
            if (reference.distribution == "unif") 
                Xnew <- apply(X, 2, simgap)
            else if (reference.distribution == "pc") 
                Xnew <- pcsim(X, d, centrotypes)
            else stop("Wrong reference distribution type")
            if (bb == 1) {
                pp <- cl
                if (ClassNr == length(cl)) 
                  pp2 <- 1:ClassNr
                else if (method == "k-means") 
                { set.seed(1)
                  pp2 <- kmeans(Xnew, ClassNr, 100)$cluster
                }
                else if (method == "single" || method == "complete" || 
                  method == "average" || method == "ward.D2" || 
                  method == "mcquitty" || method == "median" || 
                  method == "centroid"|| method=="ward.D") 
                  pp2 <- cutree(hclust(dist(Xnew), method = method), 
                    ClassNr)
                else stop("Wrong clustering method")
                if (ClassNr > 1) {
                  for (zz in (1:ClassNr)) {
                    Xuse <- X[pp == zz, ]
                    Wk0 <- Wk0 + sum(diag(var(Xuse))) * (length(pp[pp == 
                      zz]) - 1)/(dim(X)[1] - ClassNr)
                    Xuse2 <- Xnew[pp2 == zz, ]
                    WkB[1, bb] <- WkB[1, bb] + sum(diag(var(Xuse2))) * 
                      (length(pp2[pp2 == zz]) - 1)/(dim(X)[1] - 
                      ClassNr)
                  }
                }
                if (ClassNr == 1) 
                {
                  Wk0 <- sum(diag(var(X)))
                  WkB[1, bb] <- sum(diag(var(Xnew)))
                }
             }
             if (bb > 1) {
                if (ClassNr == length(cl)) 
                  pp2 <- 1:ClassNr
                else if (method == "k-means")
                {
                  set.seed(1)
                  pp2 <- kmeans(Xnew, ClassNr, 100)$cluster
                }
                else if (method == "single" || method == "complete" || 
                  method == "average" || method == "ward.D2" || 
                  method == "mcquitty" || method == "median" || 
                  method == "centroid"||method == "ward.D") 
                  pp2 <- cutree(hclust(dist(Xnew), method = method), 
                    ClassNr)
                else stop("Wrong clustering method")
                if (ClassNr > 1) {
                  for (zz in (1:ClassNr)) {
                    Xuse2 <- Xnew[pp2 == zz, ]
                    WkB[1, bb] <- WkB[1, bb] + sum(diag(var(Xuse2))) * 
                      length(pp2[pp2 == zz])/(dim(X)[1] - ClassNr)
                  }
                }
                if (ClassNr == 1) {
                  WkB[1, bb] <- sum(diag(var(Xnew)))
                }
            }
        }
        Sgap <- mean(log(WkB[1, ])) - log(Wk0)
        Sdgap <- sqrt(1 + 1/B) * sqrt(var(log(WkB[1, ]))) * sqrt((B - 
            1)/B)
        resul <- list(Sgap = Sgap, Sdgap = Sdgap)
        resul
    }
    if (sum(c("centroids", "medoids") == centrotypes) == 0) 
        stop("Wrong centrotypes argument")
    if ("medoids" == centrotypes && is.null(d)) 
        stop("For argument centrotypes = 'medoids' d can not be null")
    if (!is.null(d)) {
        if (!is.matrix(d)) {
            d <- as.matrix(d)
        }
        row.names(d) <- row.names(x)
    }
    X <- as.matrix(x)
    gap1 <- GAP(X, clall[, 1], reference.distribution, B, method, 
        d, centrotypes)
    gap <- gap1$Sgap
    gap2 <- GAP(X, clall[, 2], reference.distribution, B, method, 
        d, centrotypes)
    diffu <- gap - (gap2$Sgap - gap2$Sdgap)
    resul <- list(gap = gap, diffu = diffu)
    resul

}
    
    
  
    
    ########################################################################
    #                                                                      #
    #                              SD, sdbw, dunn                          #
    #                                                                      #
    ########################################################################   
    
    
    
    
    Index.sdindex<-function(x, clmax, cl)
    {  
      x <- as.matrix(x)
      Alpha<-Dis(clmax,x)
      Scatt<-Average.scattering(cl,x)$scatt
      Dis0<-Dis(cl,x)
      SD.indice<-Alpha*Scatt + Dis0
      return(SD.indice)
    }
    
    Index.SDbw<-function(x, cl)
    {
      x <- as.matrix(x)
      Scatt<-Average.scattering(cl,x)$scatt
      Dens.bw<-density.bw(cl,x)
      SDbw<-Scatt+Dens.bw
      return(SDbw)
    }    
    
  
    
    ########################################################################
    #                                                                      #
    #                              D index                                 #
    #                                                                      #
    ######################################################################## 
   
    
    
    
    Index.Dindex<- function(cl, x)
    {
      x <- as.matrix(x)
      distance<-density.clusters(cl, x)$distance
      n<-length(distance)
      S<-0
      for(i in 1:n)
        S<-S+distance[i]
      inertieIntra<-S/n
      return(inertieIntra)
    }    
    
    
    #####################################################################
    #                                                                   #
    #                            Dunn index                             #
    #                                                                   #
    #####################################################################
    

    
    Index.dunn <- function(md, clusters, Data=NULL, method="euclidean")
    {
      
      distance <- as.matrix(md)
      nc <- max(clusters)
      interClust <- matrix(NA, nc, nc)
      intraClust <- rep(NA, nc)
      
      for (i in 1:nc) 
      {
        c1 <- which(clusters==i)
        for (j in i:nc) {
          if (j==i) intraClust[i] <- max(distance[c1,c1])
          if (j>i) {
            c2 <- which(clusters==j)
            interClust[i,j] <- min(distance[c1,c2])
          }
        }
      }
      dunn <- min(interClust,na.rm=TRUE)/max(intraClust)
      return(dunn)
    }
    


################


 for (nc in min_nc:max_nc)
 {  
      
	   if (any(method == 1) || (method == 2) || (method == 3) || (method == 4) || 
		  (method == 5) || (method == 6) || (method == 7)||(method == 9)) 
      {
	      cl1 <- cutree(hc, k=nc)
	      cl2 <- cutree(hc, k=nc+1)
        clall <- cbind(cl1, cl2)
        clmax <- cutree(hc, k=max_nc) 
      
           
	      if (nc >= 2)
		    {
		      cl0 <- cutree(hc, k=nc-1)
          clall1 <- cbind(cl0, cl1, cl2)
	      }
	      if (nc == 1)
		    {
		      cl0 <-rep(NA,nn)
		      clall1 <- cbind(cl0, cl1, cl2)
		    }
	    }
	   
  	if (method == 8) 
    {
      set.seed(1)
		  cl2 <- kmeans(jeu,nc+1)$cluster
      set.seed(1)
		  clmax <- kmeans(jeu,max_nc)$cluster
      if (nc > 2)
		  {
        set.seed(1)
		    cl1 <- kmeans(jeu,nc)$cluster
		    clall <- cbind(cl1, cl2)
        set.seed(1)
		    cl0 <- kmeans(jeu,nc-1)$cluster
		    clall1 <- cbind(cl0, cl1, cl2)
		  }
	    if (nc == 2)
		  {
        set.seed(1)
	      cl1 <- kmeans(jeu,nc)$cluster
		    clall <- cbind(cl1, cl2)
		    cl0 <- rep(1,nn)
		    clall1 <- cbind(cl0, cl1, cl2)
	 	   }
	    if (nc == 1)
		  {
	      stop("Number of clusters must be higher than 2")
  	  }

	 }

	j <- table(cl1)  # table uses the cross-classifying factors to build a contingency table of the counts at each combination of factor levels.
	s <- sum(j==1)    
	j2 <- table(cl2)
	s2 <- sum(j2==1)
 
 
  ########### Indices.Traces-hartigan - 3e colonne de res ############ 
 	if (any(indice == 3) || (indice == 31) || (indice == 32))
	{    
	  res[nc-min_nc+1,3] <- Indices.Traces(jeu, md, clall1, index = "hart")	
	} 
  
   ########### Cubic Clustering Criterion-CCC  - 4e colonne de res ############
  if (any(indice == 4) || (indice == 31) || (indice == 32))
	{    	  
	  res[nc-min_nc+1,4] <- Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$ccc
	}

  ########### Scott and Symons - 5e colonne de res ############
	if (any(indice == 5) || (indice == 31) || (indice == 32))
	{     	  
	  res[nc-min_nc+1,5] <- Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$scott
	}

	########### Marriot - 6e colonne de res ############
	if (any(indice == 6) || (indice == 31) || (indice == 32))
	{   	  
	  res[nc-min_nc+1,6] <- Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$marriot
	}	
	
	########### Trace Cov W - 7e colonne de res ############
	if (any(indice == 7) || (indice == 31) || (indice == 32))
	{   	 
	  res[nc-min_nc+1,7] <- Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$trcovw	  
	}

  ########### Trace W - 8e colonne de res ############
	if (any(indice == 8) || (indice == 31) || (indice == 32))
	{  	  
	  res[nc-min_nc+1,8] <- Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$tracew
	}
	
	########### Friedman - 9e colonne de res ############
 	if (any(indice == 9) || (indice == 31) || (indice == 32))
	{     	  
	  res[nc-min_nc+1,9] <- Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$friedman
	}
          
  ########### Rubin - 10e colonne de res ############
   if (any(indice == 10) || (indice == 31) || (indice == 32))
	{     	  
    res[nc-min_nc+1,10] <- Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$rubin
	}
                    
 
  ########### Indices.WKWL-duda - 14e colonne de res ############
	if (any(indice == 14) || (indice == 31) || (indice == 32))
	{  
	  res[nc-min_nc+1,14] <- Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$duda	
	}
	
	
  ########### Indices.WKWL-pseudot2 - 15e colonne de res ############
	if (any(indice == 15) || (indice == 31) || (indice == 32))
	{   	  
      res[nc-min_nc+1,15] <- Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$pseudot2	
	}
  
  ########### Indices.WKWL-beale - 16e colonne de res ############
	if (any(indice == 16) || (indice == 31) || (indice == 32))
	{             	  
	  res[nc-min_nc+1,16] <- beale <- Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$beale
	}

  ########### Indices.WKWL- duda or pseudot2 or beale ############   
  if (any(indice == 14) || (indice == 15) || (indice == 16) || (indice == 31) || (indice == 32))
	{      	  
	  NM <- Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$NM
  	NK <- Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$NK
	  NL <- Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$NL
  	zz <- 3.20 # Best standard score in Milligan and Cooper 1985
  	zzz <- zz*sqrt(2*(1-8/((pi^2)*pp))/(NM*pp))
    

    
	  if (any(indice == 14) || (indice == 31) || (indice == 32))
	  {
	     	resCritical[nc-min_nc+1,1] <- critValue <- 1-(2/(pi*pp))-zzz
	  }
    
	  if ((indice == 15)|| (indice == 31) || (indice == 32))
	  {
	      critValue <- 1-(2/(pi*pp))-zzz
	      resCritical[nc-min_nc+1,2] <- ((1-critValue)/critValue)*(NK+NL-2)
	    
	  }
    
    
	  if (any(indice == 16) || (indice == 31) || (indice == 32))
	  {
	  	df2 <- (NM-2)*pp
	  	resCritical[nc-min_nc+1,3] <- 1-pf(beale,pp,df2)
	  }
	}

  ########### Indices.TracesL-ball - 18e colonne de res ############
	if (any(indice == 18) || (indice == 31) || (indice == 32))
	{        	  
	  res[nc-min_nc+1,18] <- Indices.Traces(jeu, md, clall1, index = "ball")
	}

  ########### Indice.Point-Biserial - 19e colonne de res ############ 
	if (any(indice == 19) || (indice == 31) || (indice == 32))
	{       
	  res[nc-min_nc+1,19] <- Indice.ptbiserial(x=jeu, md=md, cl1=cl1)     
	}
   
  ########### Gap - 20e colonne de res ############       	
	if (any(indice == 20) || (indice == 32))
	{
	  
	  if (method == 1) {
		resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "ward.D2", d = NULL, centrotypes = "centroids")
		}
	  if (method == 2) {
	    	resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "single", d = NULL, centrotypes = "centroids")
		}
	  if (method == 3) {
	    	resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "complete", d = NULL, centrotypes = "centroids")
		}
	  if (method == 4) {
	    	resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "average", d = NULL, centrotypes = "centroids")
		}
	  if (method == 5) {
	    	resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "mcquitty", d = NULL, centrotypes = "centroids")
		}
	  if (method == 6) {
	    	resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "median", d = NULL, centrotypes = "centroids")
		}
	  if (method == 7) {
	    	resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "centroid", d = NULL, centrotypes = "centroids")
		}
		if (method == 9) {
		  resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "ward.D", d = NULL, centrotypes = "centroids")
		}
	  if (method == 8) {
		resultSGAP <- Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "k-means", d = NULL, centrotypes = "centroids")
		}
    	res[nc-min_nc+1,20] <- resultSGAP$gap
	resCritical[nc-min_nc+1,4] <- resultSGAP$diffu
	}

	if (nc >=2)
	{
	   ########### Indices.Traces-kl - 1e colonne de res ############ 	  
 	   if (any(indice == 1) || (indice == 31) || (indice == 32)) 
	   {	 
        res[nc-min_nc+1,1] <- Indices.Traces(jeu, md, clall1, index = "kl")
	   }

     ########### Indices.Traces-ch - 2e colonne de res ############
 	   if (any(indice == 2) || (indice == 31) || (indice == 32)) 
	   {		   
	      res[nc-min_nc+1,2] <- Indices.Traces(jeu, md, clall1, index = "ch")
	   }
	   
	   ########### Indice.cindex - 11e colonne de res ############
	   if (any(indice == 11) || (indice == 31) || (indice == 32)) 
	   {	  	   
	      res[nc-min_nc+1,11] <- Indice.cindex(d=md, cl=cl1)
	   }

     ########### Indice.DB  - 12e colonne de res ############
 	   if (any(indice == 12) || (indice == 31) || (indice == 32)) 
	   {		   
	      res[nc-min_nc+1,12] <- Indice.DB(x=jeu, cl=cl1, d = NULL, centrotypes = "centroids", p = 2, q = 2)$DB
	   }                         

     ########### Silhouette - 13e colonne de res ############
 	   if (any(indice == 13) || (indice == 31) || (indice == 32)) 
	   {		   
	      res[nc-min_nc+1,13] <- Indice.S(d=md, cl=cl1)
	   }
	   
	   ########### Indices.Traces-ratkowsky- 17e colonne de res ############
	   	if (any(indice == 17) || (indice == 31) || (indice == 32))
	   {  	  
	      res[nc-min_nc+1,17] <- Indices.Traces(jeu, md, clall1, index = "ratkowsky")
     }
     
     ########### Indice.Frey - 21e colonne de res ############
      if (any(indice == 21) || (indice == 31) || (indice == 32))
	   {      
        res[nc-min_nc+1,21] <- Index.15and28(cl1=cl1,cl2=cl2,md=md)$frey
	   }

     ########### Indice.McClain - 22e colonne de res ############
	   if (any(indice == 22) || (indice == 31) || (indice == 32))
	   {  	     
	     res[nc-min_nc+1,22] <- Index.15and28(cl1=cl1,cl2=cl2,md=md)$mcclain
	   }
	    
	   ########### Indice.Gamma - 23e colonne de res ############ 
	     if (any(indice == 23) || (indice == 32))
	   {            
       
	       res[nc-min_nc+1,23] <- Index.sPlussMoins(cl1=cl1,md=md)$gamma
	   }

     ########### Indice.Gplus- 24e colonne de res ############
	   if (any(indice == 24) || (indice == 32))
	   {    	     
	     res[nc-min_nc+1,24] <- Index.sPlussMoins(cl1=cl1,md=md)$gplus
	   }

     ########### Indice.Tau  - 25e colonne de res ############
	   if (any(indice == 25) || (indice == 32))
	   {   	     
	     res[nc-min_nc+1,25] <- Index.sPlussMoins(cl1=cl1,md=md)$tau
	   } 
     
     ########### Indices.Dunn  - 26e colonne de res ############
     if (any(indice == 26 ) || (indice == 31) || (indice == 32))
	   {    	    
	     res[nc-min_nc+1,26] <- Index.dunn(md, cl1, Data=jeu, method=NULL)	
	   } 
    
     ########### Indices.Hubert - 27e colonne de res ############
     if (any(indice == 27 ) || (indice == 31) || (indice == 32))
	   {         	     
	     res[nc-min_nc+1,27] <- Index.Hubert(jeu, cl1)	
	   }	 
	   
	   ########### Indices.SD - 28e colonne de res ############
     if (any(indice == 28 ) || (indice == 31) || (indice == 32))
	   {	    
	    res[nc-min_nc+1,28] <- Index.sdindex(jeu, clmax, cl1)
	   }	
	   
	   ########### Indices.Dindex - 29e colonne de res ############ 
	   	if (any(indice == 29 ) || (indice == 31) || (indice == 32))
	    {        	        
	        res[nc-min_nc+1,29] <- Index.Dindex(cl1, jeu)	
      }  
	   
	    ########### Indices.SDbw - 30e colonne de res ############
	   if (any(indice == 30 ) || (indice == 31) || (indice == 32))
	   { 	      
	       res[nc-min_nc+1,30] <- Index.SDbw(jeu, cl1)	
	   }      	
 	   
	}

	else
  {
	res[nc-min_nc+1,1] <- NA
	res[nc-min_nc+1,2] <- NA
	res[nc-min_nc+1,11] <- NA
	res[nc-min_nc+1,12] <- NA
	res[nc-min_nc+1,13] <- NA
	res[nc-min_nc+1,17] <- NA
	res[nc-min_nc+1,21] <- NA
	res[nc-min_nc+1,22] <- NA
	res[nc-min_nc+1,23] <- NA
	res[nc-min_nc+1,24] <- NA
	res[nc-min_nc+1,25] <- NA
	res[nc-min_nc+1,26] <- NA
	res[nc-min_nc+1,27] <- NA
	res[nc-min_nc+1,28] <- NA
	res[nc-min_nc+1,29] <- NA
	res[nc-min_nc+1,30] <- NA
	}
}

#########################################################################################################
#######################                      Best Number of Clusters                  ###################
#########################################################################################################



   nc.KL<-indice.KL<-0
   if (any(indice == 1) || (indice == 31) || (indice == 32)) 
	 {  
     # KL - The value of u, which maximizes KL(u), is regarded as specifying the number of clusters [ClusterSim package].
     nc.KL <- (min_nc:max_nc)[which.max(res[,1])]
     indice.KL <- max(res[,1],na.rm = TRUE)
     best.nc<-nc.KL
   }
  
   nc.CH<-indice.CH<-0
   if (any(indice == 2) || (indice == 31) || (indice == 32)) 
	 {
     # CH - The value of u, which maximizes CH(u), is regarded as specifying the number of clusters [ClusterSim package].
     nc.CH <- (min_nc:max_nc)[which.max(res[,2])]
     indice.CH <- max(res[,2],na.rm = TRUE)
     best.nc<-nc.CH
   }
  
  nc.CCC<-indice.CCC<-0
  if (any(indice == 4) || (indice == 31) || (indice == 32))
	{
    # CCC - The maximum value accross the hierarchy levels is used to indicate the optimal number of clusters in data [29].
    nc.CCC <- (min_nc:max_nc)[which.max(res[,4])]
    indice.CCC <- max(res[,4],na.rm = TRUE)
    best.nc<-nc.CCC
  }
    
  nc.DB<-indice.DB<-0 
  if (any(indice == 12) || (indice == 31) || (indice == 32)) 
	{
    # DB - The value of u, which minimizes DB(u), is regarded as specifying the number of clusters [clusterSim package].
    nc.DB <- (min_nc:max_nc)[which.min(res[,12])]
    indice.DB <- min(res[,12],na.rm = TRUE)
    best.nc<-nc.DB
  }
  
  nc.Silhouette<-indice.Silhouette<-0
  if (any(indice == 13) || (indice == 31) || (indice == 32)) 
	{
    # SILHOUETTE - The value of u, which maximizes S(u), is regarded as specifying the number of clusters [ClusterSim package].
    nc.Silhouette <- (min_nc:max_nc)[which.max(res[,13])]
    indice.Silhouette <- max(res[,13],na.rm = TRUE)
    best.nc<-nc.Silhouette
  }

  nc.Gap<-indice.Gap<-0
  # GAP - Choose the number of clusters via finding the smallest q such that: Gap(q)=Gap(q+1)-Sq+1 (q=1,\u{85},n-2) [ClusterSim package].
  if (any(indice == 20) || (indice == 32))
  {
    found <- FALSE
    for (ncG in min_nc:max_nc){
      if ((resCritical[ncG-min_nc+1,4] >=0) && (!found)){
          ncGap <- ncG
     	  indiceGap <- res[ncG-min_nc+1,20]
    	  found <- TRUE
    	  }
     }
     if (found){
  	 nc.Gap <- ncGap
  	 indice.Gap <- indiceGap
     best.nc<-nc.Gap
      }else{
  	    nc.Gap <- NA
  	    indice.Gap <- NA
      }

  }

  nc.Duda<-indice.Duda<-0
  # DUDA - Choose the number of clusters via finding the smallest q such that: duda >= critical_value [Duda and Hart (1973)].
 
  
   if (any(indice == 14) || (indice == 31) || (indice == 32))
	{
        
        foundDuda <- FALSE
        for (ncD in min_nc:max_nc)
        {
           
           if ((res[ncD-min_nc+1,14]>=resCritical[ncD-min_nc+1,1]) && (!foundDuda))
           {
             ncDuda <- ncD
     	       indiceDuda <- res[ncD-min_nc+1,14]
    	       foundDuda <- TRUE
    	     }
        }
        if (foundDuda)
        {
  	      nc.Duda <- ncDuda
  	      indice.Duda <- indiceDuda
          best.nc<-nc.Duda
        }
        else
        {
  	        nc.Duda <- NA
  	        indice.Duda <- NA
        }
       
     
   }  
    
  nc.Pseudo<-indice.Pseudo<-0  
  # PSEUDOT2 - Chooses the number of clusters via finding the smallest q such that: pseudot2 <= critical_value [SAS User's guide].
	if (any(indice == 15) || (indice == 31) || (indice == 32))
	{
     
     foundPseudo <- FALSE
     for (ncP in min_nc:max_nc)
       {

      if ((res[ncP-min_nc+1,15]<=resCritical[ncP-min_nc+1,2]) && (!foundPseudo))
        {
          ncPseudo <- ncP
     	  indicePseudo <- res[ncP-min_nc+1,15]
    	  foundPseudo <- TRUE
    	  }
     }
      if (foundPseudo)
        {
  	 nc.Pseudo <- ncPseudo
  	 indice.Pseudo <- indicePseudo
     best.nc<-nc.Pseudo
      }
     else
        {
  	    nc.Pseudo <- NA
  	    indice.Pseudo <- NA
      }
    }
  
  
  nc.Beale<-indice.Beale<-0
  	if (any(indice == 16) || (indice == 31) || (indice == 32))
	{
  # BEALE - Chooses the number of clusters via finding the smallest q such that: Fvalue_beale >= 0.1 [Gordon (1999)].
     foundBeale <- FALSE
     for (ncB in min_nc:max_nc)
       {
     
      if ((resCritical[ncB-min_nc+1,3]>=alphaBeale) && (!foundBeale)){
          ncBeale <- ncB
     	  indiceBeale <- res[ncB-min_nc+1,16]
    	  foundBeale <- TRUE
    	  }
     }
       if (foundBeale){
  	 nc.Beale <- ncBeale
  	 indice.Beale <- indiceBeale
     best.nc<-nc.Beale
      }
     else
       {
  	    nc.Beale <- NA
  	    indice.Beale <- NA
      }
  }
 
  
  nc.ptbiserial<-indice.ptbiserial<-0
  if (any(indice == 19) || (indice == 31) || (indice == 32))
	{
    # POINT-BISERIAL - The maximum value was used to suggest the optimal number of clusters in the data [29].
    nc.ptbiserial <- (min_nc:max_nc)[which.max(res[,19])]
    indice.ptbiserial <- max(res[,19],na.rm = TRUE)
    best.nc<-nc.ptbiserial
  }

   foundNC<-foundIndice<-numeric(0)
   nc.Frey<-indice.Frey<-0
   if (any(indice == 21) || (indice == 31) || (indice == 32))
	 {
  # FREY AND VAN GROENEWOUD - The best results occured when clustering was continued until the ratio fell below 1.00 for the last 
  #			      series of times. At this point, the cluster level before this series was taken as the optimal partition. 
  #			      If the ration never fell below 1.00, a one cluster solution was assumed [29].
  
     foundFrey <- FALSE
     i<-1
     for (ncF in min_nc:max_nc)
     {          
             
          if (res[ncF-min_nc+1,21] < 1) 
          {
                       
              ncFrey <- ncF-1               
     	        indiceFrey <- res[ncF-1-min_nc+1,21]
     	        foundFrey <- TRUE
              foundNC[i]<-ncFrey
              foundIndice[i]<-indiceFrey
              i<-i+1
  	        
    	     }
       
     }
     if (foundFrey)
     {
  	   nc.Frey <- foundNC[1]
  	   indice.Frey <- foundIndice[1]
       best.nc<-nc.Frey
     }
      else 
      {
  	    nc.Frey <- NA
  	    indice.Frey <- NA
  	    print(paste("Frey index : No clustering structure in this data set"))
      }
     
      
   }  
      
      
   nc.McClain<-indice.McClain<-0
   if (any(indice == 22) || (indice == 31) || (indice == 32))
	{
  # MCCLAIN AND RAO - The minimum value of the index was found to give the best recovery information [29].
  nc.McClain <- (min_nc:max_nc)[which.min(res[,22])]
  indice.McClain <- min(res[,22],na.rm = TRUE)
  best.nc<-nc.McClain
  
  }
  
   nc.Gamma<-indice.Gamma<-0
   if (any(indice == 23) || (indice == 31) || (indice == 32))
	{
       # GAMMA - Maximum values were taken to represent the correct hierarchy level [29].
       nc.Gamma <- (min_nc:max_nc)[which.max(res[,23])]
       indice.Gamma <- max(res[,23],na.rm = TRUE)
       best.nc<-nc.Gamma
       
  }

   nc.Gplus<-indice.Gplus<-0
   if (any(indice == 24) || (indice == 31) || (indice == 32))
	{
  # GPLUS - Minimum values were used to determine the number of clusters in the data [29].
  nc.Gplus  <- (min_nc:max_nc)[which.min(res[,24])]
  indice.Gplus <- min(res[,24],na.rm = TRUE)
  best.nc<-nc.Gplus
  }

   nc.Tau<-indice.Tau<-0
   if (any(indice == 25) || (indice == 31) || (indice == 32))
	{
  # TAU - The maximum value in the hierarchy sequence was taken as indicating the correct number of clusters [29].
  nc.Tau <- (min_nc:max_nc)[which.max(res[,25])]
  indice.Tau <- max(res[,25],na.rm = TRUE)
  best.nc<-nc.Tau
  }
     
  
#Some indices need to compute difference between hierarchy levels to identify relevant number of clusters
  
 
  if((indice==3)||(indice == 5)||(indice == 6)||(indice == 7)||(indice == 8)||(indice == 9)||(indice == 10)||(indice == 18)||(indice == 27)||(indice == 11)||(indice == 29)||(indice == 31)||(indice == 32))
  {
   
    DiffLev <- array(0, c(max_nc-min_nc+1,12))
    DiffLev[,1] <- min_nc:max_nc
       for (nc3 in min_nc:max_nc)      
      {
        if (nc3==min_nc)
        {    
	       DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-NA)   # Hartigan
	       DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-NA)   #Scott
	       DiffLev[nc3-min_nc+1,4] <- abs(res[nc3-min_nc+1,6]-NA)   # Marriot
	       DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-NA)   #Trcovw
	       DiffLev[nc3-min_nc+1,6] <- abs(res[nc3-min_nc+1,8]-NA)   #Tracew
	       DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-NA)   #Friedman
	       DiffLev[nc3-min_nc+1,8] <- abs(res[nc3-min_nc+1,10]-NA)  #Rubin
	       DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-NA)  # Ball
         DiffLev[nc3-min_nc+1,10] <- abs(res[nc3-min_nc+1,27]-NA) # Hubert   
         DiffLev[nc3-min_nc+1,12] <- abs(res[nc3-min_nc+1,29]-NA) # D index
         
         
	      }
        else
        {  
          if(nc3==max_nc)
          { 
            DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-res[nc3-min_nc,3])
            DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-res[nc3-min_nc,5])
            DiffLev[nc3-min_nc+1,4] <- abs(res[nc3-min_nc+1,6]-NA) # Marriot
            DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-res[nc3-min_nc,7])  # trcovw
            DiffLev[nc3-min_nc+1,6] <- abs(res[nc3-min_nc+1,8]-NA) #traceW            
            DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-res[nc3-min_nc,9])
            DiffLev[nc3-min_nc+1,8] <- abs(res[nc3-min_nc+1,10]-NA) #Rubin
            DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-res[nc3-min_nc,18])
            DiffLev[nc3-min_nc+1,10] <- abs(res[nc3-min_nc+1,27]-NA)
            DiffLev[nc3-min_nc+1,12] <- abs(res[nc3-min_nc+1,29]-NA) # D index  

	 
	         }
          else      
          {
        
           DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-res[nc3-min_nc,3]) # Hartigan              
	         DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-res[nc3-min_nc,5]) 
           DiffLev[nc3-min_nc+1,4] <- ((res[nc3-min_nc+2,6]-res[nc3-min_nc+1,6])-(res[nc3-min_nc+1,6]-res[nc3-min_nc,6]))
           DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-res[nc3-min_nc,7])
           DiffLev[nc3-min_nc+1,6] <- ((res[nc3-min_nc+2,8]-res[nc3-min_nc+1,8])-(res[nc3-min_nc+1,8]-res[nc3-min_nc,8]))
           DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-res[nc3-min_nc,9])
           DiffLev[nc3-min_nc+1,8] <- ((res[nc3-min_nc+2,10]-res[nc3-min_nc+1,10])-(res[nc3-min_nc+1,10]-res[nc3-min_nc,10]))
           DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-res[nc3-min_nc,18])  
           DiffLev[nc3-min_nc+1,10] <- abs((res[nc3-min_nc+1,27]-res[nc3-min_nc,27]))             
           DiffLev[nc3-min_nc+1,12] <-((res[nc3-min_nc+2,29]-res[nc3-min_nc+1,29])-(res[nc3-min_nc+1,29]-res[nc3-min_nc,29])) #Dindex     
          
          }         
       }         
    }
   }

  nc.Hartigan<-indice.Hartigan<-0
  if (any(indice == 3) || (indice == 31) || (indice == 32))
	{
	# HARTIGAN - The maximum differences between hierarchy levels were taken as indicating the correct number of clusters in the data [29].
	 nc.Hartigan <- DiffLev[,1][which.max(DiffLev[,2])]
	 indice.Hartigan <- max(DiffLev[,2],na.rm = TRUE)
   best.nc<-nc.Hartigan
  }
  
  nc.Ratkowsky<-indice.Ratkowsky<-0
  if (any(indice == 17) || (indice == 31) || (indice == 32))
	{
  # RATKOWSKY - The optimal number of groups is taken as the level where this criterion exhibits its maximum value [29].
    nc.Ratkowsky <- (min_nc:max_nc)[which.max(res[,17])]
    indice.Ratkowsky <- max(res[,17],na.rm = TRUE)
    best.nc<-nc.Ratkowsky
  }
    
    nc.cindex<-indice.cindex<-0
    if (any(indice == 11) || (indice == 31) || (indice == 32)) 
    {
      #CINDEX - The minimum value across the hierarchy levels was used to indicate the optimal number of clusters [29].
      nc.cindex <- (min_nc:max_nc)[which.min(res[,11])]
      indice.cindex <- min(res[,11],na.rm = TRUE)
      best.nc<-nc.cindex
    }  
  
  nc.Scott<-indice.Scott<-0
  if (any(indice == 5) || (indice == 31) || (indice == 32))
	{
 	# SCOTT - The maximum difference between hierarchy levels was used to suggest the correct number of partitions [29].
	 nc.Scott <- DiffLev[,1][which.max(DiffLev[,3])]
	 indice.Scott <- max(DiffLev[,3],na.rm = TRUE)
   best.nc<-nc.Scott
  }
  
  nc.Marriot<-indice.Marriot<-0
  if (any(indice == 6) || (indice == 31) || (indice == 32))
	{
	# MARRIOT - The maximum difference between successive levels was used to determine the best partition level [29].
	 nc.Marriot <- DiffLev[,1][which.max(DiffLev[,4])]
	 round(nc.Marriot, digits=1)
	 indice.Marriot <- max(DiffLev[,4],na.rm = TRUE)
   best.nc<-nc.Marriot
  }
  
  nc.TrCovW<-indice.TrCovW<-0
  if (any(indice == 7) || (indice == 31) || (indice == 32))
	{
	nc.TrCovW <- DiffLev[,1][which.max(DiffLev[,5])]
	indice.TrCovW <- max(DiffLev[,5],na.rm = TRUE)
	best.nc<-nc.TrCovW
  }
  
  
  nc.TraceW<-indice.TraceW<-0
  if (any(indice == 8) || (indice == 31) || (indice == 32))
	{
  	# TRACE W - To determine the number of clusters in the data, maximum difference scores were used [29].
	  nc.TraceW <- DiffLev[,1][which.max(DiffLev[,6])]
  	indice.TraceW <- max(DiffLev[,6],na.rm = TRUE)
	  best.nc<-nc.TraceW
   }
   
  nc.Friedman<-indice.Friedman<-0
  if (any(indice == 9) || (indice == 31) || (indice == 32))
	{
  	# FRIEDMAN - The maximum difference in values of trace W-1B criterion was used to indicate the optimal number of clusters [29].
  	nc.Friedman <- DiffLev[,1][which.max(DiffLev[,7])]
  	indice.Friedman <- max(DiffLev[,7],na.rm = TRUE)
  	best.nc<-nc.Friedman
	}
	
	nc.Rubin<-indice.Rubin<-0
  if (any(indice == 10) || (indice == 31) || (indice == 32))
	{
	  # RUBIN - The difference between levels was used [29].
	  nc.Rubin <- DiffLev[,1][which.min(DiffLev[,8])]
	  indice.Rubin <- min(DiffLev[,8],na.rm = TRUE)
	  best.nc<-nc.Rubin
  }
  
  nc.Ball<-indice.Ball<-0
  if (any(indice == 18) || (indice == 31) || (indice == 32))
	{
	  # BALL - The largest difference between levels was used to indicate the optimal solution [29].
	  nc.Ball <- DiffLev[,1][which.max(DiffLev[,9])]
	  indice.Ball <- max(DiffLev[,9],na.rm = TRUE)
	  best.nc<-nc.Ball
  }
  

   nc.Dunn<-indice.Dunn<-0
   if (any(indice == 26) || (indice == 31) || (indice == 32)) 
	 {
     # Dunn - 
     nc.Dunn <- (min_nc:max_nc)[which.max(res[,26])]
     indice.Dunn <- max(res[,26],na.rm = TRUE)
     best.nc<-nc.Dunn
   }
  
  
   nc.Hubert<-indice.Hubert<-0
   if (any(indice == 27) || (indice == 31) || (indice == 32)) 
	 {       
	   # Hubert - 
     nc.Hubert  <- 0.00
     indice.Hubert  <- 0.00
     #x11()
     par(mfrow = c(1,2))
     plot(x_axis,res[,27], tck=0, type="b", col="red", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Hubert Statistic values")))
     plot(DiffLev[,1],DiffLev[,10], tck=0, type="b", col="blue", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Hubert statistic second differences")))
     cat(paste ("*** : The Hubert index is a graphical method of determining the number of clusters.
                In the plot of Hubert index, we seek a significant knee that corresponds to a 
                significant increase of the value of the measure i.e the significant peak in Hubert
                index second differences plot.", "\n", "\n"))
     }
  
   nc.sdindex<-indice.sdindex<-0
   if (any(indice == 28) || (indice == 31) || (indice == 32)) 
	 {
     # SD - 
     nc.sdindex <- (min_nc:max_nc)[which.min(res[,28])]
     indice.sdindex<- min(res[,28],na.rm = TRUE)
     best.nc<-nc.sdindex
   }
  
    
    nc.Dindex<-indice.Dindex<-0
    if (any(indice == 29) || (indice == 31) || (indice == 32)) 
	  {

     nc.Dindex <- 0.00
     indice.Dindex<- 0.00
     #x11()
     par(mfrow = c(1,2))
     plot(x_axis,res[,29], tck=0, type="b", col="red", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Dindex Values")))
     plot(DiffLev[,1],DiffLev[,12], tck=0, type="b", col="blue", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Second differences Dindex Values")))
     cat(paste ("*** : The D index is a graphical method of determining the number of clusters. 
                In the plot of D index, we seek a significant knee (the significant peak in Dindex
                second differences plot) that corresponds to a significant increase of the value of
                the measure.", "\n", "\n"))
    }
  
    nc.SDbw<-indice.SDbw<-0
    if (any(indice == 30) || (indice == 31) || (indice == 32)) 
	   {
       # SDbw - 
       nc.SDbw <- (min_nc:max_nc)[which.min(res[,30])]
       indice.SDbw<- min(res[,30],na.rm = TRUE)  
       best.nc<-nc.SDbw
     } 
  
  

######################################################################################################################
########################                Displaying results             #########################################
######################################################################################################################
    
 if (indice < 31)
 {
     res <- res[,c(indice)]
        
     if (indice == 14) { resCritical <- resCritical[,1]  }
     if (indice == 15) { resCritical <- resCritical[,2] }
     if (indice == 16) { resCritical <- resCritical[,3] }
     if (indice == 20) { resCritical <- resCritical[,4] }        
   
 }

 if (indice == 31) 
  { 
      res <- res[,c(1:19,21:22,26:30)]
		  resCritical <- resCritical[,c(1:3)]        
      		  
  }



 if (any(indice == 20) || (indice == 23) || (indice == 24) || (indice == 25) || (indice == 32))
 {
  
  results <- c(nc.KL, indice.KL, nc.CH, indice.CH, nc.Hartigan, indice.Hartigan, nc.CCC, indice.CCC, nc.Scott, indice.Scott,
		nc.Marriot, indice.Marriot, nc.TrCovW, indice.TrCovW, nc.TraceW, indice.TraceW, nc.Friedman, 
		indice.Friedman, nc.Rubin, indice.Rubin, nc.cindex, indice.cindex, nc.DB, indice.DB, nc.Silhouette,
		indice.Silhouette, nc.Duda, indice.Duda, nc.Pseudo, indice.Pseudo, nc.Beale, indice.Beale, nc.Ratkowsky,
		indice.Ratkowsky, nc.Ball, indice.Ball, nc.ptbiserial, indice.ptbiserial, nc.Gap, indice.Gap, 
		nc.Frey, indice.Frey, nc.McClain, indice.McClain, nc.Gamma, indice.Gamma, nc.Gplus, indice.Gplus,
		nc.Tau, indice.Tau, nc.Dunn, indice.Dunn, nc.Hubert, indice.Hubert, nc.sdindex, indice.sdindex, nc.Dindex, indice.Dindex, nc.SDbw, indice.SDbw)
  results1<-matrix(c(results),nrow=2,ncol=30)
  resultats <- matrix(c(results),nrow=2,ncol=30,dimnames = list(c("Number_clusters","Value_Index"),
		     c("KL","CH","Hartigan","CCC", "Scott", "Marriot", "TrCovW",
		       "TraceW", "Friedman", "Rubin", "Cindex", "DB", "Silhouette",
			"Duda","PseudoT2", "Beale", "Ratkowsky", "Ball", "PtBiserial",
			"Gap", "Frey", "McClain", "Gamma", "Gplus", "Tau", "Dunn", 
      "Hubert", "SDindex", "Dindex", "SDbw")))
   }
 else
  {
    
    results <- c(nc.KL, indice.KL, nc.CH, indice.CH, nc.Hartigan, indice.Hartigan, nc.CCC, indice.CCC, nc.Scott, indice.Scott,
		nc.Marriot, indice.Marriot, nc.TrCovW, indice.TrCovW, nc.TraceW, indice.TraceW, nc.Friedman, indice.Friedman, 
    nc.Rubin, indice.Rubin, nc.cindex, indice.cindex, nc.DB, indice.DB, nc.Silhouette, indice.Silhouette,
    nc.Duda, indice.Duda, nc.Pseudo, indice.Pseudo, nc.Beale, indice.Beale, nc.Ratkowsky, indice.Ratkowsky,
    nc.Ball, indice.Ball, nc.ptbiserial, indice.ptbiserial, nc.Frey, indice.Frey, nc.McClain, indice.McClain, 
    nc.Dunn, indice.Dunn, nc.Hubert, indice.Hubert, nc.sdindex, indice.sdindex, nc.Dindex, indice.Dindex, nc.SDbw, indice.SDbw 
    )
    results1<-matrix(c(results),nrow=2,ncol=26)
    resultats <- matrix(c(results),nrow=2,ncol=26,dimnames = list(c("Number_clusters","Value_Index"),
		c("KL","CH","Hartigan","CCC", "Scott", "Marriot", "TrCovW",
		"TraceW", "Friedman", "Rubin", "Cindex", "DB", "Silhouette",
		 "Duda","PseudoT2", "Beale", "Ratkowsky", "Ball", "PtBiserial",
		"Frey", "McClain", "Dunn", 		"Hubert", "SDindex", "Dindex", "SDbw")))
   
   }
 
   
 if (any(indice <= 20)||(indice == 23)||(indice == 24)||(indice == 25)) 
 {   
   resultats <- resultats[,c(indice)]   
 }
 
 if (any(indice == 21)|| (indice == 22)) 
 {  
  indice3 <-indice-1
  resultats <- resultats[,c(indice3)]   
 }
 
 if (any(indice == 26) || (indice == 27) || (indice == 28) || ( indice == 29)|| ( indice == 30)) 
 { 
  indice4 <- indice-4     
  resultats <- resultats[,c(indice4)] 
 }
 
    
  resultats<-round(resultats, digits=4)
  res<-round(res, digits=4)
  resCritical<-round(resCritical, digits=4)

#  if (numberObsAfter != numberObsBefore) 
#  {
#     cat(paste(numberObsAfter,"observations were used out of", numberObsBefore ,"possible observations due to missing values."))
#  }
  
#  if (numberObsAfter == numberObsBefore) 
#  {
#     cat(paste("All", numberObsAfter,"observations were used.", "\n", "\n"))
#  }
  
  
    
    ######################## Summary results #####################################
    
    
    if(any(indice == 31) || (indice == 32))
    {
      cat("*******************************************************************", "\n")
      cat("* Among all indices:                                               ", "\n")
      BestCluster<-results1[1,]
      c=0
      for(i in min.nc:max.nc)
      {
        vect<-which(BestCluster==i)
        if(length(vect)>0)
        cat("*",length(vect), "proposed", i,"as the best number of clusters", "\n")
      
        if(c<length(vect))
        { 
          j=i 
          c<-length(vect)
        }
      }
    
        cat("\n","                  ***** Conclusion *****                           ", "\n", "\n")
        cat("* According to the majority rule, the best number of clusters is ",j , "\n", "\n", "\n")
        cat("*******************************************************************", "\n")
      
        
      ########################## The Best partition    ###################
    
        if (any(method == 1) || (method == 2) || (method == 3) || (method == 4) || 
          (method == 5) || (method == 6) || (method == 7)||(method == 9))         
            partition<- cutree(hc, k=j)
    
        else
        {
            set.seed(1)
            partition<-kmeans(jeu,j)$cluster
        }
    
    }
    
    
    if (any(indice==1)||(indice==2)||(indice==3)||(indice==4)||(indice==5)||(indice==6)||(indice==7)
        ||(indice==8)||(indice==9)||(indice==10)||(indice==11)||(indice==12)||(indice==13)||(indice==14)
        ||(indice==15)||(indice==16)||(indice==17)||(indice==18)||(indice==19)||(indice==20)
        ||(indice==21)||(indice==22)||(indice==23)||(indice==24)||(indice==25)||(indice==26)
        ||(indice==28)||(indice==30))
    {
      if (any(method == 1) || (method == 2) || (method == 3) || (method == 4) || 
            (method == 5) || (method == 6) || (method == 7) || (method == 9)) 
      
        partition<- cutree(hc, k=best.nc)
      
      else
      {
        set.seed(1)
        partition<-kmeans(jeu,best.nc)$cluster
      }
            
    }
      
    
    #########################  Summary results   ############################
    
    
    
    if ((indice == 14)|| (indice == 15)|| (indice == 16)|| (indice == 20)|| (indice == 31)|| (indice == 32))
    { 
      results.final <- list(All.index=res,All.CriticalValues=resCritical,Best.nc=resultats, Best.partition=partition)
    }
    
    if ((indice == 27)|| (indice == 29))
       results.final <- list(All.index=res)
    
    if (any(indice==1)||(indice==2)||(indice==3)||(indice==4)||(indice==5)||(indice==6)||(indice==7)
        ||(indice==8)||(indice==9)||(indice==10)||(indice==11)||(indice==12)||(indice==13)
        ||(indice==17)||(indice==18)||(indice==19)||(indice==21)||(indice==22)||(indice==23)||(indice==24)
        ||(indice==25)||(indice==26)||(indice==28)||(indice==30))  
         
      results.final <- list(All.index=res,Best.nc=resultats, Best.partition=partition)
    
       
      
    return(results.final)
    
   
}
