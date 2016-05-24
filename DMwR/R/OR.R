# ======================================================================
# OUTLIER RANKING
#
# (c) Luis Torgo, 2004-2008
# ======================================================================
# Using hierarchical clustering to obtain a ranking of outlierness for a
# set of "test" cases, using also a set of "tranining" cases.
# This separation in test and train is artificial and only necessary if the
# test cases are few...
# If the test.data is NULL the ranking respects all data given in the first
# parameter (the data), otherwise the obtained ranking respects the test.data
# The user can also supply a distance matrix instead in the data parameter.
# In that case (only usefull for optimization issues when trying different
# clustering algorithms for the same distance function) the ranking is the
# same as if test.data was NULL.
#
# The ranking can be obtained by different methods:
#
#
# ----------------------------------------------------------------------
outliers.ranking <- function(data,test.data=NULL,
                             method='sizeDiff',
                             method.pars=NULL,
                             clus=list(dist='euclidean',alg='hclust',
                               meth='ward'),
                             power=1,
                             verb=F) {



  ##require(cluster)
  

  # -------- Check the type of task
  
  # T0: the rankings are going to be on the "data" but the distance matrix is given
  if (inherits(data,'dist')) {
    test.pos <- 1
    N <- (1+sqrt(1+4*length(data)*2))/2
  # T1: the rankings are going to be on the "training data"
  } else if (is.null(test.data) ) {
    N <- NROW(data)
    test.pos <- 1
    
  # T2: the rankings are for the "test data"
  } else {
    train.sz <- NROW(data)
    test.sz <- NROW(test.data)
    data <- rbind(data,test.data)
    N <- train.sz+test.sz
    test.pos <- train.sz+1
  }


  
  if (verb) cat('OR:: Distance calculation...')
  # ------- Distance Calculation
  if (!inherits(data,'dist')) { # we were given a data set
    if (power > 1) 
      dist.mtrx <- dist(data,method=clus$dist)^power
    else
      dist.mtrx <- dist(data,method=clus$dist)
  } else dist.mtrx <- data    # we were already given a distance matrix

  

  if (verb) cat('\nOR:: Clustering...')
  # ------- Hierarchical Clustering
  if (clus$alg != 'diana') {
    h <- do.call(clus$alg,list(dist.mtrx,method=clus$meth))
  } else {
    h <- do.call(clus$alg,list(dist.mtrx))
  }


  
  if (verb) cat('\nOR:: Ranking...')
  # ------- Ranking
  # This is the major step, obtain rankings based on clustering results

  # This vector will hold the ranking score of each data point
  rk <- rep(0,N)


  # This is the linear method
  if (method=='linear') {
    #out.lim <- method.pars$sz.perc*N
    out.lim <- max(round(method.pars$sz.perc*N,0),2)
    szs <- rep(0,NROW(h$merge))
    mb <- list()
    for(ln in 1:length(szs)) {
      x <- sapply(1:2,function(p) ifelse(h$merge[ln,p] < 0,1,szs[h$merge[ln,p]]))
      szs[ln] <- sum(x)
      mb[[ln]] <- c(g1 <- if (h$merge[ln,1] < 0) -h$merge[ln,1] else mb[[h$merge[ln,1]]],
                    g2 <- if (h$merge[ln,2] < 0) -h$merge[ln,2] else mb[[h$merge[ln,2]]])

      h.fac <- ln/(N-1)
      
      if (x[1] < x[2]) {
        sc <- if (x[1] > out.lim) 0 else (1-(x[1]-1)/(N-2))*h.fac
        rk[g1] <- ifelse(rk[g1] > sc,rk[g1],sc)
      } else {
        sc <- if (x[2] > out.lim) 0 else (1-(x[2]-1)/(N-2))*h.fac
        rk[g2] <- ifelse(rk[g2] > sc,rk[g2],sc)
      }
    }

    rk.outliers <- order(rk[test.pos:N],decreasing=T)
    pb.outliers <- rk[test.pos:N]



  # This is the sigmoidal
  } else  if (method=='sigmoid') {
    #out.lim <- method.pars$sz.perc*N
    out.lim <- 2*max(round(method.pars$sz.perc*N,0),2)
    #cTemp <- (7*out.lim^2)/(1-out.lim)^2 
    szs <- rep(0,NROW(h$merge))
    mb <- list()
    for(ln in 1:length(szs)) {
      x <- sapply(1:2,function(p) ifelse(h$merge[ln,p] < 0,1,szs[h$merge[ln,p]]))
      szs[ln] <- sum(x)
      mb[[ln]] <- c(g1 <- if (h$merge[ln,1] < 0) -h$merge[ln,1] else mb[[h$merge[ln,1]]],
                    g2 <- if (h$merge[ln,2] < 0) -h$merge[ln,2] else mb[[h$merge[ln,2]]])

      h.fac <- exp( (-2)*( ((ln-(N-1))^2) / ((N-1)^2) ) )
      if (x[1] < x[2]) {
        sc <- (x[1]< out.lim)* (1 - exp( (-4) * ((x[1]-out.lim)^2 / (out.lim)^2 ))) * h.fac
        rk[g1] <- ifelse(rk[g1] > sc,rk[g1],sc)
      } else {
        sc <- (x[2]< out.lim)* (1 - exp( (-4) * ((x[2]-out.lim)^2 / (out.lim)^2 ))) * h.fac
        rk[g2] <- ifelse(rk[g2] > sc,rk[g2],sc)
      }
    }

    rk.outliers <- order(rk[test.pos:N],decreasing=T)
    pb.outliers <- rk[test.pos:N]
    

  # This is the group size differences method
  } else if (method=='sizeDiff') {

    szs <- rep(0,NROW(h$merge))
    mb <- list()
    for(ln in 1:length(szs)) {
      x <- sapply(1:2,function(p) ifelse(h$merge[ln,p] < 0,1,szs[h$merge[ln,p]]))
      szs[ln] <- sum(x)
      mb[[ln]] <- c(g1 <- if (h$merge[ln,1] < 0) -h$merge[ln,1] else mb[[h$merge[ln,1]]],
                    g2 <- if (h$merge[ln,2] < 0) -h$merge[ln,2] else mb[[h$merge[ln,2]]])


      if (x[1] < x[2]) {
        sc <- (x[2]-x[1])/(x[1]+x[2])
        rk[g1] <- ifelse(rk[g1] > sc,rk[g1],sc)
      } else {
        sc <- (x[1]-x[2])/(x[1]+x[2])
        rk[g2] <- ifelse(rk[g2] > sc,rk[g2],sc)
      }
    }

    rk.outliers <- order(rk[test.pos:N],decreasing=T)
    pb.outliers <- rk[test.pos:N]
 
  }
  names(rk.outliers) <- names(pb.outliers) <- row.names(data)[test.pos:N]
  

  # ---- Now build up the list that will be returned
  if (verb) cat('\n')
  list(rank.outliers=rk.outliers,
       prob.outliers=pb.outliers, # this is in the natural order and not
       hie=h,                     # outlierness order. To get the latter do
       dist=dist.mtrx)            # res$prob.outliers[res$rank.outliers]
}




