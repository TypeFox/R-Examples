
# Mean shift base function
meanshift<-function(X, x, h){
   g1 <- function(xi, x, h){ # 1-d profile
     1/2*exp(-1/2*((x-xi)/h)^2)
   }
   gd <- function(Xi,x,h){# Multi-d profile
     d<-length(x)
     k<-1
     for (j in 1:d){
        k<- k* g1(Xi[,j],x[j],h[j])
     }
     k
   }
   d <- dim(X)[2]
   if (length(h) == 1) {
        h <- rep(h, d)
    }
   if(is.vector(X)){
       X<-matrix(X,nrow=length(X))
   }
   x <- as.numeric(x)
   g <- gd(X, x, h)
   ms <- NULL
   for (j in 1:d){
                 ms[j]<-sum(X[,j]*g)/sum(g)
                 }
   ms
   }

# Mean shift iterative function (until convergence ...)
ms.rep <- function (X, x, h, plotms=1, thresh= 0.00000001, iter=200) {    
          s  <- 0
          th <- rep(0,iter)
          M  <-matrix(0, iter, length(x))
          x0 <- x
          d  <- dim(X)[2]
          if (length(h) == 1) {
              h <- rep(h, d)
          }
         
          for (j in 1: iter){
           
            m     <- meanshift(X, x, h)        
            M[j,] <- m
            th[j] <- enorm(m-x)/enorm(x)          
            if (th[j]<thresh){
                  s<-j;                 
                  break
                }
            x     <- m
            if(plotms==1 && dim(X)[2]==2){points(x[1],x[2], col=1, pch=15)}
          }
          if (plotms==2 && dim(X)[2]==2){lines(rbind(x0,M[1:s,]))}
          #print(s)
       return(list("Meanshift.points"=M[1:s,], "Threshold.values"=th[1:s], "iterations"=s, "start"=as.numeric(x0),  "final"=m))
       }

# Mean shift clustering

ms<-function (X, h, subset, thr = 0.0001, scaled = TRUE, iter=200, plotms = 2, 
    or.labels = NULL, ...) 
{
    n <- dim(X)[1]
    d <- dim(X)[2]
    if (missing(subset)) {
        subset <- 1:n
    }
   s1       <- apply(X, 2, function(dat){ diff(range(dat))})  # range 
   if (missing(h)){
        if (!scaled){h   <- s1/10 } else { h<- 0.1}
   } # bandwidth by default: 10% of range
   if (length(h)==1){h <- rep(h,d)}
   if (scaled){        # scales the data to lie by its range
        X <- sweep(X, 2, s1, "/")
   } 

  
   if (d == 2 && plotms > 0) {
        if (missing(or.labels)) {            
            plot(X, col = "grey70", ...)
        }
        else {          
            plot(X, col = or.labels, ...)
        }
    }
    finals <- matrix(0, n, d)
    ncluster <- 0
    savecluster <- matrix(0, 0, d)
    cluster.label <- closest.label <- rep(0, n)
    if (length(h) == 1) {
        h <- rep(h, d)
    }
    for (i in subset) {
        temp.ms <- ms.rep(X, X[i, ], h, plotms = 0, thresh = 1e-08, iter)
        finals[i, ] <- temp.ms$final
        cluster.dist <- rep(0, ncluster)
        if (ncluster >= 1) {
            for (j in 1:ncluster) {
                cluster.dist[j] <- enorm(savecluster[j, ] - finals[i, 
                  ])/enorm(savecluster[j, ])
            }
        }
       
        if (ncluster == 0 || min(cluster.dist) > thr) {
            ncluster <- ncluster + 1
            savecluster <- rbind(savecluster, finals[i, ])
            cluster.label[i] <- ncluster
        } else {
            cluster.label[i] <- which(cluster.dist == min(cluster.dist))
        }
        if (d == 2 && plotms == 1) {
            lines(rbind(temp.ms$start, temp.ms$Meanshift.points), 
                col = "grey30")
        }
        if (d == 2 && plotms == 2) {
            lines(rbind(temp.ms$start, temp.ms$Meanshift.points), 
                col = cluster.label[i] + 1)
        }
    }
   # print(finals)
    for (i in subset){
         closest.label[i] <- mindist(savecluster, X[i,])$closest.item
         #closest.coords[i,]<- object$cluster.center[closest.center,]
       }
    if (d == 2 && plotms == 1) {
        points(finals, pch = 15, col = 2)
    }
    if (d == 2 && plotms == 2) {
        points(finals, pch = 15)
    }
    if (d > 2 && plotms > 1) {
        pairs(rbind(as.matrix(X), savecluster), col = c(cluster.label + 
            1, rep(1, dim(savecluster)[1])), pch = c(rep(20, 
            dim(X)[1]), rep(24, dim(savecluster)[1])), ...)
    }
    dimnames(savecluster) <- list(1:ncluster, NULL)

    fit <- list(
                cluster.center = savecluster,
                cluster.label = cluster.label,
                closest.label = closest.label,
                h=h,
                data = X,              
                scaled= scaled, 
                scaled.by =  if (scaled) s1 else rep(1,d)
                )
    class(fit) <- "ms"
    return(fit)
}


# Mean shift clustering bandwidth selection
ms.self.coverage <-
function (X, taumin = 0.02, taumax = 0.5, gridsize = 25, thr = 0.0001, 
    scaled = TRUE, cluster = FALSE,  plot.type = "o", or.labels = NULL, print=FALSE, 
    ...) 
{

    if (gridsize <10){stop("The minimum gridzise is 10.")} 
  
    X <- as.matrix(X)
    if ((!scaled) && taumax < 1) {
        warning("Please adjust the range (taumin, taumax) of tube widths by hand, as the data are not scaled.")
    }
    Pm <- NULL
    h0 <- taumin
    h1 <- taumax
    h <- seq(h0, h1, length = gridsize)
    n <-dim(X)[1]
    cover <- matrix(0, gridsize, 2)
    for (i in 1:gridsize) {
        new.h0 <- h[i]
        fit <- ms(X, new.h0,  thr = thr, scaled = scaled, plotms = 0,  or.labels=or.labels)
        #set <-   sample(n,n %/% (1/draw))
        #find <- as.numeric(names(table(fit$cluster.label[set])))
        find <- as.numeric(which(table(fit$cluster.label)>2))  # changed 23/05/11
        Pm[[i]] <- fit$cluster.center[find,] 
        if (!is.matrix(Pm[[i]])){Pm[[i]]<- matrix(Pm[[i]],nrow=1, dimnames=list(dimnames(fit$cluster.center)[[1]][find] ,NULL))}
        if (!cluster) {   
            cover[i, ] <- as.numeric(coverage.raw(fit$data, Pm[[i]],  new.h0, plot.type = 0, print=print)[1:2])
        } else {
            cover[i, ] <- as.numeric(coverage.raw(fit$data, Pm[[i]], new.h0, plot.type = 0, label = fit$cluster.label, print=print)[1:2])
        }
        
    }
    select <- select.self.coverage(self = cover, 
        smin = 1/3, plot.type = plot.type)
    result <- list(self.coverage.curve = cover, select = select, 
        type = "ms")
    class(result) <- "self"
    result
}
