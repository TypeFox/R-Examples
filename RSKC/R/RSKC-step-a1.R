

RSKC.trimkmeans <-function(d,ncl,trim,runs=1,points=Inf,maxit=nrow(as.matrix(d))*2)
  {
    ## points must be either Inf (i.e., no specification) or k by p matrix where each row represents
    ## the cluster center 
    n <- nrow(d)
    p <- ncol(d)
    nin <- ceiling((1-trim)*n)
    ## If the starting point is initialized, trimmed k-means run only once
    if (sum(is.finite(points)) > 0) nstart <- 1
    
    re <- .Call("RSKC_trimkmeans",
                as.numeric(c(d)), ## c(d) = c(d[,1],d[,2],...d[,ncol(d)])
                as.integer(n),
                as.integer(p),
                as.integer(ncl),
                as.integer(nin),
                as.integer(runs),
                as.numeric(c(t(points))),  ## c(t(points)) = (cluster center 1, cluster center 2, ...)
                as.integer(maxit),
                PACKAGE="RSKC"
                )
    names(re) <- c("labels","means","oW","WSS","classification")
    ## re$means[ik,ip] = re$means[ik+k*ip] (from C)
    
    re$means <- matrix(re$means,ncl,p,byrow=FALSE)
    
    if (nin == n)
      {
        re$oW <- "undefined"
      }else{
        ## C to R adjustment
        re$oW <- re$oW + 1
      }
    re$labels <- re$labels + 1
    re$classification <- re$classification + 1
    re$classification[re$classification==0] <- ncl + 1
    return (re)
  }



RSKC.trimkmeans.missing <-function(d,ncl,w=rep(1,nrow(d)),trim=0.1,runs=1,points=Inf,maxit=nrow(as.matrix(d))*2)
  {
    ## points must be either Inf (i.e., no specification) or k by p matrix where each row represents
    ## the cluster center 
    n <- nrow(d)
    p <- ncol(d)
    nin <- ceiling((1-trim)*n)

    ##cat("\n\n nin",nin)
    ## If the starting point is initialized, trimmed k-means run only once
    if (sum(is.finite(points)) > 0) nstart <- 1
    
    re <- .Call("RSKC_trimkmeans_missing",
                as.numeric(c(d)), ## c(d) = c(d[,1],d[,2],...d[,ncol(d)])
                as.integer(n),
                as.integer(p),
                as.integer(ncl),
                as.integer(nin),
                as.integer(runs),
                as.numeric(c(t(points))),  ## c(t(points)) = (cluster center 1, cluster center 2, ...)
                as.integer(maxit),
                as.numeric(w),
                PACKAGE="RSKC"
                )
    names(re) <- c("labels","means","oW","WSS","classification")
    ## re$means[ik,ip] = re$means[ik+k*ip] (from C)
    
    re$means <- matrix(re$means,ncl,p,byrow=FALSE)
    
    if (nin == n)
      {
        re$oW <- "undefined"
      }else{
        ## C to R adjustment
        re$oW <- re$oW + 1
      }
    re$labels <- re$labels + 1
    re$classification <- re$classification + 1
    re$classification[re$classification==0] <- ncl + 1

    ## cat("\n\n re$classification",re$classification,"length",length(re$classification))
    ## cat("\n\n re$oW",re$oW,"length",length(re$oW))
    
    return (re)
  }



WDISC <-
  function(D,mu,ncl,n,w,sumW)
{
  ## D, mu and w must be in reduced dimention	
  ## ncol(D)==length(w)
  dist<-matrix(NA,nrow=n,ncol=ncl);  
  for ( i in 1:ncl){
      s <- sc <- scale(D, center=mu[i,], scale=FALSE)^2
      s[ is.na(sc) ] <- 0;vec<-rowSums(s)
      
      adjust<-sumW/((!is.na(sc))%*%w) # n by 1 adjusted sclaers 
      dist[,i]<-vec*adjust 
    }      
  return(dist)
}



temp <- function(vec1,vec2)
{
  .Call("tempC",as.numeric(vec1),as.numeric(vec2),PACKAGE="RSKC")
}


## ====== Old R codes ======

## RSKC.trimkmeans <- function(data,k,trim=0.1, scaling=FALSE,
##                             runs=100, points=NULL,
##                             countmode=runs+1,
##                             maxit=2*nrow(as.matrix(data))){
##   ## the code of trimkmeans is modified so that
##   ## for alpha = 0, no case is trimmed. i.e. k-means clustering is performed                     	
##   data <- as.matrix(data)
##   n <- nrow(data); nc<-ncol(data)
##   nin <- ceiling((1-trim)*n)
##   if (scaling) data <- scale(data)
##   crit <- Inf
##   oldclass <- iclass <- optclass <- rep(0,n)
##   disttom <- rep(0,n)
##                                         #  optmeans <- data[sample(n,k),,drop=FALSE]
##   for (i in 1:runs)
##     {
##       if ((i/countmode)==round(i/countmode)) cat("Iteration ",i,"\n")
##       if (is.null(points))
##         means <- data[sample(n,k),,drop=FALSE]
##       else
##         means <- points
##       wend <- FALSE
##       itcounter <- 0
##       while(!wend){
##         itcounter <- itcounter+1
##         reF<-.Fortran(
##                       "disttom_iclass",
##                       as.double(data),as.integer(n),as.integer(nc),
##                       as.double(means),as.integer(k),
##                       iclass=as.integer(iclass),
##                       disttom=as.double(disttom),
##                       PACKAGE="RSKC"
##                       )            
##         disttom<-reF$disttom
##         iclass<-reF$iclass
        
##                                         # define outliers **it could be that some clusters do not have obs**
##         iclass0<-iclass;iclass0[order(disttom)[(nin+1):n]] <- 0
##         if (trim!=0) iclass<-iclass0
##                                         #    newcrit <- sum(disttom[iclass>0])
##                                         #    cat("Iteration ",i," criterion value ",newcrit,"\n")

##                                         # stopping criteria if the class labels are the same then stop the iteration..
##         if (itcounter>=maxit | identical(oldclass,iclass)) wend <- TRUE
##         else{
##           for (l in 1:k){
##                                         # if a cluster is empty then cluster center is the first outlier obs ??
##             if (sum(iclass==l)==0) means[l,] <- data[iclass0==0,,drop=FALSE][1,]
##             else{
##               if (sum(iclass==l)>1){
##                 if (dim(means)[2]==1)
##                   means[l,] <- mean(data[iclass==l,])
##                 else
##                   means[l,] <- colMeans(data[iclass==l,])
##               }
##               else means[l,] <- data[iclass==l,]
##             }
##           }
##           oldclass <- iclass
##         }
##       }
##       newcrit <- sum(disttom[iclass>0])
##       if (newcrit<=crit){
##         optclass <- iclass
##         crit <- newcrit
##         optmeans <- means
##       }
##     }
##   optclass[optclass==0] <- k+1
##   out <- list(classification=optclass,means=optmeans,
##               criterion=crit/nin,disttom=disttom,ropt=sort(disttom)[nin],
##               k=k,trim=trim,runs=runs,scaling=scaling)
##                                         #class(out) <- "tkm"
##   out
## }



## RSKC.trimkmeans.missing <- function(data,k,w,
##                                     sumW=sum(w),trim=0.1,scaling=FALSE, runs=100, points=NULL,
##                                     countmode=runs+1, printcrit=FALSE,
##                                     maxit=2*nrow(as.matrix(data))){
##   ## w must be; length(w)==ncol(data)
##   data <- as.matrix(data);
##   n <- nrow(data);p <-ncol(data)
##   nin <- ceiling((1-trim)*n)
##   if (scaling) data <- scale(data)
##   crit <- Inf
##   oldclass <- iclass <- optclass <- rep(0,n)
##   disttom <- rep(0,n)
  
##   ##  optmeans <- data[sample(n,k),,drop=FALSE]
##   for (i in 1:runs){
##     if ((i/countmode)==round(i/countmode)) cat("Iteration ",i,"\n")
##     if (is.null(points))
##       means <- data[sample(n,k),,drop=FALSE]
##     else
##       means <- points
##     wend <- FALSE
##     itcounter <- 0
##     while(!wend){
##       itcounter <- itcounter+1
##       reF <- wrapper.fort(data,nrow=n,ncol=p,mu=means,k=k,W=w,sumW=sumW)
##       disttom <- reF$disttom
##       iclass <- reF$iclass

##       iclass0<-iclass;iclass0[order(disttom)[(nin+1):n]] <- 0
##       if (trim!=0) iclass<-iclass0
      
##       ## stopping criteria if the class labels are the same
##       ## then stop the iteration..
##       if (itcounter>=maxit | identical(oldclass,iclass)) wend <- TRUE
##       else{
##         for (l in 1:k){
##           ## if a cluster is empty then cluster center is the farthest case from 
##           ## its cluster center
##           ## if this is true, then cluster center can contains missing value..
##           if (sum(iclass==l)==0) means[l,] <- data[iclass0==0,,drop=FALSE][1,]
##           else{
##             if (sum(iclass==l)>1){
##               if (dim(means)[2]==1) # if only one feature in data
##                 means[l,] <- mean(data[iclass==l,],na.rm=TRUE)
##               else
##                 means[l,] <- colMeans(data[iclass==l,],na.rm=TRUE)
##             }
##             ## if the l^th cluster cotains only one element,
##             ## cluster center can contains missing values..
##             else means[l,] <- data[iclass==l,] 
##           }
##         }
##         oldclass <- iclass
##       }
##     }
##     newcrit <- sum(disttom[iclass>0])
##     if (printcrit) cat("Iteration ",i," criterion value ",
##                        newcrit/nin,"\n")
##     if (newcrit<=crit){
##       optclass <- iclass
##       crit <- newcrit
##       optmeans <- means
##     }
##   }
##   optclass[optclass==0] <- k+1
##   out <- list(classification=optclass,means=optmeans,
##               criterion=crit/nin,disttom=disttom,ropt=sort(disttom)[nin],
##               k=k,trim=trim,runs=runs,scaling=scaling)
##   ##class(out) <- "tkm"
##   out
## }



## wrapper.fort <- function(data,nrow,ncol,mu,k,W,sumW)
##   {
##     miss_data <- matrix(1,nrow,ncol)
##     miss_data[is.na(data)] <- 0
##     data[miss_data==0] <- 0
    
##     miss_mu <- matrix(1,k,ncol)
##     miss_mu[is.na(mu)] <- 0
##     mu[miss_mu==0] <- 0
##                                         # outputs
##     iclass <- rep(0,nrow)
##     disttom <- rep(0,nrow)
    
##     result<-.Fortran("disttom_iclass_missing",
##                      as.double(data),
##                      as.integer(miss_data),
##                      as.integer(nrow),as.integer(ncol),
##                      as.double(mu),
##                      as.integer(miss_mu),
##                      as.integer(k),
##                      iclass=as.integer(iclass),
##                      disttom=as.double(disttom),
##                      as.double(W),
##                      as.double(sumW),
##                      PACKAGE="RSKC")
##     if (sum(is.na(result$disttom))!=0)
##       stop("L1 is too small (sparse) for a dataset with missing values!")
##     return(list(disttom=result$disttom,iclass=result$iclass))
##   }
