trimkmeans <- function(data,k,trim=0.1, scaling=FALSE, runs=100, points=NULL,
                       countmode=runs+1, printcrit=FALSE,
                       maxit=2*nrow(as.matrix(data))){
  data <- as.matrix(data)
  n <- nrow(data)
  nin <- ceiling((1-trim)*n)
  if (scaling) data <- scale(data)
  crit <- Inf
  oldclass <- iclass <- optclass <- rep(0,n)
  disttom <- rep(0,n)
#  optmeans <- data[sample(n,k),,drop=FALSE]
  for (i in 1:runs){
    if ((i/countmode)==round(i/countmode)) cat("Iteration ",i,"\n")
    if (is.null(points))
      means <- data[sample(n,k),,drop=FALSE]
    else
      means <- points
    wend <- FALSE
    itcounter <- 0
    while(!wend){
      itcounter <- itcounter+1
      for (j in 1:n){
        dj <- rep(0,k)
        for (l in 1:k)
          dj[l] <- sum((data[j,]-means[l,])^2)
        iclass[j] <- which.min(dj)
        disttom[j] <- min(dj)
      }      
      iclass[order(disttom)[(nin+1):n]] <- 0
#    newcrit <- sum(disttom[iclass>0])
#    cat("Iteration ",i," criterion value ",newcrit,"\n")
      if (itcounter>=maxit | identical(oldclass,iclass)) wend <- TRUE
      else{
        for (l in 1:k){
          if (sum(iclass==l)==0) means[l,] <- data[iclass==0,,drop=FALSE][1,]
          else{
            if (sum(iclass==l)>1){
              if (dim(means)[2]==1)
                means[l,] <- mean(data[iclass==l,])
              else
                means[l,] <- colMeans(data[iclass==l,])
            }
            else means[l,] <- data[iclass==l,]
          }
        }
        oldclass <- iclass
      }
    }
    newcrit <- sum(disttom[iclass>0])
    if (printcrit) cat("Iteration ",i," criterion value ",
                       newcrit/nin,"\n")
    if (newcrit<=crit){
      optclass <- iclass
      crit <- newcrit
      optmeans <- means
    }
  }
  optclass[optclass==0] <- k+1
  out <- list(classification=optclass,means=optmeans,
              criterion=crit/nin,disttom=disttom,ropt=sort(disttom)[nin],
              k=k,trim=trim,runs=runs,scaling=scaling)
  class(out) <- "tkm"
  out
}

print.tkm <- function(x,...){
  cat("* trimmed k-means *\n")
  cat("trim=",x$trim,", k=",x$k,"\n")
  cat("Classification (trimmed points are indicated by ",x$k+1,"):\n")
  print(x$classification)
  cat("Means:\n")
  print(x$means)
  cat("Trimmed mean squares: ",x$criterion,"\n")
  invisible(x)
}

plot.tkm <- function(x,data,...){
  p <- dim(as.matrix(data))[2]
  if (p==1){
    require(fpc)
    plotcluster(data,x$classification, ...)
  }
  if (p==2){
# Create the graphs for summarizing results
    plot(data,type="n",main=paste("k = ", x$k, " and trim = ",x$trim ))
    points(data[x$classification==x$k+1,],col=1,pch=3)
    t <- seq(0,2*pi,length=1000)
    j <- 0
    for (i in 1:x$k){
      j <- j+1
      points(data[x$classification==j,],col=j+1)
      circ <- sqrt(x$ropt)*cbind(cos(t),sin(t));
      lines(rep(1,1000)%*%t(x$means[j,])+circ,col=1)
    }
  }
  if (p>2){
    require(fpc)
    cv <- x$classification<x$k+1
    dcx <- discrcoord(data[cv,],x$classification[cv])
    dproj <- data %*% dcx$units
    pchs <- c(sapply(1:x$k,toString),"T")
    plot(dproj,col=x$classification,pch=pchs[x$classification],
         main="Discriminant coordinates of clusters (T=trimmed)", ...)
  }
  invisible(x)
}
 
  
  
