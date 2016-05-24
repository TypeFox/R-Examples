dbscan <- function(
  data
, eps
, MinPts    = 5
, scale     = FALSE
, method    = c("hybrid","raw","dist")
# , no.check  = FALSE
, seeds     = TRUE
, showplot  = FALSE
, countmode = NULL #c(1,2,3,5,10,100,1000,5000,10000,50000)
)
{
#  if (!require(distpatch))
  distcomb <- function(x,data){
    data <- t(data)
    temp <- apply(x, 1, function(x){
      sqrt(colSums((data-x)^2))
    })
    if (is.null(dim(temp)))
      matrix(temp, nrow(x), ncol(data))
    else
      t(temp)
  }
  method <- match.arg(method)
  data <- as.matrix(data)
  n <- nrow(data)
  if (scale)
    data <- scale(data)
  classn <- cv <- integer(n)
  isseed <- logical(n)
  cn <- integer(1)
  for (i in 1:n){
    if (i %in% countmode)
      cat("Processing point ", i," of ",n, ".\n")
    unclass <- (1:n)[cv<1]
    if (cv[i]==0){
      if (method=="dist"){
        reachables <- unclass[data[i,unclass]<=eps]
      }else{
        reachables <- unclass[as.vector(distcomb(data[i,, drop=FALSE],data[unclass,, drop=FALSE]))<=eps]
      }
      if (length(reachables)+classn[i]<MinPts)
        cv[i] <- (-1)
      else{
        cn <- cn+1
        cv[i] <- cn
        isseed[i] <- TRUE
        reachables <- setdiff(reachables, i)
        unclass <- setdiff(unclass, i)
        classn[reachables] <- classn[reachables]+1
        while (length(reachables)){
          if (showplot)
            plot(data,  col=1+cv, pch=1+isseed)
          cv[reachables] <- cn
          ap <- reachables
          reachables <- integer()
          if (method=="hybrid"){
            tempdist <- distcomb(data[ap, , drop=FALSE], data[unclass, , drop=FALSE])
            frozen.unclass <- unclass
          }
          for (i2 in seq(along=ap)){
            j <- ap[i2]
            if (showplot>1)
              plot(data, col=1+cv, pch=1+isseed)
            if (method=="dist"){
              jreachables <- unclass[data[j,unclass]<=eps]
            }else if (method=="hybrid"){
              jreachables <- unclass[tempdist[i2,match(unclass, frozen.unclass)]<=eps]
            }else{
              jreachables <- unclass[as.vector(distcomb(data[j,, drop=FALSE], data[unclass,, drop=FALSE]))<=eps]
            }
            if (length(jreachables)+classn[j]>=MinPts){
              isseed[j] <- TRUE
              cv[jreachables[cv[jreachables]<0]] <- cn
              reachables <- union(reachables, jreachables[cv[jreachables]==0])  # isseed for these new reachables tested at next while loop
            }
            # must be after querying classn, otherwise we count j itself twice
            classn[jreachables] <- classn[jreachables]+1
            unclass <- setdiff(unclass, j)
          } # for j
        } # while sum reachables>0
      } # else (sum reachables + ... >= MinPts)
    } # if cv==0
    if (!length(unclass))
      break
  } # for i
  rm(classn)
  if (any(cv==(-1))){
    cv[cv==(-1)] <- 0
  }
  if (showplot)
    plot(data,  col=1+cv, pch=1+isseed)
  out <- list(
    cluster  = cv
  , eps      = eps
  , MinPts   = MinPts
  )
  if (seeds && cn>0){
    out$isseed <- isseed
  }
  class(out) <- "dbscan"
  out
} # dbscan


print.dbscan <- function(x, ...){
  cat("dbscan Pts=", length(x$cluster), " MinPts=", x$MinPts, " eps=", x$eps, "\n", sep="")
  if (is.null(x$isseed))
    tab <- table(x$cluster)
  else{
    tab <- table(c("seed", "border")[2-x$isseed], cluster=x$cluster)
    if (is.null(dim(tab))){
      tab <- cbind(tab)
      colnames(tab) <- unique(x$cluster)
    }
    tab <- rbind(tab, total=colSums(tab))
  }
  print(tab, ...)
}

plot.dbscan <- function(x, data, ...)
{
  plot(data, col=1+x$cluster, pch=1+x$isseed, ...)
}


predict.dbscan <- function(
  object
, data
, newdata     = NULL
, predict.max = 1000
# , no.check    = FALSE
, ...
)
{
  if (is.null(newdata)){

    return(object$cluster)

  }else{

    if (is.null(object$isseed))
      stop("no seeds to predict")

    dmax <- object$eps
    data <- data[object$isseed, , drop=FALSE]
    out <- object$cluster[object$isseed]

#    if (!require(distpatch))
    distpair <- function(x,data){
      sqrt(rowSums((x-data)^2))
    }

#    require(class)
    batchpredict <- function(newdata){
      w <- as.integer(knn1(data, newdata, 1:n.orig))
      newout <- out[w]
      if (!is.null(dmax)){
        d <- distpair(data[w,,drop=FALSE], newdata)
        newout[d>dmax] <- 0
      }
      return(newout)
    }
    n <- nrow(newdata)
    n.orig <- nrow(data)
    if (n>predict.max){
      i <- 1:n
      ret <- do.call("c", lapply(split(i, (i-1)%/%predict.max), function(i)batchpredict(newdata[i, , drop=FALSE])))
    }else{
      ret <- batchpredict(newdata)
    }
    return(ret)
  }
}



# if (FALSE){
# 
#   x <- t(t(sort(c(rnorm(20), 1:10))))
#   ds1 <- dbscan1(x, MinPts=5, eps=2, showplot=1)
#   ds <- dbscan(x, MinPts=5, eps=2, showplot=1)
# 
#   par(mfrow=c(2, 1))
#   plot(x, col=1+ds1$classification)
#   plot(ds, x)
#   ds1
#   ds
#   par(mfrow=c(1, 1))
# 
# }
