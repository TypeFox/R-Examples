## This file (R/truthTable.R) is part of QCA3 package
## copyright: HUANG Ronggui 2008-2012

consistency <- function(x,...){
    UseMethod("consistency")
}

consistency.default <- function(x,y,alternative=c("less","greater"), ...){
  ## consistency(x<=y) when alternative="less"
  ## consistency(x>=y) when alternative="greater"
  ## x and y must in [0,1]
  allvalues <- !is.na(x) & !is.na(y)
  x<-x[allvalues]
  y<-y[allvalues]
  alternative <- match.arg(alternative)
  Sum <- switch(
                alternative,
                less=sum(x),
                greater=sum(y),
                )
  Min <- pmin(x,y)
  ans <- sum(Min)/Sum
  return(ans)
}

overlap <- function(x,y){
  ## x and y must in [0,1]
  allvalues <- !is.na(x) & !is.na(y)
  x<-x[allvalues]
  y<-y[allvalues]
  Min <- pmin(x,y)
  ans <- sum(Min)
  return(ans)
}

coverage <- function(x,...){
    UseMethod("coverage")
}

coverage.default <- function(x,y,alternative=c("less","greater"),...){
  ## coverage(x<=y) when alternative="less"
  ## coverage(x>=y) when alternative="greater"
  ## x and y must in [0,1]
  allvalues <- !is.na(x) & !is.na(y)
  x<-x[allvalues]
  y<-y[allvalues]
  alternative <- match.arg(alternative)
  Sum <- switch(
                alternative,
                less=sum(y),
                greater=sum(x),
                )
  Min <- pmin(x,y)
  ans <- sum(Min)/Sum
  return(ans)
}


## coverage.QCA <- function(x,traditional=TRUE,...){
##     explain <- x$call$explain
##     if (is.null(explain)) explain <- "positive" ## will be null if default
##     truthTable <- x$truthTable
##     Cases <- truthTable[rownames(truthTable) %in% rownames(x$explained), "Cases"]
##     OUT <- truthTable[rownames(truthTable) %in% rownames(x$explained), "OUT"]
##     ## make use of rownames here
##     if (pmatch(explain,"positive",0)==1) N_explained <- truthTable[rownames(truthTable) %in% rownames(x$explained), "freq1"]
##     if (pmatch(explain,"negative",0)==1) N_explained <- truthTable[rownames(truthTable) %in% rownames(x$explained), "freq0"]
##     N_total <- sum(truthTable["NCase"])
##     N_positive <- sum(truthTable["freq1"])
##     N_negative <- sum(truthTable["freq0"])
##     N_T_explained <- sum(N_explained)
##     N_explainedMatrix <- matrix(N_explained,
##                                 nrow=nrow(x$explained),
##                                 ncol=nrow(x$solutionsIDX)
##                                 )
##     idxMatrix <- lapply(x$solutions,
##            function(eachSolution) {
##                apply(eachSolution,1,
##                      function(im) {
##                          rownames(x$explained) %in% QCA3:::subSet(im)
##                      })
##            }
##            )## end of idxMatrix. It can constructed through PIChart and solutionIDX too.
##     ## a list of logic matrix (ncol=nrow(object$explained), nrow=number of implicants in a solution)
##     ## each column represents one implicant
##     ## TRUE suggest the implicant covers a
##     implicantString <- lapply(x$solutions,
##                               function(eachSolution){
##                                 apply(eachSolution,1,toString,traditional=traditional,
##                                       nlevels=x$nlevels,name=names(x$solutions[[1]]))
##                               }) ## end of implicantString
##     coverage <- vector("list",length(x$solutions))
##     for (i in 1:length(x$solutions)){
##         NMatrix <- N_explainedMatrix
##         NMatrix[which(!idxMatrix[[i]])] <- 0
##         raw <- apply(NMatrix,2,sum)
##         unique <- c()
##         for (j in 1:ncol(N_explainedMatrix)) {
##             idxUnique <- apply(NMatrix[,-j,drop=FALSE],1,FUN=function(eachRow) all(eachRow==0))
##             if (any(idxUnique)) unique <- c(unique,sum(NMatrix[,j,drop=TRUE][idxUnique])) else
##             unique <- c(unique,0)
##         }
##         coverage[[i]] <- data.frame(implicant=implicantString[[i]],
##                                     raw=raw,unique=unique)
##       }
##     cat(sprintf("Total number of cases: %i\n",N_total))
##     cat(sprintf("Number of cases [1]: %i\n",N_positive))
##     cat(sprintf("Number of cases [0]: %i\n",N_negative))
##     cat(sprintf("Number of cases to explain: %i\n\n", N_T_explained))
##     cat("Coverage of each solution.\n\n")
##     coverage
##   }

suffnec <- function(x, use=c("complete","pairwise")){
  consistency_fn <- function(x,y,alternative=c("xley","ylex"),...){
    ## helper function
    allvalues <- !is.na(x) & !is.na(y)
    x<-x[allvalues]
    y<-y[allvalues]
    Min <- pmin(x,y)
    alternative <- match.arg(alternative)
    Sum <- switch(
                  alternative,
                  xley=sum(x),
                  ylex=sum(y),
                  )
    ans <- sum(Min)/Sum
    return(ans)
  }

  minmax <- range(x,na.rm=TRUE)
  if (minmax[1] < 0 || minmax[2] > 1) stop("All the values of 'x' should be range from 0 to 1.")
  use <- match.arg(use)
  if (use=="complete") x= na.exclude(x)
  Nvar <- ncol(x)
  ans <- matrix(numeric(0),nrow=Nvar,ncol=Nvar)
  index <- t(combn(Nvar,2)) # the first col is the column index
  nindex <- nrow(index)

  for (i in 1:nindex) {
    rindex <- index[i,][2]
    cindex <- index[i,][1] ## row change fast -> fill lower matrix first.
    ans[rindex,cindex] <- consistency_fn(x=x[,rindex],y=x[,cindex],"xley")
    ## sufficient condiction
    ans[cindex,rindex] <- consistency_fn(x=x[,rindex],y=x[,cindex],"ylex")
    ## necessary condiction
  }
  dimnames(ans)= list(X=names(x),Y=names(x))
  diag(ans) <- 1
  ans2 <- t(ans)
  dimnames(ans2)= list(X=names(x),Y=names(x))
  result <- list(suff=ans,nec=ans2)
  class(result) <- "suffnec"
  result
}

print.suffnec <- function(x,digits=3,...)
{
  x<-unclass(x)
  cat("\nNecessity Scores Matrix:\n'X is necessary condition of Y'\n")
  print(x$nec,digits=digits,na.print=" ",quote = FALSE,...)
  cat("\nSufficiency Scores Matrix:\n'X is sufficient condition of Y'\n")
  print(x$suff,digits=digits,na.print=" ",quote = FALSE,...)
}


tnec <- function(x, y){
    ## x must greater or equal to y?
    ## how to handle y=1?
    ans <- mean((1 - x) / (1 - y))
    ans
}
