## this is part of the QCA3 project
## by Ronggui Huang 2009-2010

excludeCSA <- function(object,csa){
  call <- match.call()
  nlevels <- object$nlevels
  conditions <- names(object$explained)
  superSets1 <- apply(object$explained, 1, superSet,nlevels=nlevels)
  dim(superSets1) <- NULL
  superSets1 <- unique(superSets1)
  superSets0 <-  apply(id2Implicant(object$idExclude,nlevels),
                       1,superSet,nlevels=nlevels)
  dim(superSets0) <- NULL
  superSets0 <- unique(superSets0)
  superSetsCSA <- apply(csa$solutions[[1]],1, superSet,nlevels=nlevels)
  dim(superSetsCSA) <- NULL
  superSetsCSA <- unique(superSetsCSA)
  primesId <- sort(setdiff(superSets1, unique(c(superSets0,superSetsCSA))))
  primesId <- ereduce1(primesId, nlevels = nlevels)
  primeImplicants <- id2Implicant(primesId, nlevels = nlevels, names = conditions)
  PIChart <- PIChart(primeImplicants, object$explained)
  sl <- solvePIChart(PIChart)
  solutions <- apply(sl, 2, function(idx) primeImplicants[idx,])
  commonSolutions <- apply(sl, 1, function(idx) {
    if (length(id <- unique(idx)) == 1)
      id
  })
  ans <- list(solutions = solutions, commonSolutions = commonSolutions,
              solutionsIDX = sl, primeImplicants = primeImplicants,
              truthTable = object$truthTable, explained = object$explained,
              idExclude = object$idExclude,
              nlevels = nlevels, PIChart = PIChart,call=call)
  class(ans) <- c("QCA","noCSA")
  ans
}

primeImplicants <- function(object,traditional=TRUE){
    ## extract the prime implicants and print it in a pretty way
    nlevels <- object$nlevels
    var_names <- names(object$primeImplicants)
    PIs <- apply(object$primeImplicants, 1, toString, traditional = traditional,
                 nlevels = nlevels, name = var_names)
    PI <- paste(PIs, collapse = " + ")
    writeLines(strwrap(PI))
}

consistency.QCA <- function(x, data, which=1, ...){
    ## x is a fsQCA solution, data is the original data,outcome is the outcome of QCA
    if (max(x$nlevels)>2) stop("It is not applicable for mvQCA.")
    if (which>length(x$solutions)) stop("'which' is too large.")
    sol <- x$solutions[[which]]
    outcome <- x$outcome
    ## only conduct for one solution indicated by which.
    idx1 <- which(sol==1,arr.ind=TRUE)
    idx0 <- which(sol==0,arr.ind=TRUE)
    Nimplicant <- nrow(sol)
    conds <- names(sol)
    ans <- data.frame(consistency=numeric(Nimplicant+1))
    solX <- matrix(numeric(nrow(data)*Nimplicant),ncol=Nimplicant)
    for (i in seq(Nimplicant)) {
        cond1 <- conds[idx1[idx1[,1]==i,2]]
        dat1 <- data[,cond1,drop=FALSE]
        cond0 <- conds[idx0[idx0[,1]==i,2]]
        if (length(cond0)>0) dat0 <- 1-data[,cond0,drop=FALSE] else dat0 <- data[,cond0,drop=FALSE] ## empty data frame
        if (ncol(dat1)>0 & ncol(dat0)>0) {
            soli <- cbind(dat1,dat0)
        } else if(ncol(dat1)==0) {
            soli <- dat0
        } else {
            soli <- dat1
        }
        fzx <- apply(soli,1,min)
        solX[,i] <- fzx
        ans[i,"consistency"] <- consistency(x=fzx,y=data[,outcome])
    }
    ans[i+1,"consistency"] <- consistency(x=apply(solX,1,max),y=data[,outcome])
    implicantName <- apply(sol,1,function(obj) toString(obj,traditional=TRUE,nlevels=x$nlevels,conds))
    rownames(ans) <- c(implicantName,"[solution]")
    ans
}

coverage.QCA <- function(x, data, type=c("raw","unique"), which=1, ...){
    if (max(x$nlevels)>2) stop("It is not applicable for mvQCA.")
    type <- match.arg(type)
    ans <- switch(type,
                  raw= rawCoverageQCA(x,data,which),
                  unique= uniqueCoverageQCA(x,data,which)
                  )
    ans
}

rawCoverageQCA <- function(x, data, which=1){
    ## x is a fsQCA solution, data is the original data,outcome is the outcome of QCA
    if (max(x$nlevels)>2) stop("It is not applicable for mvQCA.")
    if (which>length(x$solutions)) stop("Which is too large.")
    sol <- x$solutions[[which]]
    outcome <- x$outcome
    ## only conduct for one solution indicated by which.
    idx1 <- which(sol==1,arr.ind=TRUE)
    idx0 <- which(sol==0,arr.ind=TRUE)
    Nimplicant <- nrow(sol)
    conds <- names(sol)
    ans <- data.frame(rawCoverage=numeric(Nimplicant+1))
    solX <- matrix(numeric(nrow(data)*Nimplicant),ncol=Nimplicant)
    for (i in seq(Nimplicant)) {
        cond1 <- conds[idx1[idx1[,1]==i,2]]
        dat1 <- data[,cond1,drop=FALSE]
        cond0 <- conds[idx0[idx0[,1]==i,2]]
        if (length(cond0)>0) dat0 <- 1-data[,cond0,drop=FALSE] else dat0 <- data[,cond0,drop=FALSE] ## empty data frame
        if (ncol(dat1)>0 & ncol(dat0)>0) {
            soli <- cbind(dat1,dat0)
        } else if(ncol(dat1)==0) {
            soli <- dat0
        } else {
            soli <- dat1
        }
        fzx <- apply(soli,1,min)
        solX[,i] <- fzx
        ans[i,"rawCoverage"] <- coverage(x=fzx,y=data[,outcome])
    }
    ans[i+1,"rawCoverage"] <- coverage(x=apply(solX,1,max),y=data[,outcome])
    implicantName <- apply(sol,1,function(obj) toString(obj,traditional=TRUE,nlevels=x$nlevels,conds))
    rownames(ans) <- c(implicantName,"[solution]")
    ans
}

uniqueCoverageQCA <- function(x, data, which=1){
    ## x is a fsQCA solution, data is the original data,outcome is the outcome of QCA
    if (max(x$nlevels)>2) stop("It is not applicable for mvQCA.")
    if (which>length(x$solutions)) stop("Which is too large.")
    sol <- x$solutions[[which]]
    ## only conduct for one solution indicated by which.
    outcome <- x$outcome
    idx1 <- which(sol==1,arr.ind=TRUE)
    idx0 <- which(sol==0,arr.ind=TRUE)
    Nimplicant <- nrow(sol)
    conds <- names(sol)
    ans <- data.frame(uniqueCoverage=numeric(Nimplicant+1))
    solX <- matrix(numeric(nrow(data)*Nimplicant),ncol=Nimplicant)

    for (i in seq(Nimplicant)) {
        cond1 <- conds[idx1[idx1[,1]==i,2]]
        dat1 <- data[,cond1,drop=FALSE]
        cond0 <- conds[idx0[idx0[,1]==i,2]]
        if (length(cond0)>0) dat0 <- 1-data[,cond0,drop=FALSE] else dat0 <- data[,cond0,drop=FALSE] ## empty data frame
        if (ncol(dat1)>0 & ncol(dat0)>0) {
            soli <- cbind(dat1,dat0)
        } else if(ncol(dat1)==0) {
            soli <- dat0
        } else {
            soli <- dat1
        }
        fzx <- apply(soli,1,min)
        solX[,i] <- fzx
    }
    ans[Nimplicant+1,"uniqueCoverage"] <- coverage(x=apply(solX,1,max),y=data[,outcome])
    if (Nimplicant==1){
      ## only with one recipe
      ans[1,"uniqueCoverage"] <- coverage(x=solX[,1], y=data[,outcome])
    } else {
      ## with mutiple recipes
      for (i in seq(Nimplicant)){
        notifz <- apply(solX[, -i, drop=FALSE],1, max)
        ans[i,"uniqueCoverage"] <- ans[Nimplicant+1,"uniqueCoverage"] - coverage(x=notifz, y=data[,outcome])
      }
    }
    implicantName <- apply(sol,1,function(obj) toString(obj,traditional=TRUE,nlevels=x$nlevels,conds))
    rownames(ans) <- c(implicantName,"[solution]")
    ans
}

