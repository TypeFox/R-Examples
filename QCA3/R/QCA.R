## This file (R/QCA.R) is part of QCA3 package
## copyright: HUANG Ronggui 2008-2010

## Reference
## Dusa. 2007. Enhancing Quine-McCluskey. COMPASSS Working Paper.
##   When remainders are included, this method is better. Otherwise, classic QM method will be used

## id2Implicant_old(), subSet_old() and complement1() can be found in rev 23.

#########################
### notations
### -9: indicator of dontcare case, it is assigned in a constant dontcareValue
### -1: some functions used -1 to indicate don't care case, for the sake of convenience
### NA: missing data
#########################

allGroup <- function(nlevels,names=NULL){
    ## Ragin(2000:127), see the calculation of groupings (not combinations)
    ## The total number of groupings are prod(nlevels+1); for csQCA, that is 3^k
    if (is.null(names)) names <- paste("var.",seq_len(length(nlevels)),sep="")
    ## exp <- sprintf("c(NA,0:%i)",nlevels-1)
    exp <- sprintf("c(%i,0:%i)",dontcareValue,nlevels-1)
    ans <- eval(parse(text = sprintf("expand.grid(%s)",paste(names,"=",exp,sep="",collapse=","))))
    ans
}

allCombination <- function(nlevels,names=NULL){
    ## Ragin(2000:127), see the calculation of combinations.
    ## The total number of groupings are prod(nlevels); for csQCA, that is 2^k
    if (is.null(names)) names <- paste("var.",seq_len(length(nlevels)),sep="")
    exp <- sprintf("c(0:%i)",nlevels-1)
    ans <- eval(parse(text = sprintf("expand.grid(%s)",paste(names,"=",exp,sep="",collapse=","))))
    ans
}

dontcareValue <-  -9 ## assign a constant

is.dontcare <- function(x){
    ## for a numeric vector or data frame x, elements with value of value (default=-9) is don't care case
    ## return value is logical index, TRUE for don't care case
    if (any(is.na(x))) warning("x contains missing data.") ## think about how to handle missing data
    idx <- (x == dontcareValue)
    idx
}

implicant2Id <- function(implicant,nlevels){
    ## convert an implicant to a line number, which is discussed in Dusa (2007).
    ## implicant[is.na(implicant)] <- -1
    implicant[is.dontcare(implicant)] <- -1
    ans <- sum((implicant+1) * c(1,cumprod(nlevels[-length(nlevels)]+1)) ) + 1
    ans
}

id2Implicant <- function(id,nlevels,names=NULL,to.data.frame=TRUE){
    if (is.null(names)) names <- paste("var.",seq_len(length(nlevels)),sep="")
    idx <- cumprod(nlevels+1)/(nlevels+1)
    nid <- id -1
    ans <- rep(nid,each=length(nlevels)) %/% idx %% (nlevels+1) -1
    ##  ans[ans==-1] <- NA
    ans[ans==-1] <- dontcareValue
    ans <- matrix(ans,byrow=TRUE,ncol=length(nlevels))
    colnames(ans) <- names
    rownames(ans) <- as.character(id)
    if (to.data.frame) ans <- as.data.frame(ans)
    ans
}

superSet <- function(implicant, include.itself=TRUE,rowId=TRUE,nlevels=rep(2,length(implicant))){
    ## superSet(c(1,0,1,1),nlevel=rep(2,4))
    Nvar <- length(implicant)
    index <- eval(parse(text = (sprintf("expand.grid(%s)",
                                        paste(rep("0:1",sum(!is.na(implicant))), sep = "", collapse = ",")))))
    if (include.itself) index <- index[-1,] else index <- index[-c(1,2^Nvar),]
    ans <- matrix(rep(unlist(implicant),nrow(index)),byrow=TRUE,ncol=Nvar)
    ## ans[index==0] <- NA
    ans[index==0] <- dontcareValue
    if (!is.null(colnames(implicant))) colnames(ans) <- colnames(implicant)
    if (!rowId) {
        ans
    } else {
        ans <- apply(ans,1,implicant2Id,nlevels=nlevels)
        ans
    }
}

subSet <- function(implicant,include.itself=TRUE,nlevels=rep(2,length(implicant))){
    ## new version of subSet()
    ## example: subSet(c(1,0,1,NA),nlevel=rep(2,4))
    ## idx  <- which(is.na(implicant))
    idx  <- which(is.dontcare(implicant))
    IDX <- cumprod(nlevels+1)/(nlevels+1)
    id <-  implicant2Id(implicant,nlevels=nlevels)
    if (length(idx)>0){
        nn <- nlevels[idx]
        IDX <- IDX[idx]
        exp <- sprintf("c(0:%i)",nlevels[idx])
        dat<-eval(parse(text = sprintf("expand.grid(%s)",paste(exp,sep="",collapse=","))))
        ans <- apply(dat,1,function(each) sum(IDX*each)) + id
    } else ans <- id
    if (!include.itself) ans <- ans[-1]
    ans
}


esubSet <- function(implicant,include.itself=TRUE,nlevels=rep(2,length(implicant))){
    ##enhanced version of subSet(), by using math regularity between ids.
    ## Speed improved by a factor of 2 at least.
    id <-  implicant2Id(implicant,nlevels=nlevels)
    ##idx1 <- which(is.na(implicant)) ## index of NA
    idx1 <- which(is.dontcare(implicant)) ## index of don't care case
    N <- prod(nlevels[idx1]+1)
    ans <- vector(mode = "numeric", length = N-1)
    ## if there is no dontcare case, skip to the result.
    ## use "numeric" to keep it consistent with subSet
    if ((N-1)>0) {
        idx2 <- idx1-1
        if (idx2[1]==0) {
            incr1 <- c(1,cumprod(nlevels+1)[idx2[-1]])
        } else  incr1 <- cumprod(nlevels+1)[idx2]
        incr2 <- incr1*nlevels[idx1]
        incr2 <- c(0,incr2[1:length(incr2)-1])
        incr <- incr1 - cumsum(incr2)
        idx3 <-  cumprod(nlevels[idx1]+1)
        ans[1:(idx3[1]-1)] <- incr[1]
        if (length(idx1) > 1) {
            for (i in 2:length(idx1)){
                ans[idx3[i-1]:(idx3[i]-1)] <- c(incr[i],ans[seq_len(idx3[i-1]-1)])
            }
        }
        ans <- id + cumsum(ans)
    }
    if (include.itself) ans <- c(id, ans)
    ans
}


complement1Id <- function(id,nlevels){
  IDX <- cumprod(nlevels+1)/(nlevels+1)
  idx <- (id -1) %/% IDX %% (nlevels+1)
  NApos <- which(idx==0)
  pos <- which(idx!=0)
  ans <- sapply(pos,function(sidx){
    nidx <- idx[sidx] - 1
    universal <-  sequence(nlevels[sidx])-1
    nn <- (universal + 1)  - which(universal %in% nidx)
    nn <- nn[nn!=0]
    res <- nn*IDX[sidx] + id
  })
  ans <- unlist(ans)
  attr(ans,"NApos") <- NApos
  ans
}

reduce1 <- function(IDs,nlevels){
  ## Dusa(2007: part 4) about "finding the most parsimonious prime implicants"
  IDs <- sort(IDs)
  if (length(IDs) >1){
    stop <- FALSE
    i <- 1
    while(!stop){
      IDs <- setdiff(IDs,subSet(id2Implicant(IDs[i],nlevels=nlevels),include.itself=FALSE,nlevels=nlevels))
      if (length(IDs)==i) stop <- TRUE ## should be in frot of i <- i +1
      i <- i +1
    }
    IDs
  }
}

ereduce1 <- function(IDs,nlevels){
  ## enchaned version of reduce1, use esubSet rather than subSet
  ## Dusa(2007: part 4) about "finding the most parsimonious prime implicants"
  IDs <- sort(IDs)
  if (length(IDs) >1){
    stop <- FALSE
    i <- 1
    while(!stop){
      IDs <- setdiff(IDs,esubSet(id2Implicant(IDs[i],nlevels=nlevels),include.itself=FALSE,nlevels=nlevels))
      if (length(IDs)==i) stop <- TRUE ## should be in frot of i <- i +1
      i <- i +1
    }
    IDs
  }
}

reduceByOne <- function(IDs,nlevels){
##reduceByOne(c(42,99),c(3,2,2,2))
  can_reduce <- rep(FALSE,length(IDs))
  ans_IDs <- vector("list",length(IDs))
  Ngroup <- nlevels - 1
  k <- sequence(Ngroup)
  cSum <- cumsum(Ngroup)
  for (i in seq_len(length(IDs))){
    idx1 <- complement1Id(IDs[i],nlevels=nlevels)
    NApos <- attr(idx1,"NApos")
    ## exists <- rep(FALSE,length(nlevels))
    exists <- rep(FALSE,sum(nlevels-1))
    if (length(NApos)==0)  {
      exists <- idx1 %in% IDs
    } else {
      nnn <- cumsum(nlevels-1) ## NApos may be length >1
      nnn_to = nnn[NApos]; nnn_length=nlevels[NApos] -1
      nnn <- unlist(sapply(seq_along(nnn_to),function(x) seq(to=nnn_to[x],length=nnn_length[x])))
      exists[ -nnn ] <- idx1 %in% IDs ## idx1 will be shortened
    }
    group_exists <- rep(FALSE,length(nlevels))
    jrange <- seq_along(nlevels) [!seq_along(nlevels) %in% NApos]
    for (j in jrange){
      ridx <- seq(from=cSum[j]- Ngroup[j] +1,to=cSum[j])
      group_exists[j] <- all(exists[ridx])
    }
    can_reduce[i] <- any(group_exists)
    ans_IDs[[i]] <- which(group_exists)
  }
  res <- list(reducable=can_reduce,index=ans_IDs)
  res
}


reduce2 <- function(IDs,nlevels){
    ## IDs MUST be that of commbination?
    ## otherwise, use subCommbination to discover them first.?
  reduced <- function(IDs,nlevels){
    ## helper function
    index=reduceByOne(IDs,nlevels=nlevels)
    reduce_index <- which(index$reducable)
    unreducible <- IDs[which(!index$reducable)]
    ans <- sapply(reduce_index, function(each){
      IDX <- cumprod(nlevels+1)/(nlevels+1)
      ID <- IDs[each] - 1
      reducedIDs <- IDs[each] - (IDX)*(ID %/% IDX %% (nlevels+1))
      reducedIDs <- reducedIDs[index$index[[each]]]
      ##double checked this function
    }
                  )
    ans <- unique(unlist(ans))
    res <- list(newIDs=ans, unreducible=unreducible)
    res
  } ## end of reduced()

  IDs <- sort(unique(IDs))
  stop <- FALSE
  ## final <- c()
  while(!stop){
    ans2 <- reduced(IDs=IDs,nlevels=nlevels)
    newIDs <- sort(unique(c(ans2$newIDs,ans2$unreducible)))
    ## each configuration can used more than one times (Ragin and Strand 2008:435).
    if (identical(IDs,newIDs)) stop <- TRUE
    if (!stop) IDs <- newIDs
    ##if (length(ans2$unreducible)>0) final <- c(final,ans2$unreducible)
    ## IDs <- ans2$newIDs
    ##if (is.null(ans2$newIDs)) {
    ##    stop <- TRUE
    ##}
}
  ## final <- sort(unique(final))
  ## final
  newIDs
}

PIChart <- function(primeImplicants,explained=NULL){
  ## primeImplicants with attr of "explained" if explained is NULL
  ## rows are primeimplicants and columns are original combinations
  toName <- function(x){
    NROW <- nrow(x)
    ans <- c()
    nlevels <- rep(2, ncol(x))
    nm <-names(x)
    for (i in 1:NROW){
      ans <- c(ans, toString(x[i, ], TRUE, nlevels, nm))
    }
    ans
  }
  
  if (is.null(explained)){
    explained <- attr(primeImplicants,"explained")
  }
  nr <- nrow(primeImplicants)
  nc <- nrow(explained)
  ans <- matrix(logical(0),nrow=nr,ncol=nc)
  for (i in seq_len(nr)){
    for (j in seq_len(nc)){
      ## idx <- !is.na(primeImplicants[i,])
      if (FALSE) { ## comment out the old method
        idx <- !is.dontcare(primeImplicants[i,])
        ans[i,j] <- isTRUE(all.equal(primeImplicants[i,][idx],explained[j,][idx]))
      }
      ans[i,j] <- isSuperSet(primeImplicants[i,], explained[j,])
    }
  }
  rownames(ans) <- toName(primeImplicants)
  #colnames(ans) <- toName(explained)
  colnames(ans) <- rownames(explained)
  class(ans) <- "PIChart"
  ans
}

print.PIChart <- function(x, ...){
  chart <- unclass(x)
  chart[which(x==TRUE)] <- "x"
  chart[which(x==FALSE)] <- "-"
  print(chart, quote=FALSE)
}

solvePIChart <- function(chart, all.sol=TRUE){
  ## chart: Prime implicants x Primitive Expression
  if (all(dim(chart) > 1)) {
    nPrime <- nrow(chart)
    lp  <- make.lp(0, nPrime)
    apply(chart, 2, function(xt) add.constraint(lp, xt, type = ">=", rhs=1))
    set.type(lp,1:nPrime, "binary")
    set.objfn(lp,rep(1,nPrime))
    name.lp(lp, "Simplifying a PIChart")
    # http://lpsolve.r-forge.r-project.org/
    conv <- solve(lp)
    if (conv!=0) stop("optimal solution not found")
    ans <- get.variables(lp)
    if (all.sol) {
    k <- sum(ans)
    combos <- combn(nrow(chart), k)
    sol.matrix <- combos[, apply(combos, 2, function(idx) all(colSums(chart[idx,,drop = FALSE])>0)),drop=FALSE]
    } else {
      sol.matrix <- matrix(which(ans==1), ncol=1)   
    }
  }
  else {
    sol.matrix <- matrix(seq_len(nrow(chart)),ncol=1)
  }
  sol.matrix ## now always return a matrix
}

reduce <- function(x,...){
    call <- match.call()
    UseMethod('reduce')
}

reduce.default <-
reduce.truthTable <- function(x,
                              explain=c("positive","negative"),
                              remainders=c("exclude","include"),
                              contradictions=c("remainders","positive","negative"),
                              dontcare=c("remainders","positive","negative"),
                              cdontcare=c("remainders","positive","negative"),
                              keepTruthTable=TRUE, all.sol = FALSE, ...)
{
    call <- match.call()
    call[[1]] <- as.name("reduce")
    if (!"truthTable"  %in% class(x) ) stop("x is not a truthTable.")
    mydata <- x$truthTable
    conditions <- x$conditions
    nlevels <- x$nlevels
    explain <- match.arg(explain)
    remainders <- match.arg(remainders)
    contradictions <- match.arg(contradictions)
    dontcare <- match.arg(dontcare) ## dontcare in outcome
    conddontcare <- match.arg(cdontcare) ## dontcare in condition
    if (keepTruthTable) {
        truthTable <- mydata[mydata[["OUT"]]!="?",] ## subset(mydata,OUT!="?")
        ## to avoid unbined global variable of OUT, do not use subset(mydata, OUT...)
    } else {truthTable <- NULL }
    if (dontcare=="remainders") mydata <- mydata[mydata[["OUT"]]!="-9",] ## subset(mydata,OUT!="-9" )
    if (dontcare=="positive") mydata[['OUT']][mydata[['OUT']]=="-9"] <- "1"
    if (dontcare=="negative") mydata[['OUT']][mydata[['OUT']]=="-9"] <- "0"
    dat1 <- mydata[mydata[["OUT"]]=="1",conditions] ## subset(mydata,OUT=="1",conditions)
    dat0 <- mydata[mydata[["OUT"]]=="0",conditions] ## subset(mydata,OUT=="0",conditions)
    datC <- mydata[mydata[["OUT"]]=="C",conditions] ## subset(mydata,OUT=="C",conditions)
    if (contradictions=="positive") dat1 <- rbind(dat1,datC)
    if (contradictions=="negative") dat0 <- rbind(dat0,datC)
    if (explain=="positive") {
      explained <- dat1
      idExclude <- apply(dat0,1,implicant2Id,nlevels=nlevels)
    }
    if (explain=="negative") {
      explained <- dat0
      idExclude <- apply(dat1,1,implicant2Id,nlevels=nlevels)
    }
    if (nrow(explained)==0) stop("No configuration is associated with the explained outcome.")

    if (remainders=="include"){
        ## if necessary conditons -> add some remainders to dat0
        superSets1 <- apply(dat1, 1, superSet,nlevels=nlevels)
        dim(superSets1) <- NULL
        ## set dim to NULL rather than use as.vector to speed it up.
        superSets1 <- unique(superSets1)
        superSets0 <- apply(dat0, 1, superSet,nlevels=nlevels)
        dim(superSets0) <- NULL
        superSets0 <- unique(superSets0)
        if (explain=="positive") primesId <- setdiff(superSets1,superSets0)
        if (explain=="negative") primesId <- setdiff(superSets0,superSets1)
        primesId <- ereduce1(primesId,nlevels=nlevels)
    } else if (remainders=="exclude") {
        if (conddontcare=="remainders") { ## added on 30/3/2012
            ## how to handle dontcare in conditions?
            dontcare.idx <- apply(explained, 1, function(x) any(is.dontcare(x)))
            dontcare.case <- explained[dontcare.idx,]
            primesIdP <- apply(explained[!dontcare.idx,], 1, implicant2Id, nlevels=nlevels)
            ## explained of processed cases
            dontcare.IDs <- apply(dontcare.case,1,subCombination,nlevels=nlevels)
            dontcare.IDs <- unique(unlist(dontcare.IDs))
            primesId <- unique(union(primesIdP,dontcare.IDs))
        } else if (conddontcare=="positive"){
            dat1 <- explained
            dat1[is.dontcare(explained)] <- 1
            primesId <- apply(dat1, 1, implicant2Id, nlevels = nlevels)
        } else if (conddontcare=="negative"){
            dat0 <- explained
            dat0[is.dontcare(explained)] <- 0
            primesId <- apply(dat0, 1, implicant2Id, nlevels = nlevels)
        }
        primesId <- reduce2(primesId,nlevels=nlevels)
    }

    primeImplicants <- id2Implicant(primesId ,nlevels=nlevels,names=conditions)
    PIChart <- PIChart(primeImplicants,explained)
    sl <- solvePIChart(PIChart, all.sol=all.sol)
    solutions <- apply(sl,2,function(idx)primeImplicants[idx,])
    commonSolutions <- apply(sl,1,function(idx) {if (length(id <- unique(idx))==1) id })
    ans <- list(solutions=solutions,commonSolutions=commonSolutions,
                solutionsIDX=sl,primeImplicants=primeImplicants,
                truthTable=truthTable,explained=explained,outcome=x$outcome,
                idExclude=idExclude,nlevels=nlevels,
                PIChart=PIChart, call=call, data=x$data)
    class(ans) <- c("QCA")
    ans
}

reduce.formula <- function(x, data,
                           explain=c("positive","negative"),
                           remainders=c("exclude","include"),
                           contradictions=c("remainders","positive","negative"),
                           dontcare=c("remainders","positive","negative"),
                           cdontcare=c("remainders","positive","negative"),
                           preprocess=c("cs_truthTable","fs_truthTable","mv_truthTable"),
                           keepTruthTable=TRUE, all.sol = FALSE, ...
                           )
{
    ## x is a formula
    ## note that data is mandatory
    call <- match.call()
    call[[1]] <- as.name("reduce")
    if (missing(data)) stop("argument data is missing.")
    term <- terms(x)
    if (attr(term,"response")==0) {stop("formula in the lef hand side is empty.")}
    outcome <- all.vars(attr(term,"variables")[[attr(term,"response")+1]])
    if (length(outcome)!=1) {stop("only one outcome variable is allowed.")}
    conditions <- setdiff(all.vars(x),outcome)
    if (length(conditions)<=1) stop("more conditions are needed.")
    explain <- match.arg(explain)
    remainders <- match.arg(remainders)
    contradictions <- match.arg(contradictions)
    dontcare <- match.arg(dontcare)
    conddontcare <- match.arg(cdontcare) ## dontcare in condition
    preprocess <- match.arg(preprocess)
    dots <- list(...)
    x <- do.call(preprocess,c(list(mydata=data, outcome=outcome,conditions=conditions),dots))
    ans <- do.call(reduce.truthTable,list(x=x,explain=explain,remainders=remainders,
                                          contradictions=contradictions,dontcare=dontcare,cdontcare=cdontcare,
                                          keepTruthTable=keepTruthTable, all.sol = FALSE, dots))
    ans$call <- call
    ans
}


reduce.data.frame <- function(x, outcome, conditions,
                              explain=c("positive","negative"),
                              remainders=c("exclude","include"),
                              contradictions=c("remainders","positive","negative"),
                              dontcare=c("remainders","positive","negative"),
                              cdontcare=c("remainders","positive","negative"),
                              preprocess=c("cs_truthTable","fs_truthTable","mv_truthTable"),
                              keepTruthTable=TRUE, all.sol = FALSE, ...)
{
    call <- match.call()
    call[[1]] <- as.name("reduce")
    explain <- match.arg(explain)
    remainders <- match.arg(remainders)
    contradictions <- match.arg(contradictions)
    dontcare <- match.arg(dontcare)
    conddontcare <- match.arg(cdontcare) ## dontcare in condition
    preprocess <- match.arg(preprocess)
    dots <- list(...)
    x <- do.call(preprocess,c(list(mydata=x, outcome=outcome,conditions=conditions),dots))
    ans <- do.call(reduce.truthTable,list(x=x,explain=explain,remainders=remainders,
                                          contradictions=contradictions,dontcare=dontcare,cdontcare=cdontcare,
                                          keepTruthTable=keepTruthTable,all.sol=all.sol, dots))
    ans$call <- call
    ans
}

reduceOld <- function(mydata,outcome,conditions,
                   explain=c("positive","negative"),
                   remainders=c("exclude","include"),
                   contradictions=c("remainders","positive","negative"),
                   dontcare=c("remainders","positive","negative"),
                   preprocess=c("cs_truthTable","fs_truthTable","pass"),
                   keepTruthTable=TRUE,
                   ...)
{
  ## This is the original version of reduce.default, the result is accuate
  ## The new reduce.default use ereduce1, faster, but need more tests.
  call <- match.call()
  explain <- match.arg(explain)
  contradictions <- match.arg(contradictions)
  remainders <- match.arg(remainders)
  dontcare <- match.arg(dontcare)
  if (!"truthTable" %in% class(mydata)){
    preprocess <- match.arg(preprocess)
    dots <- list(...)
    mydata <- do.call(preprocess,c(list(mydata=mydata,nlevels=nlevels,outcome=outcome,conditions=conditions),dots))
    mydata <- mydata$truthTable
    colmax <- sapply(mydata[,conditions],max,na.rm=T)
    if (any(colmax+1 > nlevels)) stop("Mismatch of values of conditions and 'nlevels' argument.")
  } else mydata <- mydata$truthTable

  ##  if (keepTruthTable) truthTable <- subset(mydata,OUT!="?") else truthTable <- NULL
  if (keepTruthTable) {
    truthTable <- mydata[mydata[["OUT"]]!="?",] ## subset(mydata,OUT!="?")
    ## to avoid unbined global variable of OUT, do not use subset(mydata, OUT...)
  } else {truthTable <- NULL }
  ##if (explain=="positive") explained <- subset(mydata,OUT=="1",conditions) ## dat1
  ##if (explain=="negative") explained <- subset(mydata,OUT=="0",conditions) ## dat0
  if (dontcare=="remainders") mydata <- mydata[mydata[["OUT"]]!="-9",] ## subset(mydata,OUT!="-9" )
  if (dontcare=="positive") mydata[['OUT']][mydata[['OUT']]=="-9"] <- "1"
  if (dontcare=="negative") mydata[['OUT']][mydata[['OUT']]=="-9"] <- "0"
  dat1 <- mydata[mydata[["OUT"]]=="1",conditions] ## subset(mydata,OUT=="1",conditions)
  dat0 <- mydata[mydata[["OUT"]]=="0",conditions] ## subset(mydata,OUT=="0",conditions)
  datC <- mydata[mydata[["OUT"]]=="C",conditions] ## subset(mydata,OUT=="C",conditions)
  if (contradictions=="positive") dat1 <- rbind(dat1,datC)
  if (contradictions=="negative") dat0 <- rbind(dat0,datC)
  idExclude <- apply(dat0,1,implicant2Id,nlevels=nlevels)
  if (explain=="positive") explained <- dat1
  if (explain=="negative") explained <- dat0
  if (remainders=="include"){
    ## if necessary conditons -> add some remainders to dat0
    superSets1 <- apply(dat1, 1, superSet,nlevels=nlevels)
    dim(superSets1) <- NULL ## set dim to NULL rather than use as.vector to speed it up.
    superSets1 <- unique(superSets1)
    superSets0 <- apply(dat0, 1, superSet,nlevels=nlevels)
    dim(superSets0) <- NULL
    superSets0 <- unique(superSets0)
    if (explain=="positive") primesId <- setdiff(superSets1,superSets0)
    if (explain=="negative") primesId <- setdiff(superSets0,superSets1)
    primesId <- reduce1(primesId,nlevels=nlevels)
  } else if (remainders=="exclude") {
    if (explain=="positive") primesId <- apply(dat1,1,implicant2Id,nlevels=nlevels)
    if (explain=="negative") primesId <- apply(dat0,1,implicant2Id,nlevels=nlevels)
    primesId <- reduce2(primesId,nlevels=nlevels)
  }
  primeImplicants <- id2Implicant(primesId ,nlevels=nlevels,names=conditions)
  ##  attr(primeImplicants,"explained") <- explained ## give it to argument of PIChart directly
  PIChart <- PIChart(primeImplicants,explained)
  sl <- solvePIChart(PIChart,"combn")
  solutions <- apply(sl,2,function(idx)primeImplicants[idx,])
  commonSolutions <- apply(sl,1,function(idx) {if (length(id <- unique(idx))==1) id })
  ans <- list(solutions=solutions,commonSolutions=commonSolutions,solutionsIDX=sl,primeImplicants=primeImplicants,
              truthTable=truthTable,explained=explained,idExclude=idExclude,nlevels=nlevels,PIChart=PIChart,
              call=call)
  class(ans) <- c("QCA")
  ans
}


toString <- function(implicant, traditional,nlevels,name){
  ## turn each implicant into a string
  ## nm <- name[!is.na(implicant)]
  ## implicant <- implicant[!is.na(implicant)]
  nm <- name[!is.dontcare(implicant)]
  implicant <- implicant[!is.dontcare(implicant)]
  if (traditional && all(nlevels==2)) {
    nm[implicant==1] <- toupper(nm[implicant==1])
    nm[implicant==0] <- tolower(nm[implicant==0])
    res <- paste(nm,sep="",collapse="*")
  } else {
    res <- paste(nm,sprintf("{%s}",implicant),sep="",collapse="*")
  }
  res
}


prettyPI <- function(object,traditional=TRUE,...){
  var_names <- names(object$explained)
  nlevels <- object$nlevels
  solutions <- object$solutions

  toPI <- function(solution){
    if (is.null(solution)) {
      ans <- list(PI="",N=0)
    } else {
      PIs <- apply(solution,1,toString,traditional=traditional,nlevels=nlevels,name=var_names)
      PI <- paste(PIs,collapse=" + ")
      ans <- list(PI=PI,N=length(PIs))
    }
  }

  ans <- lapply(solutions,toPI)
  ans
} ## end of prettyPI()

print.QCA <- function(x,traditional=TRUE,show.truthTable=FALSE,...){
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    PIs <- prettyPI(x,traditional=traditional)
    Nec <- commonConfiguration(x,traditional=traditional)
    if (!is.null(truthTable <- x$truthTable) && show.truthTable){
        cat(sprintf("truthTable with %i configuration(s)\n\n",nrow(truthTable)))
        print(truthTable)
    }
    cat(sprintf("\n----------------\nExplaining %i configuration(s)\n",nrow(x$explained)))
    for (i in seq_len(length(PIs))) {
        cat("\n----------------\n")
        cat(sprintf("Prime implicant No. %i with %i implicant(s)\n\n",i,PIs[[i]]$N))
        writeLines(strwrap(PIs[[i]]$PI))
        cat(sprintf("\nCommon configuration: %s\n",Nec[[i]]))
    }
}


summary.QCA <- function(object,traditional=TRUE,show.case=TRUE,...){
    ## summary of coverage
    ## make use of rownames of truthTable and explained components to match each other
    ## but the actualy meaning of rownames doesn't matter.
    ## cases covered by multiple PIs???
    explain <- object$call$explain
    if (is.null(explain)) explain <- "positive" ## will be null if default
    truthTable <- object$truthTable
    PIs <- prettyPI(object,traditional=traditional)
    Cases <- truthTable[rownames(truthTable) %in% rownames(object$explained), "Cases"]
    OUT <- truthTable[rownames(truthTable) %in% rownames(object$explained), "OUT"]
    ##  cidx <- OUT=="C"
    ##  Cases[cidx] <- paste("(",Cases[cidx],")",sep="") ## The contraditory configuration is in (): now in cs_truthTable
    if (show.case){
        if (pmatch(explain,"positive",0)==1) NCase <- truthTable[rownames(truthTable) %in% rownames(object$explained), "freq1"]
        if (pmatch(explain,"negative",0)==1) NCase <- truthTable[rownames(truthTable) %in% rownames(object$explained), "freq0"]
        N_total <- sum(truthTable["NCase"])
        N_positive <- sum(truthTable["freq1"])
        N_negative <- sum(truthTable["freq0"])
        N <- apply(object$PIChart, 1, function(each) sum(each * NCase))
        coverage <- apply(object$solutionsIDX,2,function(each) N[each])
        if (!is.matrix(coverage)) coverage <- as.matrix(coverage)
        ## a matrix, each column represents one solution
        rownames(coverage) <- paste("PI",seq_len(nrow(coverage)),sep=".")
        colnames(coverage) <- paste("S",seq_len(ncol(coverage)),sep=".")
        prop <- coverage/N_total
    }
    cases <- apply(object$solutionsIDX,2,function(each) {
        ByNPIs <- colSums(object$PIChart[each,,drop=FALSE])
        ## cases covered by ByNPIs PIs
        N_duplicated <- sum(NCase*(ByNPIs-1))
        ## cases covered by multiple PIs
        idx <- object$PIChart[each,,drop=FALSE]
        CasesWithN <- paste("(",ByNPIs,")",Cases,sep="")
        ans <- apply(idx,1,function(idx2) paste(CasesWithN[which(idx2)],sep="", collapse=" "))
        ## group cases for each config
        ans <- paste(ans,collapse=" + ")
        res <- c(PI=ans,Ndup=N_duplicated)
        res
    })
    gof <- list()
    for (i in seq(length(object$solutions))) {
        gof[[i]] <- cbind(consistency(object, data=object$data, which=i),
                        coverage(object,data=object$data, type="raw", which=i),
                        coverage(object,object$data,type="unique", which=i))
    }
    ans <- list(N=N_total,N1=N_positive,N0=N_negative,Ndup=as.numeric(cases["Ndup",]),
                Ncoverage=coverage,prop.total=prop,PIs=PIs,call=object$call,cases=cases["PI",], gof=gof)
    class(ans) <- "summary.QCA"
    ans
}

print.summary.QCA <- function(x,digits=3,traditional=FALSE,...){
  PIs <- x$PIs
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(PIs)==0) {
    cat("\nNo solution.\n")
  } else {
    cat(sprintf("Total number of cases: %i\n",x$N))
    cat(sprintf("Number of cases [1]: %i\n",x$N1))
    cat(sprintf("Number of cases [0]: %i\n",x$N0))
    if (x$call$explain=="positive" || is.null(x$call$explain)) {
      cat("Explaining [1]\n")
      Nexplained <- x$N1
      } else  {
        cat("Explain cases [0]\n")
        Nexplained <- x$N0
        }
    for (i in seq_len(length(PIs))) {
      cat("\n----------------\n")
      cat(sprintf("Prime implicant No. %i with %i implicant(s)\n\n",i,PIs[[i]]$N))
      writeLines(strwrap(PIs[[i]]$PI))
      cat("\nGoodness of fit", fill=TRUE)
      print(x$gof[[i]])
      cat("\n")
      writeLines(strwrap(sprintf("Number of cases: %s\n",paste(x$Ncoverage[,i],collapse=" + "))))
      ## number of case (sum is not the number of explained cases owning to cases covered by multiple PIs)
      writeLines(strwrap(sprintf("Percentage of explained cases: %s\n",
                  paste(paste(round(x$Ncoverage[,i]/Nexplained*100,digits),"%",sep=""),collapse=" + ")
                                 )
                         )
                 )
      ## prop
      writeLines(strwrap(sprintf("Cases covered by multiple PIs: %s (%s%%)\n",x$Ndup[i],
                                 round(x$Ndup[i]/Nexplained * 100,digits))))
      ## dup case
      writeLines(strwrap(sprintf("Cases: %s \n",x$cases[i])))
    }
  }
}


subCombination <- function(implicant,nlevels=rep(2,length(implicant)))
{
    ## It returns a subset of subSet(), that is subCombination
    ## The difference between the two are: subSets may contain dontcare value
    ## subCombination contains NO dontcare value.
    ## if (any(na.id <- is.na(implicant))){
    if (any(na.id <- is.dontcare(implicant))){
        IDX <- cumprod(nlevels+1)/(nlevels+1)
        idx <- which(na.id)
        nn <- nlevels[idx]
        IDX <- IDX[idx]
        exp <- sprintf("c(1:%i)",nlevels[idx])
        dat <- eval(parse(text = sprintf("expand.grid(%s)",paste(exp,sep="",collapse=","))))
        ans <- apply(dat,1,function(x) sum(IDX*x))
        ids <- implicant2Id(implicant,nlevels=nlevels) + ans
        ids
    } else {
        ids <- implicant2Id(implicant,nlevels=nlevels)
        ids
    }
}

SA <- simplifyingAssumption <- function(object,...){
## object is from reduce()
  nlevels <- object$nlevels
  id_explaind <- apply(object$explained,1,implicant2Id,nlevels=nlevels)
  Var_names <- names(object$explained) ## may be improved
  ans <-  lapply(object$solutions,function(each) unique(as.vector(apply(each,1,subCombination,nlevels=nlevels))))
  IDs <- lapply(ans,function(each) {
    res <- unique(unlist(each))
    res <- res[ !res %in% id_explaind]
    res})
  if (is.list(IDs)) {
      solutions <- lapply(IDs,function(each) id2Implicant(each,nlevels=nlevels,names=Var_names))
    } else  solutions <- id2Implicant(IDs,nlevels=nlevels,names=Var_names)
   object$solutions=solutions
   object$SAIDs <- IDs
   object$commonSolutions <- object$solutionsIDX <- object$truthTable <- NULL
   class(object) <- c("SA","QCA")
   object
}

CSA <- function(object1,object0){
  ## contradictory simplifying assumptions
  ## x1: explain="positive",remaind="include"
  ## x2: explain="negative",remaind="include"
  if (length(object1$solutions)>1) {
    stop("object1 contatins more than 1 solution.")
  }
  if (length(object0$solutions)>1) {
    stop("object0 contatins more than 1 solution.")
  }
  if (! "SA" %in% class(object1)) object1 <- SA(object1)
  if (! "SA" %in% class(object0)) object0 <- SA(object0)
  id1 <- unlist(object1$SAIDs)
  id0 <- unlist(object0$SAIDs)
  ids <- intersect(id1,id0)
  nlevels <- object1$nlevels
  names <- names(object1$explained)
  if (length(ids)>0) {
    solutions <- id2Implicant(ids,nlevels=nlevels ,names=names)
    object1$solutions <- list(solutions)
  } else {
    object1$solutions <- list(NULL)
  }
   object1$commonSolutions <- object1$solutionsIDX <- object1$truthTable <- object1$SAIDs <- NULL
   object1
}

print.SA <- function (x, traditional = TRUE, ...)
{
    if (length(x$solutions[[1]])==0) cat("No Simplifying Assumption",fill = TRUE) else {
    PIs <- prettyPI(x, traditional = traditional)
    Nec <- commonConfiguration(x, traditional = traditional)
    cat("Simplifying Assumptions",fill = TRUE)
    for (i in seq_len(length(PIs))) {
        cat("\n----------------\n")
        cat(sprintf("Prime implicant No. %i with %i implicant(s)\n\n",
            i, PIs[[i]]$N))
        writeLines(strwrap(PIs[[i]]$PI))
        cat(sprintf("\nCommon configuration: %s\n", Nec[[i]]))
    }
  }
}

'[.QCA' <- function(object,which){
  if (!all(which %in% seq_len(length(object$solutions)))) stop("which is out of range.")
  old_class <- class(object)
  object <- unclass(object)
  idx <- match("solutions",names(object))
  attrs <- object[-idx]
  solutions <- list(solutions=object$solutions[which])
  ans <- c(solutions,attrs)
  class(ans) <- old_class
  ans
}

constrReduce <- function(object,exclude=NULL,include=NULL,necessary=NULL){
  ## get the intermediate solutions
  ## all arguments are data.frame with ncol=length(object$nlevels)
  if (is.null(exclude) && is.null(include) && is.null(necessary)) {
    stop("No constraint is provided.")
  }
  explained <- object$explained
  nlevels <- object$nlevels
  solutions <- object$solutions
  if (length(solutions)>1) stop("There are multiple solutions.You can use '[' to select one.")
  solution <- solutions[[1]]
  ids1 <- unique(unlist(as.vector(apply(solution,1,subCombination ,nlevels=nlevels))))
  ## idsExplained <- apply(explained,1,implicant2Id,nlevels)
  ## ids1 might include remainders, use idsExplained instead
  if (!is.null(exclude)){
    ## double check when there is NA in exclude???
    if (class(exclude)!="data.frame") stop("exclude should be a data.frame.")
    if (any(exclude > nlevels,na.rm=T)) stop("elements of exclude out of range.")
    ## idsExclude <- unlist(as.vector(apply(exclude,1,subSet,nlevels=nlevels)))
    idsExclude <- unlist(as.vector(apply(exclude,1,subCombination,nlevels=nlevels)))
    superSetIdsExclude <- unlist(as.vector(apply(exclude,1,superSet,nlevels=nlevels)))
    ids1 <- setdiff(ids1, idsExclude)
  }
  if (!is.null(include)){
    if (class(include)!="data.frame") stop("include should be a data.frame.")
    if (any(include > nlevels,na.rm=T)) stop("elements of include out of range.")
    ## idsInclude <- unique(unlist(as.vector(apply(include,1,subSet,nlevels=nlevels))))
    idsInclude <- unique(unlist(as.vector(apply(include,1,subCombination,nlevels=nlevels))))
    if (length(ids0 <- intersect(object$idExclude,idsInclude)>0)) {
      warning("Negative cases are included.")
      idsInclude <- setdiff(idsInclude,ids0)
    }
    ids1 <- union(ids1,idsInclude)
  }
  primesId  <- reduce2(ids1,nlevels=nlevels)
  primeImplicants <- id2Implicant(primesId ,nlevels=nlevels,names=names(object$explained))
  if (!is.null(necessary)){
    Var_name <- names(necessary)
    idx <- match(Var_name, names(primeImplicants))
    if (any(is.na(idx))) {warning("Some variables are not conditions.")} ## may improve.
    idx1 <- !is.na(idx)
    idx2 <- idx[!is.na(idx)] ## the position in PIs
    primeImplicants[,idx2] <- unlist(necessary[idx1])
    primeImplicants <- unique(primeImplicants) ## the rownames of PIs is irrelevant here.
    rowid <- apply( primeImplicants,1,implicant2Id,nlevels=nlevels)
    rownames(primeImplicants) <- as.character(rowid)
  }
  sl <- solvePIChart(PIChart(primeImplicants,explained))
  solutions <- apply(sl,2,function(idx)primeImplicants[idx,])
  commonSolutions <- apply(sl,1,function(idx) {if (length(id <- unique(idx))==1) id })
  object$primeImplicants <- primeImplicants
  object$solutions <- solutions
  object$commonSolutions <- commonSolutions
  object$call <- match.call()
  object
}

update.QCA <- function (object, ..., evaluate = TRUE)
{
    call <- object$call
    extras <- match.call(expand.dots = FALSE)$...
    argsList <- c("mydata", "outcome", "conditions", "cutoff1", "cutoff0", "cutoffc",
                  "complete", "weight", "show.cases", "cases", "nlevels", "ncases_cutoff",
                  "consistency_cutoff", "quiet", "explain", "remainders",
                  "contradictions", "dontcare", "preprocess", "keepTruthTable","missing")
    IDX <- pmatch(names(extras),argsList)
    if (any(is.na(IDX))) stop("multiple arguments are matched.")
    names(extras) <- argsList[IDX]
    IDX2 <- pmatch(names(call),argsList)
    names(call) <- argsList[IDX2]
    if (length(extras) > 0) {
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if (any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }
    if (evaluate)
      eval(call, parent.frame())
    else call
  }

implicantsToDF <- function(x, conditions){
    ## implicantsToDF(x="A*S*R+A*C*S*i",conditions=c("A", "C", "S", "I", "R"))
    x <- strsplit(x,"+",fixed=T)[[1]]
    x <- strsplit(x,"*",fixed=T)
    dat <- as.data.frame(matrix(-9,length(x),length(conditions)))
    names(dat) <- toupper(conditions)
    for (i in 1:length(x)){
        xx <- x[[i]]
	idx <- match(toupper(xx),toupper(conditions))
	dat[i,idx] <- as.numeric(xx == toupper(xx))
    }
    dat
}

boolIntersect <- function(..., string=TRUE){
    x <- paste(unlist(list(...)),collapse="+")
    condition <- toupper(unique(strsplit(x,"[*+]")[[1]]))
    nlevel <- rep(2, length(condition))
    DF <- implicantsToDF(x,condition)
    ids <- apply(DF,1,subCombination) ## can be matrix of a list
    if (is.list(ids)) ids <- unlist(ids) else ids <- as.vector(ids)
    tids <- table(ids)
    ids <- names(tids)[tids==nrow(DF)]
    ids <- as.numeric(ids)
    if (length(ids)>0) {
        ans <- id2Implicant(ids,nlevel,names=condition)
        if (string) {
            ans <- apply(ans,1,toString,traditional=TRUE,nlevels=nlevel,name=condition)
            ans <- paste(ans,collapse=" + ") ## ans always has only one implicant?
        }
    } else ans <- NULL
     writeLines(strwrap(ans))
    invisible(ans)
}

thresholdssetter <- function(x, nthreshold=1, method="average", dismethod="euclidean", thresholds=NULL, plot=TRUE){
  ## method -> see mehtod of hclust
    if (is.null(thresholds)){
        nx <- sort(x)
        idx <- cutree(hclust(dist(nx,method=dismethod),method=method),nthreshold+1)
        thresholds <-  sapply(seq_len(nthreshold),FUN=function(each) (max(nx[idx==each])+min(nx[idx==(each+1)]))/2)
    } else {
        if (any(thresholds >= max(x,na.rm=TRUE))) stop("Thresholds are too large.")
        if (any(thresholds < min(x,na.rm=TRUE))) stop("Thresholds are too small.")
    }
    ans <- unclass(cut(x,c(min(x,na.rm=TRUE)-1, thresholds, max(x,na.rm=TRUE)),include.lowest=T)) - 1
    ## use min-1, so it works even the thresholds equal min
    attr(ans,"levels") <- NULL
    attr(ans,"thresholds") <- thresholds
    if (plot) {
        barplot(table(ans), ylab="Number of cases")
        box()
    }
    ans
}

## EssentialPI <- function(PI){
##   ## similar with QCA:::rowDominance2, but more greedy.
##   ## Some prime implicants are redundant because they are already covered by others,
##   ## so it is better to simplly eliminate them.
##   ## This is what rowDominance2() function does: it gets rid of redundant PIs (by Adrian).
##   N <- nrow(PI)
##   sums <- colSums(PI)
##   Min <- min(sums)
##   ans <- sapply(1:N,FUN=function(ii) colSums(PI[-ii,,drop=FALSE]))
##   idx <- apply(ans,1,FUN=function(x) all (x >= Min))
##   PI[idx,,drop=FALSE]
## }

commonConfiguration <- function(object,traditional=TRUE){
    solutions <- object$solutions
    conditions <- names(object$explained)
    is.commonConfiguration <- function(x){
        if (is.null(x)) {
            ans <- "None"
        }
        else {
            ## ans <- apply(x,2,function(idx) length(unique(idx))==1 & !all(is.na(idx)))
            ans <- apply(x,2,function(idx) length(unique(idx))==1 & !all(is.dontcare(idx)))
            if (any(ans)){
                neccond <- conditions[ans]
                val <- x[1,ans] ## values of condition
                if (traditional && all(object$nlevels==2)){
                    neccond[which(val==1)] <- toupper(neccond[which(val==1)])
                    neccond[which(val==0)] <- tolower(neccond[which(val==0)])
                    ans <- paste(neccond,collapse="*")
                } else ans <- paste(neccond,"{",val,"}",sep="",collapse="*")
            } else ans <- "None"
            ans
        }
    }
    res <- lapply(solutions,is.commonConfiguration)
    res
}

which.commonConfiguration <- function(object){
    which(unlist(commonConfiguration(object)!="None"))
}


drop1.truthTable <- function(object,scope,...){
  cl <- object$call
  conds <- as.list(cl$conditions)
  Nconds <- length(conds)-1
  idx <- c()
  for (i in seq(from=2,to=Nconds)){
   condsN <- as.call(conds[-i])
   cl$conditions <- condsN
   ans <- eval(cl, parent.frame())
   if (all(ans$truthTable$OUT!="C")){idx <- c(idx,i)}
 }
  return(unlist(conds[idx]))
  ## condition which is drop can still return truthTable without contradiction
}
