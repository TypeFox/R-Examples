isat <-
function(y, mc=TRUE, ar=NULL, ewma=NULL, mxreg=NULL,
  iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE, blocks=NULL,
  ratio.threshold=0.8, max.block.size=30,
  vcov.type=c("ordinary", "white", "newey-west"),
  t.pval=0.001, do.pet=FALSE, wald.pval=t.pval, ar.LjungB=NULL,
  arch.LjungB=NULL, normality.JarqueB=NULL,
  info.method=c("sc", "aic", "hq"), include.gum=FALSE,
  include.empty=FALSE, tol=1e-07, LAPACK=FALSE, max.regs=NULL, verbose=TRUE,
  print.searchinfo=TRUE, alarm=FALSE, plot=TRUE)
{

  ##arguments:
  isat.call <- sys.call()
  vcov.type <- match.arg(vcov.type)
  info.method <- match.arg(info.method)
  mod <- arx(y, mc=mc, ar=ar, ewma=ewma, mxreg=mxreg,
    vcov.type=vcov.type, qstat.options=NULL,
    tol=tol, LAPACK=LAPACK, verbose=FALSE, plot=FALSE)
  y <- mod$aux$y
  y.n <- mod$aux$y.n
  y.index <- mod$aux$y.index
  y.index.as.char <- as.character(y.index)
  y.name <- mod$aux$y.name
  mX <- mod$aux$mX #NULL if is.null(mX)
  mXnames <- mod$aux$mXnames #NULL if is.null(mX)
  colnames(mX) <- mXnames
  mXncol <- mod$aux$mXncol
  vcov.type <- mod$aux$vcov.type
  qstat.options <- mod$aux$qstat.options
  if(is.null(mX)){ mxkeep <- NULL }else{ mxkeep <- 1:mXncol }

  ##indicator saturation matrices:
  ISmatrices <- list()

  if(iis){ #impulse indicators
    mIIS <- matrix(0,y.n,y.n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste("iis", y.index.as.char, sep="")
    ISmatrices <- c(ISmatrices,list(IIS=mIIS))
  }

  if(sis){ #step-shift indicators
    mSIS <-matrix(0,y.n,y.n)
    loop.indx <- 1:y.n
    tmp <- function(i){ mSIS[i,1:i] <<- 1 }
    tmp <- sapply(loop.indx,tmp)
    colnames(mSIS) <- paste("sis", y.index.as.char, sep="")
    mSIS <- mSIS[,-1]
    ISmatrices <- c(ISmatrices,list(SIS=mSIS))
  }

  if(tis){ #trend indicators
    mTIS <- matrix(0,y.n,y.n)
    v1n <- seq(1,y.n)
    loop.indx <- 1:y.n
    tmp <- function(i){
      mTIS[c(i:y.n),i] <<- v1n[1:c(y.n-i+1)]
    }
    tmp <- sapply(loop.indx,tmp)
    colnames(mTIS) <- paste("tis", y.index.as.char, sep="")
    mTIS <- mTIS[,-1]
    ISmatrices <- c(ISmatrices,list(TIS=mTIS))
  }

  ##user-defined indicators:
  #if uis is a matrix:
  if(!is.list(uis) && !identical(as.numeric(uis),0)){
    #check nrow(uis):
    uis <- as.matrix(coredata(as.zoo(uis)))
    if(nrow(uis) != y.n) stop("nrow(uis) is unequal to no. of observations")
    ISmatrices <- c(ISmatrices,list(UIS=uis))
  }
  #if uis is a list of matrices:
  if(is.list(uis)){
    #check nrow(uis[[i]]):
    for(i in 1:length(uis)){
      uis[[i]] <- as.matrix(coredata(as.zoo(uis[[i]])))
      if(nrow(uis[[i]]) != y.n){
        stop(paste("nrow(uis[[",i,"]]) is unequal to no. of observations",
          sep=""))
      }
    }
    if(is.null(names(uis))){
      uis.names <- paste("UIS", 1:length(uis), sep="")
      names(uis) <- uis.names
    }else{
      uis.names <- paste("UIS", 1:length(uis), sep="")
      for(i in 1:length(uis)){
        if(names(uis)[i]==""){ names(uis)[i] <- uis.names[i] }
      }
    }
    ISmatrices <- c(ISmatrices,uis)
  }

  ##check blocks:
  if(is.list(blocks)){
    if(length(ISmatrices)!=length(blocks)){
      stop("No. of IS matrices is unequal to length(blocks)")
    }
    blocks.is.list <- TRUE
    ISblocks <- blocks
  }else{
    blocks.is.list <- FALSE
    ISblocks <- list()
  }

  ##loop on ISmatrices:
  ISfinalmodels <- list()
  for(i in 1:length(ISmatrices)){

    ##blocks:
    if(!blocks.is.list){

      ncol.adj <- NCOL(ISmatrices[[i]])

      if(is.null(blocks)){
        blockratio.value <- ncol.adj/(ratio.threshold*ncol.adj - NCOL(mX))
        blocksize.value <- ncol.adj/max.block.size
        no.of.blocks <- max(2,blockratio.value,blocksize.value)
        no.of.blocks <- ceiling(no.of.blocks)
      }else{
        no.of.blocks <- blocks
      }

      blocksize <- ceiling(ncol.adj/no.of.blocks)
      partitions.t2 <- blocksize
      for(j in 1:no.of.blocks){
        if( blocksize*j <= ncol.adj ){
          partitions.t2[j] <- blocksize*j
        }
      }
      #check if last block contains last indicator:
      if(partitions.t2[length(partitions.t2)] < ncol.adj){
        partitions.t2 <- c(partitions.t2, ncol.adj)
      }
      blocksadj <- length(partitions.t2)
      partitions.t1 <- partitions.t2 + 1
      partitions.t1 <- c(1,partitions.t1[-blocksadj])

      tmp <- list()
      for(j in 1:blocksadj){
        tmp[[j]] <- partitions.t1[j]:partitions.t2[j]
      }
      ISblocks[[i]] <- tmp

    } #end if(!blocks.is.list)

    ##gets on each block:
    ISspecific.models <- list()
    #for the future?: ISgums <- list(); ISpaths <- list(); ISterminals.results <- list()
    for( j in 1:length(ISblocks[[i]]) ){

      ##print info:
      if(print.searchinfo){
        cat("\n")
        cat(names(ISmatrices)[i],
          " block ", j, " of ", length(ISblocks[[i]]), "...\n", sep="")
        cat("\n")
      }

      mXis <- cbind(mX,ISmatrices[[i]][, ISblocks[[i]][[j]] ])
#OLD:
#      mXis <- cbind(mX,ISmatrices[[i]][, tmpBlocks[[j]] ])
#      mXis <- cbind(mX,ISmatrices[[i]][,partitions.t1[j]:partitions.t2[j]])
      #future?: mXis <- dropvar(mXis, tol=tol, LAPACK=LAPACK)
      mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
        qstat.options=qstat.options, tol=tol, LAPACK=LAPACK,
        verbose=TRUE, plot=FALSE)
#future?: keep=NULL instead of mxkeep?
#      getsis <- getsm(mod, keep=NULL, t.pval=t.pval,
      getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
        do.pet=do.pet, wald.pval=wald.pval, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        info.method=info.method, include.empty=include.empty,
        max.regs=max.regs, estimate.specific=FALSE,
        print.searchinfo=print.searchinfo, plot=FALSE)

      if(is.null(getsis$specific.spec)){
        ISspecific.models[[j]] <- NULL
      }else{
        ISspecific.models[[j]] <- names(getsis$specific.spec)
#For the future?:
#        ISgums[[j]] <- getsis$gum.mean
#        ISpaths[[j]] <- getsis$paths
#        ISterminals.results[[j]] <- getsis$terminals.results
      }

    } #end for(j in 1:length(ISblocks[[i]]))

    ##print info:
    if(print.searchinfo){
      cat("\n")
      cat("GETS of union of retained ",
        names(ISmatrices)[i], " indicators... \n",
        sep="")
      cat("\n")
    }

    ##if no indicators retained from the blocks:
    if(length(ISspecific.models)==0){
      isNames <- NULL
      ISfinalmodels[[i]] <- NULL
    }

    ##when indicators retained from the blocks:
    if(length(ISspecific.models)>0){

      isNames <- NULL

      #which indicators retained?:
      for(j in 1:length(ISspecific.models)){
        #check if mean is non-empty:
        if(!is.null(ISspecific.models[[j]])){
          isNames <- union(isNames, ISspecific.models[[j]])
        }
      } #end for(j) loop
      isNames <- setdiff(isNames, mXnames)
      #isNamesAll[[i]] <- isNames

      #redo gets with union of retained indicators:
      mXisNames <- c(mXnames,isNames)
      mXis <- cbind(mX,ISmatrices[[i]][,isNames])
      colnames(mXis) <- mXisNames
      mXis <- dropvar(mXis, tol=tol, LAPACK=LAPACK)
      mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
        qstat.options=NULL, tol=tol, LAPACK=LAPACK,
        verbose=TRUE, plot=FALSE)
      getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
        do.pet=do.pet, wald.pval=wald.pval, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        info.method=info.method, include.gum=include.gum,
        include.empty=include.empty, max.regs=max.regs,
        print.searchinfo=print.searchinfo, estimate.specific=FALSE,
        plot=FALSE)
      ISfinalmodels[[i]] <- names(getsis$specific.spec)

    } #end if(length(ISspecific.models > 0)

  } #end for(i) loop (on ISmatrices)

  ##add names to ISblocks:
  names(ISblocks) <- names(ISmatrices)

  ##gets of union of retained impulses:
  if(print.searchinfo){
    cat("\n")
    cat("GETS of union of ALL retained indicators...\n")
    cat("\n")
  }

  ##no final models estimated:
  if(length(ISfinalmodels)==0){
    ISfinalmodels <- NULL
    if(is.null(mX)){ mXis <- NULL }else{
      mXis <- zoo(cbind(mX), order.by=y.index)
      colnames(mXis) <- mXnames
    }
  }

  ##final models estimated:
  if(length(ISfinalmodels)>0){

    mIS <- NULL #matrix

    #which indicators were retained?
    for(i in 1:length(ISfinalmodels)){
      isNames <- NULL
      #check if non-empty:
      if(!is.null(ISfinalmodels[[i]])){
        isNames <- setdiff(ISfinalmodels[[i]], mXnames)
      }
      if(length(isNames)>0){
        tmp <- cbind(ISmatrices[[i]][, isNames ])
        colnames(tmp) <- isNames
        mIS <- cbind(mIS, tmp)
      }
    } #end for loop

    mXis <- dropvar(cbind(mX,mIS), tol=tol, LAPACK=LAPACK)
    mXis <- zoo(mXis, order.by=y.index)
  } #end if(length(ISfinalmodels)>0)

  ##gum and gets:
  y <- zoo(y, order.by=y.index)
  mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
    qstat.options=qstat.options, tol=tol,
    LAPACK=LAPACK, verbose=TRUE, plot=FALSE)
  getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
    do.pet=do.pet, wald.pval=wald.pval, ar.LjungB=ar.LjungB,
    arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
    info.method=info.method, include.empty=include.empty,
    max.regs=max.regs, print.searchinfo=print.searchinfo,
    plot=FALSE)

  ##names of retained impulses, mX colnames:
  ISnames <- setdiff(getsis$aux$mXnames, mXnames)
  if(length(ISnames)==0){ ISnames <- NULL }
  colnames(getsis$aux$mX) <- getsis$aux$mXnames

  ##return:
  getsis$gets.type <- "isat"
  getsis$call <- isat.call
  getsis <- c(list(ISfinalmodels=ISfinalmodels,
    ISnames=ISnames), getsis)
  getsis$aux$t.pval <- t.pval #needed for biascorr
  class(getsis) <- "isat"
  if(alarm){ alarm() }
  if(plot){ plot.isat(getsis, coef.path=TRUE) }
  return(getsis)
}
