####
#### Phylogenetic ordination tools
####
#### Thibaut Jombart 2008 (tjombart@imperial.ac.uk)

################
# Function ppca
################
ppca <- function(x, prox=NULL, method=c("patristic","nNodes","oriAbouheif","Abouheif","sumDD"), a=1,
                 center=TRUE, scale=TRUE, scannf=TRUE, nfposi=1, nfnega=0){

    ## handle arguments
    ## if(!require(ade4)) stop("The package ade4 is not installed.")
    if (is.character(chk <- checkPhylo4(x))) stop("bad phylo4d object: ",chk)
    ##if (is.character(chk <- checkData(x))) stop("bad phylo4d object: ",chk) : no longer needed

    tre <- as(x, "phylo4")
    method <- match.arg(method)
    NEARZERO <- 1e-10

    ## proximity matrix
    if(is.null(prox)){ # have to compute prox
        W <- proxTips(x, tips="all", method=method, a=a, normalize="row", symmetric=TRUE)
    } else { # prox is provided
        W <- as.matrix(prox)
        if(!is.matrix(W)) stop("W is not a matrix")
        if(ncol(W) != nrow(W)) stop("W is not a square matrix")
        diag(W) <- 0
        W <- 0.5 * (t(W) + W) # re-symmetrization
    }

    N <- nTips(x)

    ## data matrix X
    X <- tdata(x, type="tip")
    X.colnames <- names(X)
    X.rownames <- row.names(X)
    temp <- sapply(X, is.numeric)
    if(!all(temp)) {
        warning(paste("non-numeric data are removed:", X.colnames[!temp]))
        X <- X[,temp]
        X.colnames <- X.colnames[!temp]
        X.rownames <- X.rownames[!temp]
    }

    ## replace NAs
    f1 <- function(vec){
        m <- mean(vec,na.rm=TRUE)
        vec[is.na(vec)] <- m
        return(vec)
    }

    if(any(is.na(X))) {
        warning("Replacing missing values (NA) by mean values")
        X <- as.data.frame(apply(X, 2, f1))
    }

    X <- scalewt(X, center=center, scale=scale) # centring/scaling of traits


    ## main computation ##

    ## make a skeleton of dudi
    res <- dudi.pca(X, center=center, scale=scale, scannf=FALSE,nf=2)
    Upca <- as.matrix(res$c1)

    ## computations of the ppca
    X <- as.matrix(X)
    decomp <- eigen( ((t(X) %*% W %*% X)/N), symmetric=TRUE)
    U <- decomp$vectors # U: principal axes
    lambda <- decomp$values

    ## remove null eigenvalues and corresponding vectors
    toKeep <- (abs(lambda) > NEARZERO)
    lambda <- lambda[toKeep]
    U <- U[, toKeep]
    p <- ncol(U)

    if(scannf){ # interactive part
        barplot(lambda)
        cat("Select the number of global axes: ")
        nfposi <- as.integer(readLines(n = 1))
        cat("Select the number of local axes: ")
        nfnega <- as.integer(readLines(n = 1))
    }

    nfposi <- max(nfposi, 1)
    nfnega <- max(nfnega, 0)
    posi.idx <- 1:nfposi
    if(nfnega<1) {
        nega.idx <- NULL
    } else {
        nega.idx <- (p-nfnega+1):p
    }

    axes.idx <- unique(c(posi.idx, nega.idx)) # index of kept axes
    U <- U[, axes.idx, drop=FALSE]

    S <- X %*% U # S: scores (=princ. components)
    LS <- W %*% S # LS: lagged scores
    A <- t(Upca) %*% U # A: pca princ. axes onto ppca princ. axes.

    ## build the output
    axes.lab <- paste("PA",axes.idx, sep="")
    scores.lab <- paste("PC",axes.idx, sep="")

    res$cent <- res$norm <- res$co <- NULL # cleaning

    res$eig <- lambda # eigenvalues
    res$nf <- NULL
    res$nfposi <- nfposi
    res$nfnega <- nfnega
    res$kept.axes <- axes.idx

    res$c1 <- as.data.frame(U) # principal axes
    names(res$c1) <- axes.lab
    row.names(res$c1) <- X.colnames

    res$li <-  as.data.frame(S) # scores (princ. components)
    names(res$li) <- scores.lab
    row.names(res$li) <- X.rownames

    res$ls <-  as.data.frame(LS) # lagged scores
    names(res$ls) <- scores.lab
    row.names(res$ls) <- X.rownames

    res$as <- as.data.frame(A) # PCA axes onto pPCA axes
    names(res$as) <- axes.lab
    row.names(res$as) <- paste("PCA axis", 1:nrow(A))

    res$tre <- as(tre,"phylo4") # tree

    res$prox <- W # proximity matrix

    res$call <- match.call() # call

    class(res) <- "ppca"

    return(res)
} # end ppca





#####################
# Function scatter.ppca
#####################
scatter.ppca <- function(x, axes=1:ncol(x$li), useLag=FALSE, ...){
    if(useLag){
        df <- as.data.frame(x$ls)
    } else{
        df <- as.data.frame(x$li)
    }

    if(any(axes < 1 | axes > ncol(x$li)) ) stop("Wrong axes specified.")
    df <- df[, axes, drop=FALSE]

    obj <- phylo4d(x$tre,df)
    args <- list(...)
    if(is.null(args$ratio.tree)){
        args$ratio.tree <- 0.5
    }
    args <- c(obj,args)
    do.call(table.phylo4d, args)

    return(invisible(match.call()))
} # end scatter.ppca





######################
# Function print.ppca
######################
print.ppca <- function(x, ...){
  cat("\t#############################################\n")
  cat("\t# phylogenetic Principal Component Analysis #\n")
  cat("\t#############################################\n")
  cat("class: ")
  cat(class(x))
  cat("\n$call: ")
  print(x$call)
  cat("\n$nfposi:", x$nfposi, "axes-components saved")
  cat("\n$nfnega:", x$nfnega, "axes-components saved")
  cat("\n$kept.axes: index of kept axes")

  cat("\nPositive eigenvalues: ")
  l0 <- sum(x$eig >= 0)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  cat("Negative eigenvalues: ")
  l0 <- sum(x$eig <= 0)
  cat(sort(signif(x$eig, 4))[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  cat('\n')
  sumry <- array("", c(1, 4), list(1, c("vector", "length",
                                        "mode", "content")))
  sumry[1, ] <- c('$eig', length(x$eig), mode(x$eig), 'eigenvalues')
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
  sumry <- array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
  sumry[1, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "principal axes: scaled vectors of traits loadings")
  sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "principal components: coordinates of taxa ('scores')")
  sumry[3, ] <- c("$ls", nrow(x$ls), ncol(x$ls), 'lag vector of principal components')
  sumry[4, ] <- c("$as", nrow(x$as), ncol(x$as), 'pca axes onto ppca axes')

  class(sumry) <- "table"
  print(sumry)

  cat("\n$tre: a phylogeny (class phylo4)")
  cat("\n$prox: a matrix of phylogenetic proximities")

  cat("\n\nother elements: ")
  if (length(names(x)) > 16)
    cat(names(x)[17:(length(names(x)))], "\n")
  else cat("NULL\n")
} #end print.ppca





###############
# summary.ppca
###############
summary.ppca <- function (object, ..., printres=TRUE) {

    ## some checks
    if (!inherits(object, "ppca"))stop("to be used with 'ppca' object")
    ## if(!require(ade4)) stop("The package ade4 is not installed.")


    norm.w <- function(X, w) {
        f2 <- function(v) sum(v * v * w)/sum(w)
        norm <- apply(X, 2, f2)
        return(norm)
    }

    resfin <- list()

    if(printres) {
        cat("\n### Phylogenetic Principal Component Analysis ###\n")
        cat("\nCall: ")
        print(object$call)
    }

    appel <- as.list(object$call)
    ## compute original pca
    X <- object$tab # transformed data
    W <- object$prox

    nfposi <- object$nfposi
    nfnega <- object$nfnega

    dudi <- dudi.pca(X, center=FALSE, scale=FALSE, scannf=FALSE, nf=nfposi+nfnega)
    ## end of pca

    Istat <-    data.frame(attributes(moran.idx(X[,1], W,TRUE)))
    row.names(Istat) <- ""
    resfin$Istat <- Istat

    if(printres) {
        cat("\n== Moran's I statistics ==\n")
        print(Istat)
    }

    ## pca scores
    nf <- dudi$nf
    eig <- dudi$eig[1:nf]
    cum <- cumsum(dudi$eig)[1:nf]
    ratio <- cum/sum(dudi$eig)
    moran <- apply(as.matrix(dudi$l1),2,moran.idx, W)
    res <- data.frame(var=eig,cum=cum,ratio=ratio, moran=moran)
    row.names(res) <- paste("Axis",1:nf)
    if(printres) {
        cat("\n== PCA scores ==\n")
        print(res)
    }

    resfin$pca <- res


    ## ppca scores
    ## ppca is recomputed, keeping all axes
    eig <- object$eig
    nfposimax <- sum(eig > 0)
    nfnegamax <- sum(eig < 0)

    listArgs <- appel[-1]
    listArgs$nfposi <- nfposimax
    listArgs$nfnega <- nfnegamax
    listArgs$scannf <- FALSE

    ppcaFull <- do.call(ppca, listArgs) # ppca with all axes

    ndim <- dudi$rank
    nf <- nfposi + nfnega
    toKeep <- c(1:nfposi,if (nfnega>0) (ndim-nfnega+1):ndim)
    varspa <- norm.w(ppcaFull$li,dudi$lw)
    moran <- apply(as.matrix(ppcaFull$li), 2, moran.idx, W)
    res <- data.frame(eig=eig,var=varspa,moran=moran)
    row.names(res) <- paste("Axis",1:length(eig))

    if(printres) {
        cat("\n== pPCA eigenvalues decomposition ==\n")
        print(res[toKeep,])
    }

    resfin$ppca <- res

    return(invisible(resfin))
} # end summary.ppca





#################
# screeplot.ppca
#################
screeplot.ppca <- function(x,...,main=NULL){

  opar <- par("las")
  on.exit(par(las=opar))

  sumry <- summary(x,printres=FALSE)

  labels <- lapply(1:length(x$eig),function(i) bquote(lambda[.(i)]))

  par(las=1)

  xmax <- sumry$pca[1,1]*1.1
  I0 <- unlist(sumry$Istat[1])
  Imin <- unlist(sumry$Istat[2])
  Imax <- unlist(sumry$Istat[3])

  plot(x=sumry$ppca[,2],y=sumry$ppca[,3],type='n',xlab='Variance',ylab="Phylogenetic autocorrelation (I)",xlim=c(0,xmax),ylim=c(Imin*1.1,Imax*1.1),yaxt='n',...)
  text(x=sumry$ppca[,2],y=sumry$ppca[,3],do.call(expression,labels))

  ytick <- c(I0,round(seq(Imin,Imax,le=5),1))
  ytlab <- as.character(round(seq(Imin,Imax,le=5),1))
  ytlab <- c(as.character(round(I0,1)),as.character(round(Imin,1)),ytlab[2:4],as.character(round(Imax,1)))
  axis(side=2,at=ytick,labels=ytlab)

  rect(0,Imin,xmax,Imax,lty=2)
  segments(0,I0,xmax,I0,lty=2)
  abline(v=0)

  if(is.null(main)) main <- ("Decomposition of pPCA eigenvalues")
  title(main)

  return(invisible(match.call()))
} # end screeplot.ppca





############
# plot.ppca
############
plot.ppca <- function(x, axes = 1:ncol(x$li), useLag=FALSE, ...){

    ## some checks
    if (!inherits(x, "ppca")) stop("Use only with 'ppca' objects.")
    if(any(axes>ncol(x$li) | axes<0)) stop("wrong axes required.")

    ## par / layout
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mar = rep(.1,4))
    layout(matrix(c(1,2,3,4,4,4,4,4,4), ncol=3))

    ## some variables
    tre <- x$tre
    n <- nrow(x$li)

    ## 1) barplot of eigenvalues
    omar <- par("mar")
    par(mar = c(0.8, 2.8, 0.8, 0.8))
    r <- length(x$eig)
    col <- rep("white", r)

    keptAxes <- c( (1:r)[1:x$nfposi], (r:1)[1:x$nfnega]) # kept axes
    if(x$nfposi==0) keptAxes <- keptAxes[-1]
    if(x$nfnega==0) keptAxes <- keptAxes[-length(keptAxes)]
    col[keptAxes] <- "grey"

    repAxes <- gsub("PC","",colnames(x$li)[axes]) # represented axes
    repAxes <- as.numeric(repAxes)
    col[repAxes] <- "black"

    barplot(x$eig, col=col)
    title("Eigenvalues", line=-1)
    par(mar=rep(.1,4))
    box()


    ## 2) decomposition of eigenvalues
    par(mar=c(4,4,2,1))
    screeplot(x,main="Eigenvalues decomposition")
    par(mar=rep(.1,4))
    box()


    ## 3) loadings
    if(length(axes)==1){ # one axis retained
        par(mar=c(2.5,4,2,1))
        dotchart(x$c1[,1], labels=row.names(x$c1), main="Loadings",
                 cex=par("cex")*.66)
        abline(v=median(x$c1[,1]), lty=2)
        par(mar=rep(.1,4))
        box()

    } else{ # at least two axes retained
        s.arrow(x$c1[,axes], sub="Loadings")
    }


    ## 4) scatter plot
    ratioTree <- .6
    cexLabel <- 1
    cexSymbol <- 1

    temp <- try(scatter(x, axes=axes, ratio.tree=ratioTree,
                        cex.lab=cexLabel, cex.sym=cexSymbol,
                        show.node=FALSE, useLag=useLag), silent=TRUE) # try default plot
    scatterOk <- !inherits(temp,"try-error")

    while(!scatterOk){
        ## clear 4th screen
        par(new=TRUE)
        plot(1, type="n",axes=FALSE)
        rect(-10,-10, 10,10,col="white")
        par(new=TRUE)
        if(ratioTree > .25 & cexSymbol <= .7) {
            ratioTree <- ratioTree - .05
        }
        if(cexLabel > .65 & cexSymbol <= .5) {
            cexLabel <- cexLabel - .05
        }
        cexSymbol <- cexSymbol - .05

        temp <- try(scatter(x, axes=axes, ratio.tree=ratioTree,
                        cex.lab=cexLabel, cex.sym=cexSymbol,
                        show.node=FALSE, useLag=useLag), silent=TRUE) # try default plot
        scatterOk <- !inherits(temp,"try-error")
    }

    return(invisible(match.call()))

} # end plot.phylo


### testing
## obj <- phylo4d(read.tree(text=mjrochet$tre),mjrochet$tab)
## x@edge.length= rep(1,length(x@edge.label))
## M = cophenetic.phylo(as(x,"phylo"))
## M = 1/M
## diag(M) <- 0


## ppca1 <- ppca(obj,scannf=FALSE,nfp=1,nfn=0)

## plot(ppca1)
