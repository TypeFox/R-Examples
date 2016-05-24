########################################################################################
#  betaRPMM - Recursively Partitioned Mixture Model for Beta Mixtures
#    
#  AUTHOR:  E. Andres Houseman, Sc.D.
#  CREATED:        July, 2008
#  LAST MODIFIED:  March, 2012
#
########################################################################################

# Required library
# library(cluster)

########################################################################################
# FUNCTION:  glcTree
#    Performs gaussian mixture modeling using 
#    recursive partitioning
#
# ARGUMENTS:
#    x:              Data matrix (n x j) on which to perform clustering
#                    Missing values are currently not supported
#    initFunctions:  Function of type "glcInitialize..." for initializing
#                    latent class model.  See glcInitializeFanny() for an
#                    example of arguments and return values
#    weight:         Weight corresponding to the indices passed (see "index").
#                    Defaults to 1 for all indices
#    index:          Row indices of data matrix to include.
#                    Defaults to all (1 to n).
#    wthresh:        Weight threshold for filtering data to children.
#                    Indices having weight less than this value will not be
#                    passed to children nodes.  Default=1E-8.
#    maxlevel:       Maximum depth to recurse.  Default=Inf.
#    verbose:        Level of verbosity.  Default=2 (too much).  0 for quiet.
#    nthresh:        Total weight in node required for node to be a candidate
#                    for splitting.  Nodes with weight less than this value
#                    will never split.  Defaults to 5.
#    ICL:            Boolean:  use ICL-BIC (i.e. consider entropy in BIC)?
#                              See Ji et al. (2005) for definition of ICL-BIC
#    env:            Object of class "glcTree" to store tree data.
#                    Defaults to a new object.  USER SHOULD NOT SET THIS.
#    nodename:       Name of object that will reprsent node in tree data object.
#                    Defaults to "root".  USER SHOULD NOT SET THIS.
#    level:          Current level.  Defaults to 0.  USER SHUOLD NOT SET THIS.
#    unsplit:        Latent class parameters from parent, to store in current node.
#                    Defaults to NULL for root.  This is used in plotting functions.
#                    USER SHOULD NOT SET THIS.
#
# RETURNS:   Object of class "glcTree"
#
# NOTES:
#   The class "glcTree" is currently implemented as an environment object with
#   nodes represented flatly, with name indicating positition in hierarchy
#   (e.g. "rLLR" = "right child of left child of left child of root")
#   This implementation is to make certain plotting and update functions simpler
#   than would be required if the data were stored in a more natural "list of list"
#   format.
#
#   The following error may appear during the course of the algorithm:
#      Error in optim(logab, betaObjf, ydata = y, wdata = w, weights = weights,  : 
#           non-finite value supplied by optim
#   This is merely an indication that the node being split is too small, in which case
#      the splitting will terminate at that node; in other words, it is nothing 
#      to worry about.
#
########################################################################################

glcTree <- function(x, initFunctions=list(glcInitializeSplitFanny(nu=1.5)),
  weight=NULL, index=NULL, wthresh=1E-8, nodename="root", 
  maxlevel=Inf, verbose=2, nthresh=5, 
  level=0, env=NULL, unsplit=NULL, splitCriterion=glcSplitCriterionBIC){

  n <- dim(x)[1]
  if(is.null(env)) env <- new.env()
  if(is.null(weight)) weight <- rep(1,n)
  if(is.null(index)) index <- 1:n

  if(verbose>0) print(nodename)  

  node <- glcSplit(x, initFunctions, weight, index, level,
     wthresh, verbose=verbose, nthresh=nthresh, splitCriterion=splitCriterion)
  
  node$unsplit <- unsplit
  node$split <- (node$split & (level<maxlevel))

  env[[nodename]] <- node

  if(node$split){
     if(nodename=="root") nodename <- "r"
     nodeleft <- paste(nodename,"L",sep="")
     noderight <- paste(nodename,"R",sep="")

     glcTree(node$x, initFunctions, node$ww[,1], node$index, wthresh,
        env=env, nodename=nodeleft, level=level+1, maxlevel=maxlevel, 
        verbose=verbose, nthresh=nthresh, 
        unsplit=list(mu=node$lco$mu[1,],sigma=node$lco$sigma[1,]), splitCriterion=splitCriterion)
     glcTree(node$x, initFunctions, node$ww[,2], node$index, wthresh,
        env=env, nodename=noderight, level=level+1, maxlevel=maxlevel, 
        verbose=verbose, nthresh=nthresh,
        unsplit=list(mu=node$lco$mu[2,],sigma=node$lco$sigma[2,]), splitCriterion=splitCriterion)
  }

  class(env) <- "glcTree"
  env
}

########################################################################################
# FUNCTION:  print.glcTree
#    Print method for objects of type glcTree
#
# ARGUMENTS:
#    tr:  Object to print
########################################################################################

print.glcTree <- function(x,...){
  cat("Recursively partitioned beta mixture model:")
  cat("\t",length(objects(x)),"nodes, ")
  trms <- unlist(glcTreeApply(x,function(u) 1, terminalOnly=TRUE))
  if(length(trms)>0) trms <- sum(trms)
  else trms <- 0
  cat(trms, "terminal nodes.\n")
}

########################################################################################
# FUNCTION:  plot.glcTree
#    Plot method for objects of type glcTree
#       Plots profiles of terminal nodes in color.
#
# MODIFIED 12-Sep-2014 to address changes in specification of graphical parameters to 'rect' function
#
# ARGUMENTS:
#    env:       Object to print
#    method:    "weight" or "binary"
#      "weight":  width of column represents weight in node
#      "binary":  width of column represents depth of node (narrower=>deeper)
#    palette:   Color palette.  Defaults to green-black-red.
#    xorder:    Order of variables (can be useful for constant ordering across plots)
#       Defaults to ordering by hierarchical clustering (Euc. metric, average linkage).
#    dimensions:  Dimensions of source data to show.  Defaults to all.  Useful for subsets.
#    labelType:  "LR" or "01".
#
########################################################################################
plot.glcTree <- function(x, ...) plotImage.glcTree(x, ...)

plotImage.glcTree <- function (env, start = "r", method = "weight", palette = colorRampPalette(c("yellow", 
    "black", "blue"), space = "Lab")(128), divcol = "red", xorder = NULL, 
    dimensions = NULL, labelType = "LR", muColorEps = 1e-08) {
    
    start2 <- ifelse(start == "r", "root", start)
    xdat <- env[[start2]]$x
    if (is.null(dimensions)) 
        dimensions <- 1:(dim(xdat)[2])
    if (is.null(xorder)) 
        xorder <- hclust(dist(t(xdat[, dimensions])), method = "average")$order
    nodes <- setdiff(objects(env), "root")
    levs <- max(sapply(nodes, nchar)) - nchar(start)
    offset <- nchar(start)
    K <- length(dimensions)
    if (method == "binary") {
        QQ <- levs - offset + 1
        QQ2 <- 2^QQ
        image(0:QQ2, 0:K, matrix(0, QQ2, K), xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", col = "white")
        placement <- function(s, tree) {
            Q <- levs
            y <- strsplit(gsub("R", "1", gsub("L", "0", s)), 
                "")[[1]][-(1:offset)]
            if (length(y) < Q) 
                y <- c(y, rep("0", Q - length(y)))
            sum(as.numeric(y) * 2^(levs:1 - 1))
        }
        places <- c(unlist(glcTreeApply(env, placement, asObject = FALSE, 
            terminalOnly = TRUE)), QQ2)
        nmplaces <- names(places)
    }
    else if (method == "weight") {
        placement <- function(s, tree) {
            sum(s$weight)
        }
        places <- c(0, unlist(glcTreeApply(env, placement, terminalOnly = TRUE)))
        QQ2 <- round(sum(places), 0)
        places <- cumsum(places)
        nmplaces <- names(places)[-1]
        image(0:QQ2, 0:K, matrix(0, QQ2, K), xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", col = "white")
        print(places)
    }
    ncol <- length(palette)
    nplaces <- length(places) - 1
    lpos <- rep(NA, nplaces)
    for (i in 1:nplaces) {
        x0 <- places[i]
        x1 <- places[i + 1]
        muo <- env[[nmplaces[i]]]$unsplit
        mu <- muo$mu[dimensions][xorder]
        muMin <- min(mu)
        muColor <- (mu - muMin + muColorEps)/(max(mu) - muMin + 
            muColorEps)
        for (k in 1:K) {
            rect(x0, k - 1, x1, k, border = NA, density = NA, 
                col = palette[ceiling(muColor[k] * ncol)])
            lpos[i] <- (x0 + x1)/2
        }
        abline(v = x1, col = divcol, lwd = 2)
    }
    if (labelType == "LR") 
        labs <- gsub("^r", "", nmplaces[1:nplaces])
    else {
        labs <- gsub("^r", "", nmplaces[1:nplaces])
        labs <- gsub("L", "0", labs)
        labs <- gsub("R", "1", labs)
    }
    axis(1, lpos, labs, las = 2)
    invisible(xorder)
}


########################################################################################
# FUNCTION:  plotTree.glcTree
#    Alternate plot method for objects of type glcTree:  plots tree.
#
# ARGUMENTS:
#    env:       Object to print
#    start:     Node to start.  Default="r" for "root".
#               Use is deprecated, subsumed by passing the result of "glcSubTree"
#    labelFunction:  Function for generating node labels.  See example.
#    labelAllNodes:  TRUE=All nodes will be labeled; FALSE=Terminal nodes only.
#    labelDigits:    ****TO DO****
#    ...:    ****TO DO****
#
########################################################################################

plotTree.glcTree <- function(env,start="r",labelFunction=NULL,
  buff=4,cex=0.9,square=TRUE,labelAllNodes=FALSE,labelDigits=1,...){
  nodes <- setdiff(objects(env),"root")
  levs <- max(sapply(nodes,nchar))-nchar(start)
  levs2 <- 2^levs
  offset <- max(nchar(start)-1,1)
  plot(c(1,2^levs)/(2^levs+1),c(1-buff,levs),type="n",xlab="",ylab="",xaxt="n",yaxt="n")

  placement <- function(s){
     Q <- nchar(s)-offset
     y <- as.numeric(strsplit(gsub("R","1",gsub("L","0",s)),"")[[1]][-(1:offset)])
     xpos <- (1+sum((2^(Q:1-1))*y))/(2^Q+1)
     c(xpos,levs-Q)
  }

  if(is.null(labelFunction)) label <- function(s)s
  else label <- function(s) {
    x <- labelFunction(env[[s]],digits=labelDigits,...)
    if(!is.null(names(x))) 
      x <- paste(paste(names(x),x,sep="="),collapse=",")
    else x <- round(x,labelDigits)
    x   
  }

  pos1 <- placement(start)
  points(pos1[1],pos1[2],pch=19)
  if(labelAllNodes) 
    text(pos1[1],pos1[2],label(ifelse(start=="r","root",start)),srt=90,adj=1,cex=cex)

  f <- function(nn,prnt){

    if(!(nn=="r"|nn==start)) {
      pos1 <- placement(prnt)
      pos2 <- placement(nn)

      if(square){
         lines(c(pos1[1],pos2[1]), c(pos1[2],pos1[2]))
         lines(c(pos2[1],pos2[1]), c(pos1[2],pos2[2]))
      }
      else {
         lines(c(pos1[1],pos2[1]), c(pos1[2],pos2[2]))
         if(env[[nn]]$split) points(pos2[1],pos2[2],pch=19)
      }

      if(!env[[nn]]$split){
         text(pos2[1],pos2[2],label(nn),srt=90,adj=1,cex=cex)
         return()
      } 
      else if(labelAllNodes) text(pos2[1],pos2[2],label(nn),srt=90,adj=1,cex=cex)
    }

    f(paste(nn,"L",sep=""),nn)
    f(paste(nn,"R",sep=""),nn)

  }
  invisible(f(start,gsub("[[:alpha:]]$","",start)))
}

########################################################################################
# FUNCTION:  glcSubTree
#    Subsets a "glcTree" object, i.e. considers the tree whose root is a given node.
#
# ARGUMENTS:
#    tr:  "glcTree" object to subset
#    node:  Name of node to make root.
#
# RETURNS:   A "glcTree" object whose root is the given node of "tr"
#
# NOTES:
#    This function creates a copy of the data
#
########################################################################################

glcSubTree <- function(tr, node){
  nodeExp <- paste("^",node,sep="")
  os <- objects(tr)[grep(nodeExp,objects(tr))]
  os2 <- gsub(nodeExp,"r",os)
  os2[os2=="r"] <- "root"

  tr2 <- new.env()
  class(tr2) <- "glcTree"
  for (i in 1:length(os)) tr2[[os2[i]]] <- tr[[os[i]]]
  tr2
}

########################################################################################
# FUNCTION:  glcTreeApply
#    Applies a function to every [terminal] node in a glcTree object
#
# ARGUMENTS:
#    tr:           "glcTree" object to recurse.
#    f:            Function to apply at every node.
#    start:        Starting node.  Default="root".
#    terminalOnly: TRUE=only terminal nodes, FALSE=all nodes.
#    asObject:     TRUE="f" accepts node as object
#                  FALSE="f" accepts node by node name and glcTree object, f(nn,tr)
#                     The latter is useful for certain operations requiring knowledge
#                     of tree depth.
#    ...           Additional arguments to pass to "f"
#
# RETURNS:   A list of results; names of elements are names of nodes.
#
# NOTES:
#
########################################################################################

glcTreeApply <- function(tr,f,start="root",terminalOnly=FALSE,asObject=TRUE,...){
  env <- new.env()
  env$result <- list()
  env$nodename <- character(0)
  env$n <- 0

  recurse <- function(nn){
    if(tr[[nn]]$split){
      if(!terminalOnly) {
        env$n <- env$n+1
        if(asObject)  env$result[[env$n]] <- f(tr[[nn]],...) 
        else  env$result[[env$n]] <- f(nn,tree=tr,...) 
        env$nodename[env$n] <- nn 
      }   
      recurse(paste(ifelse(nn=="root","r",nn),"L",sep=""))
      recurse(paste(ifelse(nn=="root","r",nn),"R",sep=""))
    }
    else {
      env$n <- env$n+1
      if(asObject) env$result[[env$n]] <- f(tr[[nn]],...) 
      else  env$result[[env$n]] <- f(nn,tree=tr,...) 
      env$nodename[env$n] <- nn
    }
  }

  recurse(start)
  names(env$result) <- env$nodename
  env$result
}

########################################################################################
# FUNCTION:  glcTreeLeafMatrix
#    Gets matrix of class membership weights for terminal nodes
#
# ARGUMENTS:
#    hmObject:     "glcTree" object to recurse.
#    rounding:     # Digits to round.  Note that if this value is too large,
#           then some subjects will have fractional weight and the function will fail.

# RETURNS:   A factor corresponding to class assignments.
#
# NOTES:
#
########################################################################################

glcTreeLeafMatrix <- function(tr, rounding=3){
  N <- length(tr$root$index)
  CC <- data.frame(glcTreeApply(tr, function(u){
    x <- rep(0,N)
    x[u$index] <- round(u$weight, rounding)
    x
  }, terminalOnly=TRUE))
  CC <- as.matrix(CC)
  CC <- (1/apply(CC,1,sum)) * CC
  CC
}


########################################################################################
# FUNCTION:  glcTreeLeafClasses
#    Gets class assignments corresponding to tree
#    (Works only if terminal class membership weights are close to zero or one)
#
# ARGUMENTS:
#    hmObject:     "glcTree" object to recurse.
#    rounding:     # Digits to round.  Note that if this value is too large,
#           then some subjects will have fractional weight and the function will fail.

# RETURNS:   A factor corresponding to class assignments.
#
# ADDITIONAL NOTES:  Use "glcTreeLeafMatrix" if subjects have fractional class membership
#
########################################################################################

glcTreeLeafClasses <- function(tr){
  W = glcTreeLeafMatrix(tr)
  a = apply(W,1,which.max) 
  as.factor(colnames(W)[a])
}

########################################################################################
# FUNCTION:  glcInitializeSplitFanny
#    Creates a function for initializing latent class model 
#     using "fanny"
#
# ARGUMENTS:
#    nu:      Initial "memb.exp" parameter of "fanny" function
#    nufac:   Factor by which to multipy "memb.exp" and retry if "fanny" fails
#    metric:  Metric to use for fanny
#
# RETURNS:   A function that accepts an n x J data matrix and returns an n x 2 weight matrix
#
# NOTES:
#    This function creates another function that performs the initialization.
#    Auxilary parameters for the initialization are passed here.
#
########################################################################################

glcInitializeSplitFanny <- function(nu=2, nufac=0.875, metric="euclidean"){
  function(xdat){

    warn0 <- options()$warn
    options(warn=-1)
    fano <- try(fanny(xdat, 2, memb.exp=nu, metric=metric),silent=TRUE)

    if(inherits(fano,"try-error")){
       wtmp <- runif(dim(xdat)[1])
       fano <- list(member=cbind(wtmp,1-wtmp))
    }
    else while (abs(fano$coeff["normalized"]) < 1e-07){
      nu <- nu*nufac
      fano <- fanny(xdat, 2, memb.exp=nu)
    }
    options(warn=warn0)

    fano$member
  }
}

########################################################################################
# FUNCTION:  glcInitializeSplitHClust
#    Creates a function for initializing latent class model 
#     using hierarchical clustering
#
# ARGUMENTS:
#    metric:  Metric to use for "dist" object passed to "hclust" function
#    method:  Clustering method used in "hclust" function
#
# RETURNS:   A function that accepts an n x J data matrix and returns an n x 2 weight matrix
#
# NOTES:
#    This function creates another function that performs the initialization.
#    Auxilary parameters for the initialization are passed here.
#
########################################################################################

glcInitializeSplitHClust <- function(metric="manhattan",method="ward"){
  function(xdat){

    hc <- hclust(dist(xdat,method=metric), method=method)
    hcc <- cutree(hc,k=2)
    1*cbind(hcc==1,hcc==2)
  }
}

########################################################################################
# FUNCTION:  glcInitializeSplitEigen
#    Creates a function for initializing latent class model 
#     using eigendecomposition
#
# ARGUMENTS:
#    eigendim:  which eigenvector to use
#    assignmentf:  function that converts eigenvector to class probability
#
# RETURNS:   A function that accepts an n x J data matrix and returns an n x 2 weight matrix
#
# NOTES:
#    This function creates another function that performs the initialization.
#    Auxilary parameters for the initialization are passed here.
#    THIS INITIALIZATION IS PREFERRED WHEN THE NUMBER OF SUBJECTS/SPECIMENS IS LARGE
#
########################################################################################

glcInitializeSplitEigen <- function(
  eigendim=1,
  assignmentf=function(s)(rank(s)-0.5)/length(s)){
  function(xdat){
    z = scale(xdat)
    R = var(z,use="pair")
    s = z %*% eigen(R)$vec
    u = assignmentf(s[,eigendim])
    cbind(u,1-u)
  }
}

########################################################################################
# SPLIT CRITERION FUNCTIONS 
#  (to be documented later)
########################################################################################

# Use BIC
glcSplitCriterionBIC = function(llike1, llike2, weight, ww, J, level){
  out = list()
  effN = sum(weight)
  out$bic1 = log(effN)*2*J-2*llike1
  out$bic2 = log(effN)*(4*J+1)-2*llike2
  out$split = (out$bic2 < out$bic1)
  out
}

# Use BIC-like quantity weighted by level
glcSplitCriterionLevelWtdBIC = function(llike1, llike2, weight, ww, J, level){
  out = list()
  effN = sum(weight)
  out$bic1 = level*log(effN)*2*J-2*llike1
  out$bic2 = level*log(effN)*(4*J+1)-2*llike2
  out$split = (out$bic2 < out$bic1)
  out
}


# Use ICL-BIC (i.e. consider entropy in BIC)
# See Ji et al. (2005) for definition of ICL-BIC
glcSplitCriterionBICICL = function(llike1, llike2, weight, ww, J, level){
  out = list()
  effN = sum(weight)
  out$bic1 = log(effN)*2*J-2*llike1
  out$bic2 = log(effN)*(4*J+1)-2*llike2
  out$entropy <- -sum(ifelse(ww==0,0,ww*log(ww)))
  out$split = (out$bic2 + 2*out$entropy < out$bic1)
  out
}

# Use likelihood ratio test p value
glcSplitCriterionLRT = function(llike1, llike2, weight, ww, J, level){
  glcSplitCriterionLRT_ALPHA = 0.05
  out = list()
  out$llike1 = llike1
  out$llike2 = llike2
  out$J = J
  out$ww = ww
  out$weight = weight
  out$degFreedom = 2*J+1
  out$chiSquareStat = 2*(llike2 - llike1)
  out$split = (pchisq(out$chiSquareStat,out$degFreedom)>1-glcSplitCriterionLRT_ALPHA)
  out
}

# Always split, but record all the information as you go
glcSplitCriterionJustRecordEverything = function(llike1, llike2, weight, ww, J, level){
  out = list()
  out$llike1 = llike1
  out$llike2 = llike2
  out$J = J
  out$ww = ww
  out$weight = weight
  out$split = TRUE
  out
}

########################################################################################
# FUNCTION:  glcSplit
#    Splits a data set into two beta mixture models
#
# ARGUMENTS:
#    x:              Data matrix (n x j) on which to perform clustering
#                    Missing values are currently not supported
#    initFunctions:  Function of type "glcInitialize..." for initializing
#                    latent class model.  See glcInitializeFanny() for an
#                    example of arguments and return values
#    weight:         Weight corresponding to the indices passed (see "index").
#                    Defaults to 1 for all indices
#    index:          Row indices of data matrix to include.
#                    Defaults to all (1 to n).
#    wthresh:        Weight threshold for filtering data to children.
#                    Indices having weight less than this value will not be
#                    passed to children nodes.  Default=1E-8.
#    ICL:            Boolean:  use ICL-BIC (i.e. consider entropy in BIC)?
#                              See Ji et al. (2005) for definition of ICL-BIC
#    verbose:        Level of verbosity.  Default=2 (too much).  0 for quiet.
#    nthresh:        Total weight in node required for node to be a candidate
#                    for splitting.  Nodes with weight less than this value
#                    will never split.  Defaults to 5.
#
# RETURNS:   A list of objects representing split
#
#
########################################################################################

glcSplit <- 
function (x, initFunctions, weight = NULL, index = NULL, level = 0, 
    wthresh = 1e-09, verbose = TRUE, nthresh = 5, splitCriterion = glcSplitCriterionBIC) 
{
    n <- dim(x)[1]
    if (is.null(weight)) 
        weight <- rep(1, n)
    if (is.null(index)) 
        index <- 1:n
    flag <- weight > wthresh
    xdat <- x[flag, , drop = FALSE]
    index <- index[flag]
    weight <- weight[flag]
    n <- dim(xdat)[1]
    obs = !is.na(xdat)
    sumweight <- apply(weight*obs,2,sum)
    msumweight <- mean(sumweight)
    if (!is.null(splitCriterion)) {
        J <- dim(xdat)[2]
        llike1 <- 0
        for (j in 1:J) {
            wY <- weight * xdat[, j]
            mu <- sum(wY, na.rm=TRUE)/sumweight[j]
            sigma <- pmax(1e-08, sqrt(sum(wY * xdat[, j], na.rm=TRUE)/sumweight[j] - 
                mu * mu))
            llike1 <- try(llike1 + sum(weight * dnorm(xdat[, 
                j], mu, sigma, log = TRUE), na.rm = TRUE),silent=TRUE)
            if (inherits(llike1, "try-error")) 
                browser()
        }
    }
    out <- list(index = index, weight = weight, x = xdat)
    lcoList <- list()
    nfits <- 0
    if (nthresh < msumweight) 
        for (s in 1:length(initFunctions)) {
            lco <- try(glc(xdat, w = initFunctions[[s]](xdat), 
                weights = weight, verbose = (verbose > 1)),silent=TRUE)
            if (!inherits(lco, "try-error")) {
                nfits <- nfits + 1
                lcoList[[nfits]] <- lco
            }
        }
    if (nfits > 0) {
        likes <- unlist(lapply(lcoList, function(u) u$llike))
        if (verbose > 0) 
            print(likes)
        lco <- lcoList[[which.max(likes)]]
    }
    else {
        out$split <- FALSE
        return(out)
    }
    out$lco <- lco
    out$split <- TRUE
    out$ww <- weight * lco$w
    if (!is.null(splitCriterion)) {
        out$splitInfo <- splitCriterion(llike1, lco$llike, weight, 
            out$ww, J, level)
        out$split <- out$splitInfo$split
    }
    out
}


########################################################################################
# FUNCTION:  glcTreeOverallBIC
#    Computes the BIC for the latent class model represented by terminal nodes.
#
# ARGUMENTS:
#    tr:           "glcTree" object for which to compute overall BIC.
#
# RETURNS:   A double numeric value representing BIC.
#
# NOTES:
#
########################################################################################

glcTreeOverallBIC <- function(tr,ICL=FALSE){

  C <- unlist(glcTreeApply(tr,function(u,tree)u, 
     asObject=FALSE,terminalOnly=TRUE))
  K <- length(C)
  J <- dim(tr$root$lco$mu)[2]
  L <- array(0,dim=c(length(tr$root$index),K,J))

  W <- array(0,dim=c(length(tr$root$index),K))
  dimnames(L) <- list(NULL, C, NULL)
  dimnames(W) <- list(NULL, C)

  f1 <- function(u,tree){
     n <- dim(tree[[u]]$x)[1]
     for(i in 1:n){
       L[tree[[u]]$index[i],u,] <<- 
        dnorm(tree[[u]]$x[i,],
            tree[[u]]$unsplit$mu,tree[[u]]$unsplit$sigma,log=TRUE)
     } 
     W[tree[[u]]$index,u] <<- tree[[u]]$weight
  }
  glcTreeApply(tr,f1,asObject=FALSE,terminalOnly=TRUE)

  N <- sum(W)
  eta <- apply(W,2,sum)/N
  N <- round(N,0)
  like <- adjConst <- rep(0,N)

  for(i in 1:N){
     if(dim(L)[2]==1) Ltmp <- matrix(L[i,,],nrow=1)
     else Ltmp <- L[i,,]
     lltmp <- apply(Ltmp,1,sum)
     adjConst[i] <- max(lltmp)
     like[i] <- eta %*% exp(lltmp-adjConst[i])
  }
  if(ICL){
    entrop <- -2*sum(ifelse(W==0,0,W*log(W)))
  }
  else entrop <- 0


  log(N)*(2*K*J+K-1)-2*sum(log(like)+adjConst) + entrop
}

########################################################################################
########################################################################################
########################################################################################

########################################################################################
# FUNCTION:  glc
#    Fits a beta mixture model for any number of classes
#
# ARGUMENTS:
#    Y:              Data matrix (n x j) on which to perform clustering
#                    Missing values are currently not supported
#    w:              Initial weight matrix representing classification
#    maxiter:        Maximum number of EM iterations
#    tol:            Convergence tolerance
#    weights:        Case weights
#    verbose:        Verbose output?
#
# RETURNS:   A list of parameters representing mixture model fit, including
#        posterior weights and log-likelihood
#
#
########################################################################################

glc <- function (Y, w, maxiter = 100, tol = 1e-06, weights = NULL, verbose = TRUE) 
{
    J <- dim(Y)[2]
    K <- dim(w)[2]
    n <- dim(w)[1]
    obs <- !is.na(Y)

    if (n != dim(Y)[1]) 
        stop("Dimensions of w and Y do not agree")
    if (is.null(weights)) 
        weights <- rep(1, n)
    mu <- sigma <- matrix(Inf, K, J)
    crit <- Inf
    for (i in 1:maxiter) {
        warn0 <- options()$warn
        options(warn = -1)
        eta <- apply(weights * w, 2, sum)/sum(weights)
        mu0 <- mu
        for (k in 1:K) {
            wt <- w[, k] * weights
            #nn <- sum(wt)
            nn = apply(wt*obs,2,sum)

            wY <- wt * Y
            wY2 <- wY * Y
            mu[k, ] <- apply(wY, 2, sum, na.rm = TRUE)/nn
            sigma[k, ] <- pmax(1e-08, sqrt(apply(wY2, 2, sum, 
                na.rm = TRUE)/nn - mu[k, ] * mu[k, ]))
        }
        ww <- array(NA, dim = c(n, J, K))
        for (k in 1:K) {
            for (j in 1:J) {
                ww[, j, k] <- dnorm(Y[, j], mu[k, j], sigma[k, 
                  j], log = TRUE)
            }
        }
        options(warn = warn0)
        w <- apply(ww, c(1, 3), sum, na.rm = TRUE)
        wmax <- apply(w, 1, max)
        for (k in 1:K) w[, k] <- w[, k] - wmax
        w <- t(eta * t(exp(w)))
        like <- apply(w, 1, sum)
        w <- (1/like) * w
        llike <- weights * (log(like) + wmax)
        crit <- max(abs(mu - mu0),na.rm=TRUE)
        if (verbose) 
            print(crit)
        if (crit < tol) 
            break
    }
    return(list(mu = mu, sigma = sigma, eta = eta, w = w, llike = sum(llike)))
}


########################################################################################
# FUNCTION:  betaEstMultiple
#    Maximum likelihood estimator for beta model on matrix of values
#    (columns having different, independent beta distributions)
#
# ARGUMENTS:
#    Y:              data matrix
#    weights:        case weights
#
# RETURNS:  list of beta parameters and BIC
#
########################################################################################

gaussEstMultiple <- function(Y, weights=NULL){
  J <- dim(Y)[2]
  n <- dim(Y)[1]
  mu <- sigma <- numeric(J)

  one <- rep(1,n)
  if(is.null(weights)) weights <- one

  like <- 0
  for(j in 1:J) {
    nn <- sum(weights)
    wY <- weights*Y[,j]
    wY2 <- wY*Y[,j]
    mu[j] <- sum(wY)/nn
    sigma[j] <- pmax(1E-8,sqrt(sum(wY2)/nn - mu[j]*mu[j]) )

    like <- like + sum(dnorm(Y[,j],mu[j],sigma[j],log=TRUE))
  }       

  bic <- 2*J*log(n) - 2*like

  list(mu=mu, sigma=sigma, bic=bic)
}


