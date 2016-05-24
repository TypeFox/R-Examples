#' Joint Random Forest for the simultaneous estimation of multiple related networks
#'
#' MAIN FUNCTION -- > JRF
#' 
#' INPUT
#' 
#' X            list object containing data for each class
#' ntree        number of trees
#' mtry         number of variables to be sampled at each node
#' genes.name   list of gene names 
#' 
#' OUTPUT: importance score of interactions
#'
#'
#' OTHER FUNCTIONS -- > importance  and  JRF_onetarget
#' 
#' importance     compute importance score for an object of class JRF 
#' (this file is a modified version of file importance contained in package randomForest, A. Liaw and M. Wiener (2002))
#' 
#' JRF_onetarget  for each class, model the expression of a target gene as a function of the expression of other genes via random forest. 
#'                class specific tree ensemble are designed to borrow information across them. 
#' (this file is a modified version of file randomForest contained in package randomForest, A. Liaw and M. Wiener (2002))
#'   
#"JRF" <-  function(X, ...)UseMethod("JRF")


importance <- function(x,  scale=TRUE) {
  # --- Function importance is a modified version of function importance from R package randomForest
  
  type=NULL;
  class=NULL;
  if (!inherits(x, "randomForest"))
        stop("x is not of class randomForest")
    classRF <- x$type != "regression"
    hasImp <- !is.null(dim(x$importance)) || ncol(x$importance) == 1
    hasType <- !is.null(type)
    if (hasType && type == 1 && !hasImp)
        stop("That measure has not been computed")
    allImp <- is.null(type) && hasImp
    if (hasType) {
        if (!(type %in% 1:2)) stop("Wrong type specified")
        if (type == 2 && !is.null(class))
            stop("No class-specific measure for that type")
    }
    
    imp <- x$importance
    if (hasType && type == 2) {
        if (hasImp) imp <- imp[, ncol(imp), drop=FALSE]
    } else {
        if (scale) {
            SD <- x$importanceSD
            imp[, -ncol(imp)] <-
                imp[, -ncol(imp), drop=FALSE] /
                    ifelse(SD < .Machine$double.eps, 1, SD)
        }
        if (!allImp) {
            if (is.null(class)) {
                ## The average decrease in accuracy measure:
                imp <- imp[, ncol(imp) - 1, drop=FALSE]
            } else {
                whichCol <- if (classRF) match(class, colnames(imp)) else 1
                if (is.na(whichCol)) stop(paste("Class", class, "not found."))
                imp <- imp[, whichCol, drop=FALSE]
            }
        }
    }
    imp<-imp[,2]
    imp
}


# --- Function JRF_onetarget is a modified version of function randomForest from R package randomForest

"JRF_onetarget" <-
  function(x, y=NULL,  xtest=NULL, ytest=NULL, ntree,
           sampsize,              
           totsize = if (replace) ncol(x) else ceiling(.632*ncol(x)),
           mtry=if (!is.null(y) && !is.factor(y))
             max(floor(nrow(x)/3), 1) else floor(sqrt(nrow(x))),
           replace=TRUE, classwt=NULL, cutoff, strata,
           nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
           maxnodes=NULL,
           importance=FALSE, localImp=FALSE, nPerm=1,
           proximity, oob.prox=proximity,
           norm.votes=TRUE, do.trace=FALSE,
           keep.forest=!is.null(y) && is.null(xtest), corr.bias=FALSE,
           keep.inbag=FALSE, nclasses, ...) {
    
    ww=1/sampsize;
    nclass=mylevels=ipi=sw=NULL
    addclass <- is.null(y)
    classRF <- addclass || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
      warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
    }
    if (classRF && !addclass && length(unique(y)) < 2)
      stop("Need at least two classes to do classification.")
    
    n <- ncol(y)           # number of samples
    p <- nrow(x)/nclasses  # number of variables
    
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    keep.forest=!is.null(y) 
    xtest=NULL; ytest=NULL
    testdat <- !is.null(xtest)
    if (testdat) {
      if (ncol(x) != ncol(xtest))
        stop("x and xtest must have same number of columns")
      ntest <- nrow(xtest)
      xts.row.names <- rownames(xtest)
    }
    
    prox <- proxts <- double(1)
    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")
    if (testdat && any(is.na(xtest))) stop("NA not permitted in xtest")
    if (any(is.na(y))) stop("NA not permitted in response")
    if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")
    
    if (is.data.frame(x)) {
      xlevels <- lapply(x, mylevels)
      ncat <- sapply(xlevels, length)
      ## Treat ordered factors as numerics.
      ncat <- ifelse(sapply(x, is.ordered), 1, ncat)
      x <- data.matrix(x)
      if(testdat) {
        if(!is.data.frame(xtest))
          stop("xtest must be data frame if x is")
        xfactor <- which(sapply(xtest, is.factor))
        if (length(xfactor) > 0) {
          for (i in xfactor) {
            if (any(! levels(xtest[[i]]) %in% xlevels[[i]]))
              stop("New factor levels in xtest not present in x")
            xtest[[i]] <-
              factor(xlevels[[i]][match(xtest[[i]], xlevels[[i]])],
                     levels=xlevels[[i]])
          }
        }
        xtest <- data.matrix(xtest)
      }
    } else {
      ncat <- rep(1, p)
      xlevels <- as.list(rep(0, p))
    }
    maxcat <- max(ncat)
    if (maxcat > 32)
      stop("Can not handle categorical predictors with more than 32 categories.")
    
    addclass <- FALSE
    
    proximity <- addclass
    
    
    
    impout <- matrix(0.0, p*nclasses, 2)
    impSD <- matrix(0.0, p*nclasses, 1)
    #  names(impSD) <- x.col.names
    
    
    
    nsample <- if (addclass) 2 * n else n
    Stratify <- length(n) > 1
    
    nodesize=5;
    nrnodes <- 2 * trunc(n/max(1, nodesize - 4)) + 1
    
    maxnodes=NULL
    if (!is.null(maxnodes)) {
      ## convert # of terminal nodes to total # of nodes
      maxnodes <- 2 * maxnodes - 1
      if (maxnodes > nrnodes) warning("maxnodes exceeds its max value.")
      nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
    }
    
    
    ## Compiled code expects variables in rows and observations in columns.
    # x <- t(x)
    storage.mode(x) <- "double"
    
    xtest <- double(1)
    ytest <- double(1)
    ntest <- 1
    labelts <- FALSE
    nt <- if (keep.forest) ntree else 1
    
    
    nPerm=1
    do.trace=F; oob.prox=F
    corr.bias=FALSE
    keep.inbag=FALSE
    impmat <- double(1)
    replace=T
    
    
    
    if (classRF) {
      cwt <- classwt
      threshold <- cutoff
      error.test <- if (labelts) double((nclass+1) * ntree) else double(1)
      
      
      ####### -- call C function to compute tree ------------------------------------------------------- ###########
      rfout <- .C("classRF",
                  x = x,
                  xdim = as.integer(c(p, n)),
                  y = as.integer(y),
                  nclass = as.integer(nclass),
                  ncat = as.integer(ncat),
                  maxcat = as.integer(maxcat),
                  sampsize = as.integer(sampsize),
                  strata = if (Stratify) as.integer(strata) else integer(1),
                  Options = as.integer(c(addclass,
                                         importance,
                                         localImp,
                                         proximity,
                                         oob.prox,
                                         do.trace,
                                         keep.forest,
                                         replace,
                                         Stratify,
                                         keep.inbag)),
                  ntree = as.integer(ntree),
                  mtry = as.integer(mtry),
                  ipi = as.integer(ipi),
                  classwt = as.double(cwt),
                  cutoff = as.double(threshold),
                  nodesize = as.integer(nodesize),
                  outcl = integer(nsample),
                  counttr = integer(nclass * nsample),
                  prox = prox,
                  impout = impout,
                  impSD = impSD,
                  impmat = impmat,
                  nrnodes = as.integer(nrnodes),
                  ndbigtree = integer(ntree),
                  nodestatus = integer(nt * nrnodes),
                  bestvar = integer(nt * nrnodes),
                  treemap = integer(nt * 2 * nrnodes),
                  nodepred = integer(nt * nrnodes),
                  xbestsplit = double(nt * nrnodes),
                  errtr = double((nclass+1) * ntree),
                  testdat = as.integer(testdat),
                  xts = as.double(xtest),
                  clts = as.integer(ytest),
                  nts = as.integer(ntest),
                  countts = double(nclass * ntest),
                  outclts = as.integer(numeric(ntest)),
                  labelts = as.integer(labelts),
                  proxts = proxts,
                  errts = error.test,
                  inbag = if (keep.inbag)
                    matrix(integer(n * ntree), n) else integer(n), as.double(sw))[-1]
      if (keep.forest) {
        ## deal with the random forest outputs
        max.nodes <- max(rfout$ndbigtree)
        treemap <- aperm(array(rfout$treemap, dim = c(2, nrnodes, ntree)),
                         c(2, 1, 3))[1:max.nodes, , , drop=FALSE]
      }
      if (!addclass) {
        ## Turn the predicted class into a factor like y.
        out.class <- factor(rfout$outcl, levels=1:nclass,
                            labels=levels(y))
        names(out.class) <- x.row.names
        con <- table(observed = y,
                     predicted = out.class)[levels(y), levels(y)]
        con <- cbind(con, class.error = 1 - diag(con)/rowSums(con))
      }
      out.votes <- t(matrix(rfout$counttr, nclass, nsample))[1:n, ]
      oob.times <- rowSums(out.votes)
      if (norm.votes)
        out.votes <- t(apply(out.votes, 1, function(x) x/sum(x)))
      dimnames(out.votes) <- list(x.row.names, levels(y))
      class(out.votes) <- c(class(out.votes), "votes")
      if (testdat) {
        out.class.ts <- factor(rfout$outclts, levels=1:nclass,
                               labels=levels(y))
        names(out.class.ts) <- xts.row.names
        out.votes.ts <- t(matrix(rfout$countts, nclass, ntest))
        dimnames(out.votes.ts) <- list(xts.row.names, levels(y))
        if (norm.votes)
          out.votes.ts <- t(apply(out.votes.ts, 1,
                                  function(x) x/sum(x)))
        class(out.votes.ts) <- c(class(out.votes.ts), "votes")
        if (labelts) {
          testcon <- table(observed = ytest,
                           predicted = out.class.ts)[levels(y), levels(y)]
          testcon <- cbind(testcon,
                           class.error = 1 - diag(testcon)/rowSums(testcon))
        }
      }
      cl <- match.call()
      cl[[1]] <- as.name("randomForest")
      out <- list(call = cl,
                  type = if (addclass) "unsupervised" else "classification",
                  predicted = if (addclass) NULL else out.class,
                  err.rate = if (addclass) NULL else t(matrix(rfout$errtr,
                                                              nclass+1, ntree,
                                                              dimnames=list(c("OOB", levels(y)), NULL))),
                  confusion = if (addclass) NULL else con,
                  votes = out.votes,
                  oob.times = oob.times,
                  classes = levels(y),
                  importance = if (importance)
                    matrix(rfout$impout, p, nclass+2,
                           dimnames = list(x.col.names,
                                           c(levels(y), "MeanDecreaseAccuracy",
                                             "MeanDecreaseGini")))
                  else matrix(rfout$impout, ncol=1,
                              dimnames=list(x.col.names, "MeanDecreaseGini")),
                  importanceSD = if (importance)
                    matrix(rfout$impSD, p, nclass + 1,
                           dimnames = list(x.col.names,
                                           c(levels(y), "MeanDecreaseAccuracy")))
                  else NULL,
                  localImportance = if (localImp)
                    matrix(rfout$impmat, p, n,
                           dimnames = list(x.col.names,x.row.names)) else NULL,
                  proximity = if (proximity) matrix(rfout$prox, n, n,
                                                    dimnames = list(x.row.names, x.row.names)) else NULL,
                  ntree = ntree,
                  mtry = mtry,
                  forest = if (!keep.forest) NULL else {
                    list(ndbigtree = rfout$ndbigtree,
                         nodestatus = matrix(rfout$nodestatus,
                                             ncol = ntree)[1:max.nodes,, drop=FALSE],
                         bestvar = matrix(rfout$bestvar, ncol = ntree)[1:max.nodes,, drop=FALSE],
                         treemap = treemap,
                         nodepred = matrix(rfout$nodepred,
                                           ncol = ntree)[1:max.nodes,, drop=FALSE],
                         xbestsplit = matrix(rfout$xbestsplit,
                                             ncol = ntree)[1:max.nodes,, drop=FALSE],
                         pid = rfout$classwt, cutoff=cutoff, ncat=ncat,
                         maxcat = maxcat,
                         nrnodes = max.nodes, ntree = ntree,
                         nclass = nclass, xlevels=xlevels)
                  },
                  y = if (addclass) NULL else y,
                  test = if(!testdat) NULL else list(
                    predicted = out.class.ts,
                    err.rate = if (labelts) t(matrix(rfout$errts, nclass+1,
                                                     ntree,
                                                     dimnames=list(c("Test", levels(y)), NULL))) else NULL,
                    confusion = if (labelts) testcon else NULL,
                    votes = out.votes.ts,
                    proximity = if(proximity) matrix(rfout$proxts, nrow=ntest,
                                                     dimnames = list(xts.row.names, c(xts.row.names,
                                                                                      x.row.names))) else NULL),
                  inbag = if (keep.inbag) rfout$inbag else NULL)
    } 
    else {
      
      rfout <- .C("regRF",
                  x,
                  y, ww,
                  as.integer(c(totsize, p)),
                  sampsize=as.integer(sampsize), as.integer(totsize),
                  as.integer(nodesize),
                  as.integer(nrnodes),
                  as.integer(ntree),
                  as.integer(mtry),
                  as.integer(c(importance, localImp, nPerm)),
                  as.integer(ncat),
                  as.integer(maxcat),
                  as.integer(do.trace),
                  as.integer(proximity),
                  as.integer(oob.prox),
                  as.integer(corr.bias),
                  ypred = double(n * nclasses),
                  impout = impout,
                  impmat = impmat,
                  impSD = impSD,
                  prox = prox,
                  ndbigtree = integer(ntree),
                  nodestatus = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                  leftDaughter = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                  rightDaughter = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                  nodepred = matrix(double(nrnodes * nt * nclasses), ncol=nt),
                  bestvar = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                  xbestsplit = matrix(double(nrnodes * nt * nclasses), ncol=nt),
                  mse = double(ntree * nclasses),
                  keep = as.integer(c(keep.forest, keep.inbag)),
                  replace = as.integer(replace),
                  testdat = as.integer(testdat),
                  xts = xtest,
                  ntest = as.integer(ntest),
                  yts = as.double(ytest),
                  labelts = as.integer(labelts),
                  ytestpred = double(ntest),
                  proxts = proxts,
                  msets = double(if (labelts) ntree else 1),
                  coef = double(2),
                  oob.times = integer(n),
                  inbag = if (keep.inbag)
                    matrix(integer(n * ntree), n) else integer(1), as.integer(nclasses))[c(16:28, 36:41)]
      #         ## Format the forest component, if present.
      if (keep.forest) {
        max.nodes <- max(rfout$ndbigtree)
        rfout$nodestatus <-
          rfout$nodestatus[1:max.nodes, , drop=FALSE]
        rfout$bestvar <-
          rfout$bestvar[1:max.nodes, , drop=FALSE]
        rfout$nodepred <-
          rfout$nodepred[1:max.nodes, , drop=FALSE]
        rfout$xbestsplit <-
          rfout$xbestsplit[1:max.nodes, , drop=FALSE]
        rfout$leftDaughter <-
          rfout$leftDaughter[1:max.nodes, , drop=FALSE]
        rfout$rightDaughter <-
          rfout$rightDaughter[1:max.nodes, , drop=FALSE]
      }
      cl <- match.call()
      cl[[1]] <- as.name("randomForest")
      #         ## Make sure those obs. that have not been OOB get NA as prediction.
      ypred <- rfout$ypred
      if (any(rfout$oob.times < 1)) {
        ypred[rfout$oob.times == 0] <- NA
      }
      out <- list(call = cl,
                  type = "regression",
                  predicted =0,
                  mse = rfout$mse,
                  rsq = 1 - rfout$mse / (var(y[1,]) * (n-1) / n),
                  oob.times = rfout$oob.times,
                  importance = if (importance) matrix(rfout$impout, p * nclasses, 2) else
                    matrix(rfout$impout, ncol=1),
                  importanceSD=if (importance) rfout$impSD else NULL,
                  localImportance = if (localImp)
                    matrix(rfout$impmat, p, n, dimnames=list(x.col.names,
                                                             x.row.names)) else NULL,
                  proximity = if (proximity) matrix(rfout$prox, n, n,
                                                    dimnames = list(x.row.names, x.row.names)) else NULL,
                  ntree = ntree,
                  mtry = mtry,
                  forest = if (keep.forest)
                    c(rfout[c("ndbigtree", "nodestatus", "leftDaughter",
                              "rightDaughter", "nodepred", "bestvar",
                              "xbestsplit")],
                      list(ncat = ncat), list(nrnodes=max.nodes),
                      list(ntree=ntree), list(xlevels=xlevels)) else NULL,
                  coefs = if (corr.bias) rfout$coef else NULL,
                  y = y,
                  test = if(testdat) {
                    list(predicted = structure(rfout$ytestpred,
                                               names=xts.row.names),
                         mse = if(labelts) rfout$msets else NULL,
                         rsq = if(labelts) 1 - rfout$msets /
                           (var(ytest) * (n-1) / n) else NULL,
                         proximity = if (proximity)
                           matrix(rfout$proxts / ntree, nrow = ntest,
                                  dimnames = list(xts.row.names,
                                                  c(xts.row.names,
                                                    x.row.names))) else NULL)
                  } else NULL,
                  inbag = if (keep.inbag)
                    matrix(rfout$inbag, nrow(rfout$inbag),
                           dimnames=list(x.row.names, NULL)) else NULL)
    }
    
    #      print(rfout$mse)
    class(out) <- "randomForest"
    return(out)
    
  }


# --- MAIN function
"JRF" <-
  function(X, ntree,mtry,genes.name) {

    nclasses<-length(X)
    sampsize<-rep(0,nclasses)
    
    for (j in 1:nclasses) sampsize[j]<-dim(X[[j]])[2]
    
    totsize<-tot<-max(sampsize);
    p<-dim(X[[1]])[1];
    imp<-array(0,c(p,length(genes.name),nclasses))
    
    imp.final<-matrix(0,p*(p-1)/2,nclasses);
    vec1<-matrix(rep(genes.name,p),p,p)
    vec2<-t(vec1)
    vec1<-vec1[lower.tri(vec1,diag=FALSE)]
    vec2<-vec2[lower.tri(vec2,diag=FALSE)]
    
    index<-seq(1,p)
    
    for (j in 1:length(genes.name)){

      covar<-matrix(0,(p-1)*nclasses,tot)             
      y<-matrix(0,nclasses,tot)             
      
      for (c in 1:nclasses)  {
        y[c,seq(1,sampsize[c])]<-as.matrix(X[[c]][j,])
        covar[seq((c-1)*(p-1)+1,c*(p-1)),seq(1,sampsize[c])]<-X[[c]][-j,]
      }
      
      jrf.out<-JRF_onetarget(x=covar,y=y,mtry=mtry,importance=TRUE,sampsize=sampsize,nclasses=nclasses,ntree=ntree)
      
      for (s in 1:nclasses) imp[-j,j,s]<-importance(jrf.out,scale=FALSE)[seq((p-1)*(s-1)+1,(p-1)*(s-1)+p-1)]  #- save importance score for net1
      
      }
      
    # --- Derive importance score for each interaction 
      for (s in 1:nclasses){ 
         imp.s<-imp[,,s]; t.imp<-t(imp.s)
         imp.final[,s]<-(imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2        
      }
    
    out<-cbind(as.character(vec1),as.character(vec2),as.data.frame(imp.final),stringsAsFactors=FALSE)
    colnames(out)<-c(paste0('gene',1:2),paste0('importance',1:nclasses))
    return(out)
    
  }

