########################################################################
# Categorical Network Class Methods
# Pathway Analysis

read.pathway <- function(file.path, pathid) {

  lines <- readLines(file.path)
  n <- length(lines)
  if(is.character(pathid)) { 
    for(ll in lines) {
      ls <- strsplit(ll,"\t")[[1]]
      if(ls[1] == pathid)
        break
    }
    if(length(ls) < 1 || ls[1] != pathid)
      stop("no valid pathway")
  }
  else if(is.numeric(pathid)) {
    pathid <- as.integer(pathid)
    if(pathid < 1 || pathid > n)
      stop("n out of range")
    ls <- strsplit(lines[[pathid]],"\t")[[1]]
    cat(ls[1],"\n")
  } 
 if(length(ls) < 3)
    stop("no valid pathway")
  pathname <- ls[1]
  path <- ls[3:length(ls)]

  return(path)
}

cnPathwayTest <- function(pathwayFile, cdata, classes, numParents=1, numCategories=3, sampleFact=10, kstep=4,model=NULL,echo=TRUE, alpha=0.05, outDir=NULL,...) {
  if(is.null(pathwayFile))
    stop("Valid pathway file is expected")
  lines <- readLines(pathwayFile)
  n <- length(lines)
  if(n <= 0)
    stop(pathwayFile, " is invalid")
  
  if(!is.matrix(cdata) && !is.data.frame(cdata))
    stop("data should be a matrix or data frame")
  if(is.data.frame(cdata))
    cdata <- t(cdata)
  if(!is.numeric(cdata))
    stop("data should be numeric")
  numnodes <- nrow(cdata)
  numsamples <- ncol(cdata)
  if(numnodes < 1 || numsamples < 1)
    stop("Invalid sample")

  lvs <- levels(classes)
  if(!is.factor(classes) || length(lvs) != 2)
    stop("Two-classes factor vector is expected")
  if(length(classes) != numsamples)
    stop("One class value for each sample instance is expected")
  classes[classes==lvs[1]] <- 1
  classes[classes==lvs[2]] <- 2
  classes <- as.integer(classes)
  numsamples1 <- sum(classes==1)
  numsamples2 <- sum(classes==2)

  data <- cnDiscretize(cdata, numCategories, "uniform")

  vpathways <- NULL
  vpvals <- NULL
  for(ll in lines) {
    
    ls <- strsplit(ll,"\t")[[1]]
    if(length(ls) < 3) {
      warning("no valid pathway ", ls)
      next
    }
    pathname <- ls[1]
    path <- ls[3:length(ls)]

    pathid <- NULL
    rnames <- rownames(data)
    for(i in 1:length(path)) {
      id <- which(rnames == path[i])
      if(length(id) != 1)
        next
      str <- rnames[id]
      if(nchar(str) > 16 || length(which(pathid==id[1]))>0)
        next
      pathid <- c(pathid, id)
    }
    path <- as.character(rnames[pathid])

    if(length(path) <= 1) {
      if(echo)
        cat("empty\n")
      next
    }
    
    if(echo)
      cat("[", length(vpvals)+1, "\\", n,"] ", ls[1]," (", length(path),") ... ")
    vpathways <- c(vpathways, pathname)
    
    nodeOrder <- 1:length(path)
    nodeCats <- lapply(1:length(path), function(i) return(1:numCategories))
    names(nodeCats) <- path

    path.data <- data[pathid, ]
    
    catnet.res <- cnSearchOrder(path.data, NULL, numParents, NULL, 0, 
                         nodeOrder=nodeOrder, nodeCats=nodeCats, 
                         parentsPool=NULL, NULL, echo=FALSE)

    if(!is.null(model)) {
      bnet <- catnet.res@nets[[1]]
      if(model == "BIC") 
        bnet <- cnFindBIC(catnet.res)
      else if(model == "AIC") 
        bnet <- cnFindAIC(catnet.res)
      net1 <- cnSetProb(bnet, path.data[,classes==1], nodeCats=nodeCats) 
      net2 <- cnSetProb(bnet, path.data[,classes==2], nodeCats=nodeCats)    
      l1 <- net1@likelihood
      l2 <- net2@likelihood
      ps1 <- cnSamples(net1, sampleFact*numsamples)
      lhisto <- cnLoglik(net1, ps1, NULL, bysample=TRUE)
      pval1 <- 2*min(sum(lhisto<l2)/length(lhisto), sum(lhisto>l2)/length(lhisto))
      ps2 <- cnSamples(net2, sampleFact*numsamples)
      lhisto <- cnLoglik(net2, ps2, NULL, bysample=TRUE)
      pval2 <- 2*min(sum(lhisto<l1)/length(lhisto), sum(lhisto>l1)/length(lhisto))
      pval <- pval1*pval2
    }
    else {
      pw <- 0
      pval <- 0
      k <- 1
      w <- catnet.res@loglik-0.5*log(numsamples)*catnet.res@complexity/numsamples
      w <- exp((w - max(w))*numsamples)
      while(k < length(catnet.res@nets)) {
        
        bnet <- catnet.res@nets[[k]]
        net1 <- cnSetProb(bnet, path.data[,classes==1], nodeCats=nodeCats) 
        net2 <- cnSetProb(bnet, path.data[,classes==2], nodeCats=nodeCats)
        
        l1 <- net1@likelihood
        l2 <- net2@likelihood
        
        ps1 <- cnSamples(net1, sampleFact*numsamples)
        lhisto <- cnLoglik(net1, ps1, NULL, bysample=TRUE)
        pval1 <- 2*min(sum(lhisto<l2)/length(lhisto), sum(lhisto>l2)/length(lhisto))
        ps2 <- cnSamples(net2, sampleFact*numsamples)
        lhisto <- cnLoglik(net2, ps2, NULL, bysample=TRUE)
        pval2 <- 2*min(sum(lhisto<l1)/length(lhisto), sum(lhisto>l1)/length(lhisto))
        
        pval <- pval + w[k]*pval1*pval2
        pw <- pw + w[k]
        k <- k + kstep
      }
      pval <- pval / pw
    }
    
    vpvals <- c(vpvals, pval)

    if(echo)
      if(pval < alpha)
         cat(pval, " *\n")
         else
         cat(pval, "\n")
      
    if(!is.null(outDir) && pval < alpha) {
      save(bnet, file=paste(outDir, "/", pathname, "_catnet.Rda", sep=""))
    }
    
  }
  names(vpvals) <- vpathways
  
  return(vpvals)
}

cnPathwayLR <- function(pathwayFile, cdata, classes, numParents=1, numCategories=3, model=NULL,echo=TRUE, alpha=0.05, outDir=NULL,...) {
  if(is.null(pathwayFile))
    stop("Valid pathway file is expected")
  lines <- readLines(pathwayFile)
  n <- length(lines)
  if(n <= 0)
    stop(pathwayFile, " is invalid")
  
  if(!is.matrix(cdata) && !is.data.frame(cdata))
    stop("data should be a matrix or data frame")
  if(is.data.frame(cdata))
    cdata <- t(cdata)
  if(!is.numeric(cdata))
    stop("data should be numeric")
  numnodes <- nrow(cdata)
  numsamples <- ncol(cdata)
  if(numnodes < 1 || numsamples < 1)
    stop("Invalid sample")

  lvs <- levels(classes)
  if(!is.factor(classes) || length(lvs) != 2)
    stop("Two-classes factor vector is expected")
  if(length(classes) != numsamples)
    stop("One class value for each sample instance is expected")
  classes[classes==lvs[1]] <- 1
  classes[classes==lvs[2]] <- 2
  classes <- as.integer(classes)
  numsamples1 <- sum(classes==1)
  numsamples2 <- sum(classes==2)

  data <- cnDiscretize(cdata, numCategories, "uniform")

  vpathways <- NULL
  vpvals <- NULL
  for(ll in lines) {
    
    ls <- strsplit(ll,"\t")[[1]]
    if(length(ls) < 3) {
      warning("no valid pathway ", ls)
      next
    }
    pathname <- ls[1]
    path <- ls[3:length(ls)]

    pathid <- NULL
    rnames <- rownames(data)
    for(i in 1:length(path)) {
      id <- which(rnames == path[i])
      if(length(id) != 1)
        next
      str <- rnames[id]
      if(nchar(str) > 16 || length(which(pathid==id[1]))>0)
        next
      pathid <- c(pathid, id)
    }
    path <- as.character(rnames[pathid])

    if(length(path) <= 1) {
      if(echo)
        cat("empty\n")
      vpvals <- c(vpvals, -1)
      next
    }
    
    if(echo)
      cat("[", length(vpvals)+1, "\\", n,"] ", ls[1]," (", length(path),") ... ")
    vpathways <- c(vpathways, pathname)
    
    nodeOrder <- 1:length(path)
    nodeCats <- lapply(1:length(path), function(i) return(1:numCategories))
    names(nodeCats) <- path

    path.data <- data[pathid, ]
    
    catnet.res <- cnSearchOrder(path.data, NULL, numParents, NULL, 0, 
                         nodeOrder=nodeOrder, nodeCats=nodeCats, 
                         parentsPool=NULL, NULL, echo=FALSE)

    bnet <- catnet.res@nets[[1]]
    if(model == "BIC") 
      bnet <- cnFindBIC(catnet.res)
    else if(model == "AIC") 
      bnet <- cnFindAIC(catnet.res)
    net1 <- cnSetProb(bnet, path.data[,classes==1], nodeCats=nodeCats) 
    net2 <- cnSetProb(bnet, path.data[,classes==2], nodeCats=nodeCats)    
    l1 <- net1@likelihood
    l2 <- net2@likelihood
    pval <- abs(l1-l2)/length(path)
      
    vpvals <- c(vpvals, pval)

    if(echo)
      if(2*pval > alpha && pval < alpha)
        cat(pval, " *\n")
      else if(pval > alpha)
        cat(pval, " **\n")
      else
        cat(pval, "\n")

    if(!is.null(outDir) && pval > alpha) {
      save(bnet, file=paste(outDir, "/", pathname, "_catnet.Rda", sep=""))
    }
  }
  names(vpvals) <- vpathways
  
  return(vpvals)
}

cnPathwayTest2 <- function(pathwayFile, cdata, classes, numParents=1, numCategories=3, complexity="BIC", bootFact=10, echo=TRUE, alpha=0.05, out.dir=NULL,...) {
  if(is.null(pathwayFile))
    stop("Valid pathway file is expected")
  lines <- readLines(pathwayFile)
  n <- length(lines)
  if(n <= 0)
    stop(pathwayFile, " is invalid")
  
  if(!is.matrix(cdata) && !is.data.frame(cdata))
    stop("data should be a matrix or data frame")
  if(is.data.frame(cdata))
    cdata <- t(cdata)
  if(!is.numeric(cdata))
    stop("data should be numeric")
  numnodes <- nrow(cdata)
  numsamples <- ncol(cdata)
  if(numnodes < 1 || numsamples < 1)
    stop("Invalid sample")

  lvs <- levels(classes)
  if(!is.factor(classes) || length(lvs) != 2)
    stop("Two-classes factor vector is expected")
  if(length(classes) != numsamples)
    stop("One class value for each sample instance is expected")
  classes[classes==lvs[1]] <- 1
  classes[classes==lvs[2]] <- 2
  classes <- as.integer(classes)
  numsamples1 <- sum(classes==1)
  numsamples2 <- sum(classes==2)

  if(complexity != "BIC" && complexity != "AIC" && !is.integer(complexity))
    stop("Invalid complexity")
  
  data <- cnDiscretize(cdata, numCategories, "uniform")

  vpathways <- NULL
  vpvals <- NULL
  for(ll in lines) {
    
    ls <- strsplit(ll,"\t")[[1]]
    if(length(ls) < 3) {
      warning("no valid pathway ", ls)
      next
    }
    pathname <- ls[1]
    path <- ls[3:length(ls)]

    pathid <- NULL
    rnames <- rownames(data)
    for(i in 1:length(path)) {
      id <- which(rnames == path[i])
      if(length(id) != 1)
        next
      str <- rnames[id]
      if(nchar(str) > 16 || length(which(pathid==id[1]))>0)
        next
      pathid <- c(pathid, id)
    }
    path <- as.character(rnames[pathid])

    if(length(path) < 1) {
      if(echo)
        cat("empty\n")
      next
    }
    
    if(echo)
      cat("[", length(vpvals)+1, "\\", n,"] ",ls[1]," ... ")
    vpathways <- c(vpathways, pathname)
    
    nodeOrder <- 1:length(path)
    nodeCats <- lapply(1:length(path), function(i) return(1:numCategories))
    names(nodeCats) <- path

    path.data <- data[pathid, ]
    
    catnet.res <- cnSearchOrder(path.data, NULL, numParents, NULL, 0, 
                         nodeOrder=nodeOrder, nodeCats=nodeCats, 
                         parentsPool=NULL, NULL, echo=FALSE)
    if(is.integer(complexity))
      bnet <- cnFind(catnet.res, complexity)
    if(complexity=="BIC")
      bnet <- cnFindBIC(catnet.res)
    if(complexity=="AIC")
      bnet <- cnFindAIC(catnet.res)

    bnet <- catnet.res@nets[[1]]
    ldiff <- 0
    for(bn in catnet.res@nets) {
      net1 <- cnSetProb(bn, path.data[,classes==1], nodeCats=nodeCats) 
      net2 <- cnSetProb(bn, path.data[,classes==2], nodeCats=nodeCats)
      diff <- abs(net1@likelihood - net2@likelihood)
      dd <- c(dd,diff)
      if(ldiff < diff) {
        ldiff <- diff
        bnet <- bn
      }
    }

    lhisto <- NULL
    for(bi in 1:(bootFact*numsamples)) {
      ind <- sample(1:numsamples, numsamples, replace=FALSE)
      net1 <- cnSetProb(bnet, path.data[,classes[ind]==1], nodeCats=nodeCats) 
      net2 <- cnSetProb(bnet, path.data[,classes[ind]==2], nodeCats=nodeCats)
      lhisto <- c(lhisto, net1@likelihood - net2@likelihood)
    }
    pval <- min(sum(lhisto<ldiff)/length(lhisto), sum(lhisto>ldiff)/length(lhisto))
    vpvals <- c(vpvals, pval)

    if(echo)
      if(pval < alpha)
         cat(pval, " *\n")
         else
         cat(pval, "\n")
      
    if(!is.null(out.dir) && pval < alpha) {
      save(bnet, file=paste(out.dir, "/", pathname, "_catnet.Rda", sep=""))
    }
    
  }
  names(vpvals) <- vpathways
  
  return(vpvals)
}

