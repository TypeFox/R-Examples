
load.gsea <- function(pathwatFile="c2.cp.kegg.v3.0.symbols.gmt", path, npars=2, ncats=3) {

  lines <- readLines(pathwayFile)
  n <- length(lines)
  if(is.character(path)) { 
    for(ll in lines) {
      ls <- strsplit(ll,"\t")[[1]]
      if(ls[1] == path)
        break
    }
    if(length(ls) < 1 || ls[1] != path)
      stop("no valid pathway")
  }
  else if(is.numeric(path)) {
    path <- as.integer(path)
    if(path < 1 || path > n)
      path <- n
    ls <- strsplit(lines[[path]],"\t")[[1]]
    cat(ls[1],"\n")
  } 
 if(length(ls) < 3)
    stop("no valid pathway")
  pathname <- ls[1]
  path <- ls[3:length(ls)]
    
  data <- read.table("Diabetes_collapsed_symbols.gct", header=TRUE, sep='\t', skip=2)

  pathid <- NULL
  levs <- levels(data[,1])
  for(i in 1:length(path)) {
    id <- which(data[,1] == path[i])
    if(length(id) != 1)
      next
    str <- levs[data[id[1], 1]]
    if(nchar(str) > 16 || length(which(pathid==id[1]))>0)
      next
    pathid <- c(pathid, id[1])
  }
  path <- as.character(data[pathid,1])

  path <- as.character(sapply(path, function(ss) sub("-",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("-",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("@",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("@",".",ss)))
  path <- as.character(sapply(path, function(ss) sub(" /// ",".",ss)))
  path <- as.character(sapply(path, function(ss) sub(" /// ",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("_",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("_",".",ss)))
 
  cdata <- data[pathid, 3:ncol(data)]
  rownames(cdata) <- path
  numnodes <- nrow(cdata)
  numsamples <- ncol(cdata)
  cdata <- as.matrix(cdata, nrow=numnodes)
    
  ## reorder according to samples
  cls <- c(rep(1,17),rep(2,17))
 
  cat(numnodes, " nodes and ", numsamples, " samples\n")
  
  nodePars <- rep(npars, numnodes)  
  nodeOrder <- 1:numnodes
  nodeCats <- lapply(1:nrow(cdata), function(i) return(1:ncats))
  names(nodeCats) <- path

  return(list(pathname=pathname, cdata=cdata, cls=cls, pathway=path, nodeOrder=nodeOrder, nodePars=nodePars, nodeCats=nodeCats))
}

load.topK <- function(K=100, npars=3, ncats=3) {

  data <- read.table("Diabetes_collapsed_symbols.gct", header=TRUE, sep='\t', skip=2)

  cls <- c(rep(1,17),rep(2,17))
  
  cdata <- as.matrix(data[, 3:ncol(data)])
  pvals <-  rep(0, nrow(cdata))
  ##ksvals <-  rep(0, nrow(cdata))
  for(j in 1:nrow(cdata)) {
    tt <- t.test(cdata[j,cls==1&!is.na(cdata[j,])], cdata[j,cls==2&!is.na(cdata[j,])])
    pvals[j] <- tt$p.value
    ##ksvals[j] <- ks.test((cdata[j,]-mean(cdata[j,]))/sqrt(var(cdata[j,])), pnorm)$p.value
  }
  qq <- quantile(pvals,(K)/nrow(cdata))
  pathid <- which(pvals<qq)
  path <- as.character(data[pathid,1])

  path <- as.character(sapply(path, function(ss) sub("-",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("-",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("@",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("@",".",ss)))
  path <- as.character(sapply(path, function(ss) sub(" /// ",".",ss)))
  path <- as.character(sapply(path, function(ss) sub(" /// ",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("_",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("_",".",ss)))
  
  cdata <- as.matrix(data[pathid, 3:ncol(data)])
  rownames(cdata) <- path
  numnodes <- nrow(cdata)
  numsamples <- ncol(cdata)

  nodePars <- rep(npars, numnodes)
  
  nodeOrder <- rank(pvals[pathid])

  nodeCats <- lapply(1:nrow(cdata), function(i) return(1:ncats))
  names(nodeCats) <- path
  
  return(list(cdata=cdata, cls=cls, pathway=path, nodeOrder=nodeOrder, nodePars=nodePars, nodeCats=nodeCats))
}


