## BIF: Interchange Format for Bayesian Networks

bif2char <- function(bstr) {
  nc <- nchar(bstr)
  if(substr(bstr,1,1)=="\"")
    bstr <- substr(bstr,2,nc)
  nc <- nc-1
  if(substr(bstr,nc,nc)=="\"")
    bstr <- substr(bstr,1,nc-1)
  if(substr(bstr,nc-1,nc-1)=="\"")
    bstr <- substr(bstr,1,nc-2)
  return(bstr)
}

bif.parse.var <- function(strvar) {
  nodename <- NA
  strvar <- unlist(strsplit(strvar, split=" "))
  nc <- which(strvar=="")
  if(length(nc)>0)
    strvar <- strvar[-nc]
  if(strvar[1]!="variable")
    return(NULL)
  strvar <- strvar[-1]
  tagnode <- bif2char(strvar[1])
  strvar <- strvar[-1]
  if(strvar[1] != "{" || strvar[length(strvar)] != "}")
    return(NULL)
  strvar <- strvar[-c(1, length(strvar))]
  if(strvar[1] =="property" && strvar[2] == "label" && strvar[3] == "=") {
    nodename <- bif2char(strvar[4])
    strvar <- strvar[-c(1,2,3,4)]
    if(length(strvar) < 1)
      return(NULL)
  }
  while(length(strvar) > 0 && strvar[1] != "type")
    strvar <- strvar[-1]
  if(substr(strvar[2],1,8)  != "discrete") 
      return(NULL)
  while(strvar[1] != "{")
    strvar <- strvar[-1]
  if(strvar[1] != "{")
    return(NULL)
  strvar <- strvar[-1]
  nc <- length(strvar)
  if(nc<2)
    return(NULL)
  while(strvar[nc] != "}" && strvar[nc] != "};")
    nc <- nc - 1
  if(nc < 2)
    return(NULL)
  strvar <- strvar[1:(nc-1)]
  for(i in 1:length(strvar)) {
    nc <- nchar(strvar[i])
    if(substr(strvar[i], nc,nc) == ",")
      strvar[i] <- substr(strvar[i], 1,nc-1)
  }
  if(is.na(nodename))
    nodename <- tagnode
  return(list(tagnode, nodename, strvar))
}

cnCatnetFromBif <- function(file) {

  lines <- tryCatch(readLines(file), error=function(e) NA, warning=function(e) NA)
  nlines <- length(lines)
  if(nlines<2) {
    stop("cannot read ", file)
    return(NULL)
  }

  nl <- 0
  while(nl < nlines) {
    nl <- nl+1
    str <- unlist(strsplit(lines[nl], split=" "))
    if(length(str) < 1)
      next
    while(str[1]=="") str <- str[-1]
    if(str[1] == "network")
      break
  }
  if(str[1] != "network")
    stop("wrong format")
  netname <- bif2char(str[2])
  while(nl < nlines) {
    str <- NULL
    while(length(str) < 1 && nl < nlines) {
      nl <- nl+1
      str <- unlist(strsplit(lines[nl], split=" "))
    }
    while(str[1]=="") str <- str[-1]
    if(str[1] == "}")
      break
  }
  
  nodes <- NULL
  cats <- NULL
  nnodes <- 0
  tagnodes <- NULL
  while(nl < nlines) {
    str <- NULL
    while(length(str) < 1 && nl < nlines) {
      nl <- nl+1
      str <- unlist(strsplit(lines[nl], split=" "))
    }
    while(str[1]=="") str <- str[-1]
    if(str[1] != "variable" && str[1] != "probability") 
      next
    if(str[1] == "probability")
      break
    if(length(which(str == "{")) < 1)
      stop("wrong format, line ", nl)
    nstrend <- which(str == "//")
    if(length(nstrend) == 1)
      str <- str[1:(nstrend-1)]
    strvar <- paste(str, collapse=" ")
    while(1) {
      str <- NULL
      while(length(str) < 1 && nl < nlines) {
        nl <- nl+1
        str <- unlist(strsplit(lines[nl], split=" "))
      }
      while(str[1]=="") str <- str[-1]
      nstrend <- which(str == "//")
      if(length(nstrend) == 1)
        str <- str[1:(nstrend-1)]
      str <- paste(str, collapse=" ")
      for(i in 1:(nchar(str)-1)){
        if(substr(str, i, i+1) == "//") {
          str <- substr(str, 1, i-1)
          break
        }
      } 
      strvar <- paste(strvar, str, collapse=" ")
      if(length(which(str == "}")) == 1)
        break
    }
    node.var <- bif.parse.var(strvar)
    if(is.null(node.var))
      stop("wrong variable, line ", nl)

    tagnodes <- c(tagnodes, node.var[[1]])
    nodes <- c(nodes, node.var[[2]])
    nnodes <- nnodes+1 
    newcats <- vector("list", nnodes)
    if(nnodes > 1)
      for(i in 1:(nnodes-1))
        newcats[[i]] <- cats[[i]]
    newcats[[nnodes]] <- node.var[[3]]
    cats <- newcats
  }

  lpars <- vector("list", nnodes)
  lprobs <- vector("list", nnodes)
  while(nl < nlines) {
    if(str[1] != "probability") {
      str <- NULL
      while(length(str) < 1 && nl < nlines) {
        nl <- nl+1
        str <- unlist(strsplit(lines[nl], split=" "))
      }
      while(str[1]=="") str <- str[-1]
      next
    }
    strnode <- str[3]
    nnode <- which(tagnodes==strnode)
    if(length(nnode)!=1)
      stop("wrong node ", strnode)
    c1 <- which(str == "|")
    c2 <- which(str == ")")
    if(length(c2)==0)
      c2 <- which(str == "){")
    if(length(c1)!=1)
      c1 <- c2
    nodepars <- NULL
    if(c1 < c2-1) {
      for(cc in (c1+1):(c2-1)) {
        spar <- str[cc]
        if(substr(spar, nchar(spar),nchar(spar)) == ",")
          spar <- substr(spar, 1,nchar(spar)-1)
        if(substr(spar, 1,1) == ",")
          spar <- substr(spar, 2,nchar(spar))
        npar <- which(tagnodes==spar)
        if(length(npar)==1)
          nodepars <- c(nodepars, npar)
      }
      lpars[[nnode]] <- nodepars
    }

    nl <- nl + 1
    k <- 0
    while(nl+k < nlines) {
      strend <- unlist(strsplit(lines[nl+k], split=" "))
      while(strend[1]=="") strend <- strend[-1]
      if(strend[1] == "property") {
        nl <- nl+1
        next
      }
      if(strend[1] == "}" || strend[1] == "};")
        break
      k <- k+1
    }
    nodeprobs <- vector("list", k)
    k <- 1
    while(nl < nlines) {
      strend <- unlist(strsplit(lines[nl], split=" "))
      while(strend[1]=="") strend <- strend[-1]
      if(strend[1] == "}" || strend[1] == "};")
        break
      nodeprobs[[k]] <- unlist(strsplit(lines[nl], split=" "))
      nl <- nl+1
      k <- k+1
    }
    lprobs[[nnode]] <- nodeprobs
    nl <- nl+1
    if(nl >= nlines)
      break
    str <- unlist(strsplit(lines[nl], split=" "))
    while(str[1]=="") str <- str[-1]
  }
  
  cn <- cnNew(nodes, cats, lpars)
  cn@meta <- netname

  for(i in 1:nnodes) {
    nodecats <- cn@cats[[i]]
    nodepars <- cn@pars[[i]]
    if(length(nodepars)<1) {
      str <- lprobs[[i]][[1]]
      ns <- length(str)
      nns <- 1
      while(nns<ns) {
        if(str[nns]=="table")
          break
        nns <- nns+1
      }
      nns <- nns+1
      cprobs <- NULL
      while(nns <= ns) {
        snum <- str[nns]
        if(substr(snum, nchar(snum),nchar(snum)) == "," || substr(snum, nchar(snum),nchar(snum)) == ";")
          snum <- substr(snum, 1,nchar(snum)-1)
        cprobs <- c(cprobs, as.numeric(snum))
        nns <- nns+1
      }
      cn@probs[[i]] <- cprobs
      next
    }
    ## nodeparents > 0
    probtable <- vector("list", length(lprobs[[i]]))
    idtable <- 0
    for(j in 1:length(lprobs[[i]])) {
      str <- lprobs[[i]][[j]]
      while(str[1]=="") str <- str[-1]
      strcat <- paste(str[1:length(nodepars)], collapse="")
      ns <- nchar(strcat)
      if(substr(strcat, 1,1) != "(" || substr(strcat, ns,ns) != ")")
        stop("wrong parents ", strcat)
      strcat <- substr(strcat, 2, ns-1)
      strcat <- unlist(strsplit(strcat, split=","))
      parcats <- NULL
      for(np in 1:length(nodepars)) {
        cc <- strcat[np]
        nc <- which(cn@cats[[nodepars[np]]] == cc)
        if(length(nc) != 1)
          stop("wrong parent categories ", cc)
        parcats <- c(parcats, as.integer(nc))
      }
      nodeprobs <- NULL
      for(nns in 1:length(nodecats)) {
        while(str[length(nodepars)+nns]==",") str <- str[-c(length(nodepars)+nns)]
        strprob <- str[length(nodepars)+nns]
        ns <- nchar(strprob)
        if(substr(strprob, ns,ns) == "," || substr(strprob, ns,ns) == ";")
          strprob <- substr(strprob, 1,ns-1)
        nodeprobs <- c(nodeprobs, as.numeric(strprob))
      }
      idtable <- idtable+1
      probtable[[idtable]] <- nodeprobs
    }
    cn@probs[[i]] <- setNodeProbTable(i, nodepars, cn@cats, 1:length(nodepars), probtable, 1)
  }

  return(cn) 
} 

setNodeProbTable <- function(idroot, ppars, pcatlist, idx, probtable, idtable) {

  if(is.null(ppars) || length(ppars) < 1 || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 1) {
      return(NULL)
    }
    ##cat(idtable, "\n")
    return(as.vector(probtable[[idtable]]))
  }
  id <- ppars[idx[1]]
  off <- 1
  if(idx[1]>1)
    for(ip in 1:(idx[1]-1))
      off <- off*length(pcatlist[[ppars[ip]]])
  ##cat("pars ", nodepars[idx], ", ", off, "\n")
  poutlist <- lapply(1:length(pcatlist[[id]]),
                     function(cat)
                     setNodeProbTable(idroot, ppars, pcatlist, idx[-1], probtable, idtable+off*(cat-1)))
  return(poutlist)
}
