
cnCatnetFromXdsl <- function(file) {

  lines <- tryCatch(readLines(file), error=function(e) NA, warning=function(e) NA)
  nlines <- length(lines)
  if(nlines<2) {
    stop("cannot read ", file)
    return(NULL)
  }

  nl <- 0
  word <- ""
  while(nl < nlines) {
    nl <- nl+1
    str <- unlist(strsplit(lines[nl], split=" "))
    if(length(str) < 1)
      next
    while(str[1]=="") str <- str[-1]
    word <- str[1]
    while(substr(word,1,1)=="\t") word <- substr(word,2,nchar(word))
    if(substr(word,1,1)=="<") word <- substr(word,2,nchar(word))
    if(word == "smile")
      break
  }
  if(word != "smile")
    stop("wrong format: ", word)
  word <- ""
  while(nl < nlines) {
    nl <- nl+1
    str <- unlist(strsplit(lines[nl], split=" "))
    if(length(str) < 1)
      next
    while(str[1]=="") str <- str[-1]
    word <- str[1]
    while(substr(word,1,1)=="\t") word <- substr(word,2,nchar(word))
    if(substr(word,1,1)=="<") word <- substr(word,2,nchar(word))
    if(substr(word,nchar(word),nchar(word))==">") word <- substr(word,1,nchar(word)-1)
    if(word == "nodes")
      break
  }
  if(word != "nodes")
    stop("wrong format: ", word)

  nnodes <- 0
  nodes <- NULL
  lpars <- NULL
  lcats <- NULL
  lprobs <- NULL
  word <- ""
  while(nl < nlines) {
    nl <- nl + 1
    str <- unlist(strsplit(lines[nl], split=" "))
    while(str[1]=="") str <- str[-1]
    word <- str[1]
    while(substr(word,1,1)=="\t") word <- substr(word,2,nchar(word))
    if(substr(word,1,1)=="<") word <- substr(word,2,nchar(word))
    if(substr(word,nchar(word),nchar(word))==">") word <- substr(word,1,nchar(word)-1)
    if(word == "/nodes")
      break
    if(word == "cpt") {
      word <- str[2]
      if(substr(word,nchar(word),nchar(word))==">") word <- substr(word,1,nchar(word)-1)
      if(substr(word,1,3) != "id=")
        next
      word <- substr(word,4,nchar(word))
      if(substr(word,1,1)=="\"") word <- substr(word,2,nchar(word))
      if(substr(word,nchar(word),nchar(word))=="\"") word <- substr(word,1,nchar(word)-1)
      node <- word
      cats <- NULL
      pars <- NULL
      probs <- NULL
      while(nl < nlines) {
        nl <- nl + 1
        str <- unlist(strsplit(lines[nl], split=" "))
        while(str[1]=="") str <- str[-1]
        word <- str[1]
        while(substr(word,1,1)=="\t") word <- substr(word,2,nchar(word))
        if(substr(word,1,1)=="<") word <- substr(word,2,nchar(word))
        if(word == "/cpt>")
          break

        if(word == "state") {
          word <- str[2]
          if(substr(word,nchar(word),nchar(word))==">") word <- substr(word,1,nchar(word)-1)
          if(substr(word,1,3) != "id=")
            next
          word <- substr(word,4,nchar(word))
          if(substr(word,1,1)=="\"") word <- substr(word,2,nchar(word))
          if(substr(word,nchar(word),nchar(word))=="\"") word <- substr(word,1,nchar(word)-1)
          cats <- c(cats, word)
        }

        if(substr(word,1,8) == "parents>") {
          word <- substr(word,9,nchar(word))
          if(substr(word,nchar(word)-9,nchar(word))=="</parents>") word <- substr(word,1,nchar(word)-10)
          pars <- word
          while(length(str)>1) {
            str <- str[-1]
            word <- str[1]
            if(substr(word,nchar(word)-9,nchar(word))=="</parents>") word <- substr(word,1,nchar(word)-10)
            pars <- c(pars,word)
          }
        }
        
        if(substr(word,1,14) == "probabilities>") {
          word <- substr(word,15,nchar(word))
          if(substr(word,nchar(word)-15,nchar(word))=="</probabilities>") word <- substr(word,1,nchar(word)-16)
          probs <- word
          while(length(str)>1) {
            str <- str[-1]
            word <- str[1]
            if(substr(word,nchar(word)-15,nchar(word))=="</probabilities>") word <- substr(word,1,nchar(word)-16)
            probs <- c(probs, word)
          }
        }
      }

      #at("node = ", node, "\n")
      #at(cats,"\n")
      #at(pars,"\n")
      #at(probs,"\n")

      nodes <- c(nodes, node)
      if(is.null(lpars))
        lpars <- list("")
      else
        lpars <- c(lpars, "")
      if(!is.null(pars))
        lpars[[length(lpars)]] <- pars
      if(is.null(lcats))
        lcats <- list("")
      else
        lcats <- c(lcats, "")
      lcats[[length(lcats)]] <- cats
      if(is.null(lprobs))
        lprobs <- list("")
      else
        lprobs <- c(lprobs, "")
      lprobs[[length(lprobs)]] <- probs
      
      nnodes <- nnodes + 1
    } ## cpt
  }

  word <- ""
  while(nl < nlines) {
    nl <- nl+1
    str <- unlist(strsplit(lines[nl], split=" "))
    if(length(str) < 1)
      next
    while(str[1]=="") str <- str[-1]
    word <- str[1]
    while(substr(word,1,1)=="\t") word <- substr(word,2,nchar(word))
    if(substr(word,1,1)=="<") word <- substr(word,2,nchar(word))
    if(substr(word,nchar(word),nchar(word))==">") word <- substr(word,1,nchar(word)-1)
    if(word == "extensions")
      break
  }
  if(word != "extensions")
    stop("wrong format, expected extensions: ", word)
  word <- ""
  while(nl < nlines) {
    nl <- nl+1
    str <- unlist(strsplit(lines[nl], split=" "))
    if(length(str) < 1)
      next
    while(str[1]=="") str <- str[-1]
    word <- str[1]
    while(substr(word,1,1)=="\t") word <- substr(word,2,nchar(word))
    if(substr(word,1,1)=="<") word <- substr(word,2,nchar(word))
    if(substr(word,nchar(word),nchar(word))==">") word <- substr(word,1,nchar(word)-1)
    if(word == "genie")
      break
  }
  if(word != "genie")
    stop("wrong format, expected genie: ", word)
  netname <- ""
  str <- str[-c(1,2)]
  word <- str[1]
  if(substr(word, 1,5)=="name=") {
    word <- substr(word, 6, nchar(word))
    if(substr(word,1,1)=="\"") word <- substr(word,2,nchar(word))
    if(substr(word,nchar(word),nchar(word))=="\"") {
      word <- substr(word,1,(nchar(word)-1))
      netname <- word
    }
    else {
      netname <- word
      while(length(str)>1) {
        str <- str[-1]
        word <- str[1]
        if(substr(word,nchar(word),nchar(word))=="\"") {
          word <- substr(word,1,(nchar(word)-1))
          netname <- paste(netname, word)
          break
        }
        netname <- paste(netname, word)
      }
    }
  }

  lnodes <- nodes
  while(nl < nlines) {
    nl <- nl + 1
    str <- unlist(strsplit(lines[nl], split=" "))
    while(str[1]=="") str <- str[-1]
    word <- str[1]
    while(substr(word,1,1)=="\t") word <- substr(word,2,nchar(word))
    if(substr(word,1,1)=="<") word <- substr(word,2,nchar(word))
    if(substr(word,nchar(word),nchar(word))==">") word <- substr(word,1,nchar(word)-1)
    if(word == "/genie")
      break
    if(word == "node") {
      word <- str[2]
      if(substr(word,nchar(word),nchar(word))==">") word <- substr(word,1,nchar(word)-1)
      if(substr(word,1,3) != "id=")
        next
      word <- substr(word,4,nchar(word))
      if(substr(word,1,1)=="\"") word <- substr(word,2,nchar(word))
      if(substr(word,nchar(word),nchar(word))=="\"") word <- substr(word,1,nchar(word)-1)
      id <- which(nodes==word)
      while(nl < nlines) {
        nl <- nl + 1
        str <- unlist(strsplit(lines[nl], split=" "))
        while(str[1]=="") str <- str[-1]
        word <- str[1]
        while(substr(word,1,1)=="\t") word <- substr(word,2,nchar(word))
        if(substr(word, 1, 6) == "<name>") {
          word <- substr(word, 7, nchar(word))
          if(substr(word,nchar(word)-6,nchar(word))=="</name>") {
            word <- substr(word, 1, nchar(word)-7)
            lnodes[id] <- word
          }
          else {
            lnodes[id] <- word
            while(length(str)>1) {
              str <- str[-1]
              word <- str[1]
              if(substr(word,nchar(word)-6,nchar(word))=="</name>")
                word <- substr(word, 1, nchar(word)-7)
              lnodes[id] <- paste(lnodes[id], word)
            }
          }
        }
        if(substr(word, 1, 7) == "</node>")
          break
      }
    }
  }
  
  set.probs <- function(idroot, ppars, pcatlist, idx, probtable, lentable) {
    ncats <- length(pcatlist[[idroot]])
    if(is.null(ppars) || length(ppars) < 1 || length(idx) < 1) {
      if(ncats < 1 || lentable < 1)
        return(NULL)
      return(as.vector(probtable[1:lentable]))
    }
    ncats <- length(pcatlist[[ppars[idx[1]]]])
    lentable <- as.integer(lentable/ncats)
    poutlist <- lapply(0:(ncats-1),
                       function(cat)
                       set.probs(idroot, ppars, pcatlist, idx[-1],
                                 probtable[(1+cat*lentable):((cat+1)*lentable)], lentable))
    return(poutlist)
  }

  lnpars <- vector("list", nnodes)
  for( i in 1:nnodes)
    if(length(lpars[[i]])>0 && lpars[[i]] != "")
      lnpars[[i]] <- sapply(length(lpars[[i]]):1, function(ip) which(nodes==lpars[[i]][ip]))

    lnprobs <- lapply(1:nnodes, function(i)
                      set.probs(i,lnpars[[i]],
                                lcats,1:length(lpars[[i]]),as.numeric(lprobs[[i]]),length(lprobs[[i]])))
  
  cn <- cnNew(lnodes, lcats, lnpars, lnprobs)
  cn@meta <- netname
  
  return(cn) 
} 
