##
## Reading / writing Bayesian networks from / to HUGIN net files
##

loadHuginNet <- function(file, description=rev(unlist(strsplit(file, "/")))[1],
                         details=0){

  xxx      <-.readHuginNet(file,details)
  yyy      <-.transformHuginNet2internal(xxx)
  universe <- .asUniverse(yyy)
  plist    <- lapply(yyy$potentialList, .hpot2cptable, universe)
  value    <- grain(compileCPT(plist))
  return(value)
}


.transformHuginNet2internal <- function(x){
  nodeList2 <- lapply(x$nodeList, .getNodeSpec)
  potentialList2 <- lapply(x$potentialList, .getPotentialSpec)

  nl <- .makeNodeNamesUnique(nodeList2)

  repeat{
    if (length(nl$nonunique)==0)
      break()
    nl <- .makeNodeNamesUnique(nl$nodeList)
  }

  nodeList2 <- nl$nodeList

  value <- structure(list(nodeList=nodeList2, potentialList=potentialList2))
  class(value)<- "huginnet"

  return(value)
}


.readHuginNet <- function(file, details=0){

  .infoPrint(details, 1, cat(".HUGIN netfile:", file,"\n"))
  nodeCount <- 0
  con <- file(file, "rb")
  repeat{
    cline <- .getLine(con);  #print(cline)
    if (!length(cline))
      break()

    if (.hasToken("node", cline)) ## Fragile if 'node' is the name of a variable...
      nodeCount <- nodeCount + 1
  }
  close(con)

  .infoPrint(details, 3, cat("...there are around", nodeCount, "nodes \n"))

  ## Data structure for holding specification (possibly too long)
  ##
  nodeList <- potentialList <- as.list(rep(NA, nodeCount))

  con <- file(file, "rb")
  currNode <- currPotential <- 1
  state<-"start"
  repeat{
    cline <- .getLine(con);  #print(cline)
    if (!length(cline))
      break()
    switch(state,
           "start"={
             if (.hasToken("net",cline)){
               state="net"
               .infoPrint(details, 2, cat("..NET action\n"))
               wline <- cline
             }
           },
           "net"={
             wline <- c(wline, cline)
             if (.hasToken("}",cline)){
               state="run1"
               .infoPrint(details,2,cat("..end NET action\n"))
             }
           },
           "run1"={
             if (.hasToken("node", cline)){
               state="node"
               .infoPrint(details, 2, cat("..NODE action\n"))
             } else {
               if (.hasToken("potential", cline)){
                 state="potential";
                 .infoPrint(details,2, cat("..POTENTIAL action\n"))
               }
             }
             wline <- cline
           },
           "node"={
             wline <- c(wline, cline)
             if (.hasToken("}",cline)){
               state="run1";
               .infoPrint(details,2,cat("..end NODE action\n"))
               nodeList[[currNode]] <- wline;
               currNode <- currNode + 1
             }
           },
           "potential"={
             wline <- c(wline, cline)
             if (.hasToken("}",cline)){
               state="run1";
               .infoPrint(details,2, cat("..end POTENTIAL action\n"))
               potentialList[[currPotential]] <- wline;
               currPotential <- currPotential + 1
             }
           }
           )
  }
  close(con)

  nodeList <- nodeList[!sapply(lapply(nodeList, is.na),all)]
  potentialList <- potentialList[!sapply(lapply(potentialList, is.na),all)]


  value <- structure(list(nodeList=nodeList, potentialList=potentialList))
  return(value)
}



.asUniverse <- function(from){
  ccshort   <-sapply(from$nodeList, function(x)x$nodeVar)
  ccnames   <-sapply(from$nodeList, function(x)x$nodeLabel)
  cclabels  <-lapply(from$nodeList, function(x)x$nodeStates)
  names(cclabels) <- ccnames
  di <- c(lapply(cclabels, length),recursive=TRUE)
  list(nodes=ccnames, short=ccshort, levels=cclabels, nlev=di)
}


.hpot2cptable <- function(cpot, universe){
  idx <- match(c(cpot[c("nodeVar","parentVar")],recursive=TRUE), universe$short)
  vpa <- universe$nodes[idx]
  v   <- vpa[1]
  cptable(vpa, values=cpot$potential, levels=universe$levels[[v]])
}




.getLine   <- function(con) {
  readLines(con, n=1)
}

.hasToken  <- function(token, cline) {
  ##print(cline)
  cline <- gsub("^ +","",cline)
  a <- unlist(strsplit(cline," "))[1]

  if (!is.na(a))
    a==token
  else
    FALSE
}


.tokenIdx <- function(token, x){
  idx <- which(as.logical(lapply(x, function(d) grep(token,d))))
  idx
}


.capWords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
                           {s <- substring(s,2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

## .toCamel <- function(s){
##   s<-gsub(" +"," ",s)
##   s<-unlist(strsplit(s, " "))
##   paste(sapply(s, .capWords),collapse='')
## }

.toCamel <- function(s){
  s<-gsub(" +"," ",s)
  s<-unlist(strsplit(s, " "))
  paste(c(s[1],sapply(s[-1], .capWords)),collapse='')
}



.getNodeSpec <- function(nodeSpec){

  tmp <- nodeSpec[.tokenIdx("node", nodeSpec)]
  nodeVar <- gsub("node +","",tmp)[1]
  nodeVar <- gsub(" +","",nodeVar)


  tmp <- nodeSpec[.tokenIdx("label", nodeSpec)]
  nodeLabel <- gsub(" +label += +","",tmp);
  nodeLabel <- gsub(";", "", nodeLabel)
  nodeLabel <- gsub('"',"", nodeLabel)

  nodeLabel <- gsub(" +"," ",nodeLabel)

  if (length(nodeLabel) && nchar(nodeLabel)>0){

    nodeLabel <- .toCamel(nodeLabel)

    nl <- gsub("[^[:alnum:]]","",nodeLabel)
    nodeLabel <- gsub("[^[:alnum:]|\\.]","",nodeLabel)

    base<-as.character(0:9)
    if(subsetof(unlist(strsplit(nl,"")), base)){
      nodeLabel <- paste("X",nodeLabel,sep='')
    }
  } else {
    ##if (nchar(nodeLabel)==0)
    nodeLabel <- nodeVar
  }

  tmp <- nodeSpec[.tokenIdx("states", nodeSpec)]
  nodeStates <- gsub(" +states += +","",tmp);
  nodeStates <- gsub("[\\(,\\);]","",nodeStates);
  nodeStates <- unlist(strsplit(nodeStates, '\\"'))
  nodeStates <- sapply(nodeStates, function(d) gsub("^ +","",d))
  nodeStates <- nodeStates[sapply(nodeStates, nchar)>0]

  nodeStates <- sapply(nodeStates, .toCamel)

  nodeStates <- gsub(" +",".", nodeStates)
  names(nodeStates)<-NULL


  value <- list(nodeVar=nodeVar, nodeLabel=nodeLabel, nodeStates=nodeStates)
  value
}

.getPotentialSpec <- function(potSpec){
  tmp <- potSpec[.tokenIdx("potential", potSpec)]
  tmp <- gsub("potential +","", tmp)
  tmp <- gsub("[\\(,\\),|]","", tmp)
  tmp <- gsub(" +"," ", tmp)
  tmp <- unlist(strsplit(tmp," "))
  tmp <- tmp[sapply(tmp, nchar)>0]

  nodeVar <- tmp[1]
  parentVar <- tmp[-1]

  sss  <- paste(potSpec,collapse="") ##; ss <<- sss
  sss2 <- gsub("^.*data[[:space:]]*=([^;]*);(.*)", "\\1", sss) ##; ss2<<-sss2

  ##sss3: ((( 0.5 1.2E-5 ) ( 3E3 0.5 )) ( 0.5 0.5 ) ( 0.5 0.5 )))
  sss3 <- gsub("\\)[^\\)]*\\(", ") (", sss2) ##; ss3<<-sss3

  ## sss4: "  0.5 1.2E-5   3E3 0.5   0.5 0.5   0.5 0.5 "s
  sss4 <- gsub("[\\(,\\),\\}]","", sss3)

  ## sss5: remove leading white space: "0.5 1.2E-5   3E3 0.5   0.5 0.5   0.5 0.5 "
  sss5 <- gsub("^[[:space:]]*","",sss4)
  ## sss6: remove trailing white space: "0.5 1.2E-5   3E3 0.5   0.5 0.5   0.5 0.5"
  sss6 <- gsub("[[:space:]]$*","",sss5)
  ## sss7: split to atoms
  sss7 <- strsplit(sss6, " +")[[1]]

  ###: Now create numerical values
  pot <- as.numeric( sss7 )

  value <- list(nodeVar=nodeVar, parentVar=rev(parentVar), potential=pot)
  value
}




.makeNodeNamesUnique <- function(nodeList2){
  nl<-t(sapply(nodeList2, function(d)unlist(d[1:2])))

  nonunique <- names(which(table(nl[,2])>1))

  if (length(nonunique)){
    cat ("Label(s): {", nonunique, "} appears mode than once in NET file\n")
    for (i in 1:length(nonunique)){
      cnu <- nonunique[i]
      idx<-which(cnu ==nl[,2])
      for (j in idx){
        a <- nodeList2[[j]]$nodeVar
        cat("  Replacing label", cnu, " with node name", a, "\n")
        nodeList2[[j]]$nodeLabel <- a
      }
    }
  }

  return(list(nodeList=nodeList2, nonunique=nonunique))
}


saveHuginNet <- function(gin, file, details=0){

  cptlist <- gin$cptlist
  gmd     <- gin$universe

  vlab <- gmd$levels
  vnam <- gmd$nodes
  nn   <- length(vlab)

  th     <- cumsum(c(0,rep(2*pi/nn, nn-1)))
  r      <- 100
  coords <- lapply(th, function(d) round(r+r*c(cos(d), sin(d))))

  con <- file(file, "wb")

  ## Write (trivial) net specification
  ##
  writeLines("net\n{", con)
  writeLines("  node_size = (100 30);", con)
  writeLines("\n}\n\n", con)

  ## Write node specification
  ##
  for (ii in 1:length(vlab)){
    st<-paste("node ", vnam[ii],"\n","{","\n",sep='')
    writeLines(st,con,sep="")
    ## cat(st)
    st <- paste("   label = \"\";","\n")
    writeLines(st,con,sep="")
    ## cat(st)
    st <- paste("   position = (", paste(coords[[ii]], collapse=' '), ");\n")
    writeLines(st,con,sep="")
    ## cat(st)

    st2 <- sapply(vlab[[ii]], function(d) paste('"',d,'"',sep=''))
    st  <- paste("   states = (", paste(st2, collapse=' '), ");\n")
    writeLines(st,con,sep="")
    ## cat(st)
    st <- paste("}\n")
    writeLines(st,con,sep="")
    ## cat(st)
  }


  for (ii in 1:length(cptlist)){

    cpot <- cptlist[[ii]]
    nam <- varNames(cpot)    ## BRIS
    lev <- valueLabels(cpot) ## BRIS
    val <- cpot              ## BRIS

    v  <- nam[1]
    pa <- nam[-1]

    lev   <- rev(lev[-1])
    wval  <- val
    if (length(lev)>0){
      for (kk in 1:length(lev)){
        ##print("splitVec:"); print(wval); print(class(wval))
        wval<-splitVec(wval,length(lev[[kk]]))
      }
    }
    ##print(wval); print(class(wval))
    plx <- printlist(wval)

    if (length(pa)){
      st <- paste("potential (",v, "|", paste(rev(pa), collapse=' '),")\n")
      writeLines(st,con,sep="")
      ## cat(st)
      st <- "{\n";
      writeLines(st,con,sep="")
      ## cat(st)
      st <- paste("   data = \n")
      writeLines(st,con,sep="")
      ## cat(st)
      ##a<-lapply(plx, cat, "\n")
      a<-lapply(plx, writeLines, con, sep="\n")
      st <- paste(";\n}\n")
      writeLines(st,con,sep="")
      ## cat(st)

    } else {
      st <- paste("potential (", v, ")\n")
      writeLines(st,con,sep="")
      ## cat(st)
      st <- "{\n";
      writeLines(st,con,sep="")
      ## cat(st)
      st <- paste("   data = ", plx, ";\n")
      writeLines(st,con,sep="")
      ## cat(st)
      st <- "}\n\n";
      writeLines(st,con,sep="")
      ## cat(st)
    }
  }

  close(con)
}









































