#  File R/EgoStat.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
# An EgoStat.* function takes the data frame of egos, a data frame of
# alters, and the arguments passed to the corresponding ERGM terms,
# and returns a matrix of h(e[i]) values, with egos in rows and
# elements of h(e[i]) in columns.

EgoStat.edges <- function(egodata){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
   

  ties<-merge(egos[egoIDcol],alters[egoIDcol],by=egoIDcol,suffixes=c(".ego",".alter"))

  alterct <- as.data.frame(table(ties[[egoIDcol]]), stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[egoIDcol],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  h <- cbind(egos$.degree)
  colnames(h) <- "edges"
  rownames(h) <- egos[[egoIDcol]]

  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]/2
}


EgoStat.nodecov <- function(egodata, attrname){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
   
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties <- data.frame(egoID=c(ties[[egoIDcol]],ties[[egoIDcol]],isolates),x=c(ties$.e,ties$.a,rep(0,length(isolates))),stringsAsFactors=FALSE)
  
  h <- cbind(sapply(tapply(ties$x,list(egoID=ties$egoID),FUN=sum),identity))
  colnames(h) <- paste("nodecov",attrname,sep=".")
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]/2
}


EgoStat.nodefactor <- function(egodata, attrname, base=1){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties <- data.frame(egoID=c(ties[[egoIDcol]],ties[[egoIDcol]],isolates),x=c(ties$.e,ties$.a,rep(0,length(isolates))),stringsAsFactors=FALSE)
  
  h <- t(sapply(tapply(ties$x, list(egoID=ties$egoID), FUN=tabulate, nbins=length(levs)),identity))
  colnames(h) <- paste("nodefactor",attrname,levs,sep=".")  

  if(length(base)==0 || base==0) h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]/2
  else h[match(egodata$egos[[egoIDcol]],rownames(h)),-base,drop=FALSE]/2
}

EgoStat.nodematch <- function(egodata, attrname, diff=FALSE, keep=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  ties$match <- ifelse(ties$.e==ties$.a, as.integer(ties$.e), 0)
  if(!is.null(keep)) ties$match[!(ties$match%in%keep)] <- 0
  if(!diff) ties$match[ties$match!=0] <- 1
  
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 

  ties <- data.frame(egoID=c(ties[[egoIDcol]],isolates), match=c(ties$match,rep(0,length(isolates))),stringsAsFactors=FALSE)

  h <- t(rbind(sapply(tapply(ties$match, list(egoID=ties$egoID), FUN=tabulate, nbins=if(diff) length(levs) else 1),identity)))

  colnames(h) <- if(diff) paste("nodematch",attrname,levs,sep=".") else paste("nodematch",attrname,sep=".")
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),if(!is.null(keep) && diff) keep else TRUE,drop=FALSE]/2
}


EgoStat.nodemix <- function(egodata, attrname, base=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  
  
  namevec <- outer(levs,levs,paste,sep=".")
  namevec <- namevec[upper.tri(namevec,diag=TRUE)]
  
  if (!is.null(base) && !identical(base,0)) {
    namevec <- namevec[-base]
  }
  
  mat <- t(apply(cbind(ties[,2:3]),1,sort))
  h <- table(ties[,1],paste(levs[mat[,1]],levs[mat[,2]],sep="."))
  
  h<- t(apply(h,1,function(x)x[namevec,drop=FALSE]))
  
  if(length(isolates)){
    isolates.mat <- matrix(0,nrow=length(isolates),ncol=length(namevec))
    rownames(isolates.mat) <- isolates	
    h <- rbind(h,isolates.mat)
  }
  
  h <- h[order(as.numeric(rownames(h))),]
  h[is.na(h)] <- 0
  colnames(h) <- paste("mix",attrname,namevec,sep=".")
  h[match(egodata$egos[[egoIDcol]],rownames(h)),]/2
}

EgoStat.absdiff <- function(egodata, attrname, pow=1){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties)<-c(egoIDcol,".e",".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties$.absdiff <- abs(ties$.e-ties$.a)^pow 
  
  ties <- data.frame(egoID=c(ties[[egoIDcol]],isolates), absdiff=c(ties$.absdiff,rep(0,length(isolates))),stringsAsFactors=FALSE)
  
  h <- cbind(sapply(tapply(ties$absdiff,list(egoID=ties$egoID),FUN=sum),identity))         
  colnames(h) <- if(pow==1) paste("absdiff",attrname,sep=".") else paste("absdiff",pow,".",attrname,sep="")
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]/2
}

EgoStat.degree <- function(egodata, d, by=NULL, homophily=FALSE){
  ## if(any(d==0)) warning("degree(0) (isolate) count statistic depends strongly on the specified population network size.")
  
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  if(!is.null(by)){
    levs <- sort(unique(c(egos[[by]],alters[[by]])))
  }
  
  ties<-merge(egos[c(egoIDcol,by)],alters[c(egoIDcol,by)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(by)) names(ties) <- c(egoIDcol,".e",".a")
  if(!is.null(by) && homophily) ties <- ties[ties$.e==ties$.a,]
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol,by)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  if(!is.null(by) && !homophily){
    bys <- rep(levs,each=length(d))
    degs <- rep(d,length(levs))
    
    h <- sapply(seq_along(bys), function(i) egos$.degree==degs[i] & egos[[by]]==bys[i])
    colnames(h) <- paste("deg",degs,".",by,bys,sep="")
  }else{
    h <- sapply(d, function(i) egos$.degree==i)
    colnames(h) <- if(homophily) paste("deg",d,".homophily.",by,sep="") else paste("degree",d,sep="")
  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}

EgoStat.degrange <- function(egodata, from=NULL, to=Inf, byarg=NULL, homophily=FALSE){
  ## if(any(from==0)) warning("degrange(0,...) (isolate) count depends strongly on the specified population network size.")
  
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  to <- ifelse(to==Inf, .Machine$integer.max, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term degrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange must have from<to.")

  if(!is.null(byarg)){
    levs <- sort(unique(c(egos[[byarg]],alters[[byarg]])))
  }

  ties<-merge(egos[c(egoIDcol,byarg)],alters[c(egoIDcol,byarg)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(byarg)) names(ties) <- c(egoIDcol,".e",".a")
  if(!is.null(byarg) && homophily) ties <- ties[ties$.e==ties$.a,]
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol,byarg)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  if(!is.null(byarg) && !homophily){
    bys <- rep(levs,each=length(from))
    froms <- rep(from,length(levs))
    tos <- rep(to,length(levs))
    
    h <- sapply(seq_along(bys), function(i) egos$.degree>=froms[i] & egos$.degree<tos[i] & egos[[byarg]]==bys[i])
    colnames(h) <-  ifelse(tos>=.Machine$integer.max,
                           paste("deg", from, "+.",          byarg, bys, sep=""),
                           paste("deg", from, "to", to, ".", byarg, bys, sep=""))

  }else{
    h <- sapply(seq_along(from), function(i) egos$.degree>=from[i] & egos$.degree<to[i])
    colnames(h) <-
      if(homophily)
        ifelse(to>=.Machine$integer.max,
               paste("deg", from,  "+",     ".homophily.", byarg, sep=""),
               paste("deg", from, "to", to, ".homophily.", byarg, sep=""))
      else
        ifelse(to>=.Machine$integer.max,
               paste("deg", from,  "+", sep=""),
               paste("deg", from, "to", to, sep=""))

  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}

EgoStat.concurrent <- function(egodata, by=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  if(!is.null(by)){
    levs <- sort(unique(c(egos[[by]],alters[[by]])))
  }

  ties<-merge(egos[c(egoIDcol,by)],alters[c(egoIDcol,by)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(by)) names(ties) <- c(egoIDcol,".e",".a")
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol,by)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  if(!is.null(by)){
    bys <- levs
    
    h <- sapply(seq_along(bys), function(i) egos$.degree>=2 & egos[[by]]==bys[i])
    colnames(h) <- paste("concurrent.", by, bys, sep="")

  }else{
    h <- cbind(egos$.degree>=2)
    colnames(h) <- "concurrent"
  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}

EgoStat.concurrentties <- function(egodata, by=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  if(!is.null(by)){
    levs <- sort(unique(c(egos[[by]],alters[[by]])))
  }

  ties<-merge(egos[c(egoIDcol,by)],alters[c(egoIDcol,by)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(by)) names(ties) <- c(egoIDcol,".e",".a")
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol,by)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  if(!is.null(by)){
    bys <- levs
    
    h <- sapply(seq_along(bys), function(i) cbind(ifelse(egos[[by]]==bys[i], pmax(egos$.degree-1,0), 0)))
    colnames(h) <- paste("concurrentties.", by, bys, sep="")

  }else{
    h <- cbind(pmax(egos$.degree-1,0))
    colnames(h) <- "concurrentties"    
  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}


EgoStat.degreepopularity <- function(egodata){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
    
  ties<-merge(egos[egoIDcol],alters[egoIDcol],by=egoIDcol,suffixes=c(".ego",".alter"))

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  h <- cbind(egos$.degree^(3/2))
  colnames(h) <- "degreepopularity"
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}
