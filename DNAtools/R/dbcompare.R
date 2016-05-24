dbCompare <- function(x, profiles = NULL, hit = 7, trace = TRUE, vector = FALSE, collapse = FALSE,
                      wildcard = FALSE, wildcard.effect = FALSE, wildcard.impose = FALSE,
                      Rallele = FALSE, threads = 2){
  
  if(all(wildcard,wildcard.effect))
    stop("Only one wildcard option is allowed to be 'TRUE' at any time")
  
  if(threads < 1){
    stop("Must have at least one thread")
  }
  
  x.cp <- x ## Keep copy of input
  
  isAmelogenin = function(x){
    ## Looks after amelogenin loci - and drop these
    amel.candidate <- unlist(lapply(head(x, n = 20), function(y){any(grepl("^(x|y|X|Y|xx|xy|XX|XY)$",paste(y)))}))
    
    return(amel.candidate)
  }
    
  x <- x[,!isAmelogenin(x)]
  
  ## 
  stripID = function(x){
    id = nid = NULL
    
    if(ncol(x) %% 2 != 0){ ## Unequal number of columns - assumes first col is identifier
      id <- x[ ,1, drop = TRUE]
      nid <- names(x)[1]
      x <- x[,-1]
    }else{
      id = 1:nrow(x)
      nid = NULL
    }
    
    return(list(x = x, id = id, nid = nid))
    
  }
  
  r = stripID(x)
  x = r$x
  id = r$id
  nid = r$nid
  
  ## number of loci
  numLoci <- ncol(x)/2
  
  toInteger = function(x){
    ## Converts all alleles to a*10, e.g. 9 -> 90 and 9.3 -> 93 (to deal with .1, .2 and .3 alleles)
    for(i in 1:ncol(x)){
      if(class(x[[i]])=="factor") x[[i]] <- paste(x[[i]])
      x[[i]] <- as.numeric(x[[i]])*10
      class(x[[i]]) <- "integer"
    }
    
    return(x)
  }
  
  x = toInteger(x)
  
  ## Specific profile(s) that should be compared to profiles in x
  if(!is.null(profiles)){
    nx <- names(x)
    if(is.null(dim(profiles))) 
      dim(profiles) <- c(1,length(profiles))
    
    if(is.matrix(profiles)) 
      profiles <- as.data.frame(profiles)
    
    #amel.candidate <- unlist(lapply(head(profiles,n=min(nrow(profiles),20)),function(y) any(grepl("^(x|y|X|Y|xx|xy|XX|XY)$",paste(y)))))
    profiles <- profiles[ ,!isAmelogenin(profiles)]
    
    if(ncol(profiles)==ncol(x)){ 
      names(profiles) <- nx
    }else{
      id <- c(profiles[,1,drop=TRUE],id)
      profiles <- structure(profiles[,-1],.Names=nx)
    }
    
#     for(i in 1:ncol(profiles)){
#       if(class(profiles[[i]])=="factor") profiles[[i]] <- paste(profiles[[i]])
#       profiles[[i]] <- as.numeric(profiles[[i]])*10
#       class(profiles[[i]]) <- "integer"
#     }
    profiles = toInteger(profiles)
    x <- rbind(profiles,x)
    single <- nrow(profiles) ## Number of profiles to compare with x // input for C++ code
  }else{ 
    single <- 0
  }
  
  orderAlleles = function(x){
    for(i in ((0:(numLoci-1))*2+1)){ ## Looks for situations where A[1]>A[2]: C++ code assumes A[1]<=A[2]
      if(any(swap.index <- x[[i]]>x[[i+1]])){
        if(Rallele) 
          swap.index[x[[i]]==990] <- FALSE ## If R-allele (coded as 99*10, cf above) don't swap (order matters)
        x[swap.index,c(i,i+1)] <- x[swap.index,c(i+1,i)]
      }
      ## Force homs to have wildcards
      if(wildcard.impose) 
         x[[i]][x[[i]]==x[[i+1]]] <- 0
    }
    return(x)
  }
  
  x = orderAlleles(x)
  
  if(!wildcard){
    ## If no wildcard is allowed, but the database contains "0" which is equivalent of "F" then these profiles are removed before comparison
    x0 <- x[(xnull <- apply(x[,-1],1,function(y) any(y==0))),]
    x <- x[!xnull,]
  }
  
  x.cp1 <- x ## Keep copy of x
  x <- do.call("paste",c(x=x,sep="\t")) ## Converts every line in the DB to a string separated by "\t"
  if(threads>1){
    res <- .Call("DNAtools_mcompare", x, as.integer(numLoci), as.integer(hit), as.integer(trace), as.integer(single), 
                                as.integer(threads), as.integer(wildcard), as.integer(wildcard.effect), as.integer(Rallele),
                                PACKAGE = "DNAtools")
  }else{
    res <- .Call("DNAtools_compare", x, as.integer(numLoci), as.integer(hit), as.integer(trace), as.integer(single), 
                               as.integer(wildcard), as.integer(wildcard.effect), as.integer(Rallele),
                               PACKAGE = "DNAtools")
  }
  
  if(wildcard.effect){
    result <- list(m = matrix(res$M, 2 * numLoci + 1, 2 * numLoci + 1, byrow = TRUE,
                              dimnames = list(Genuine = 0:(2 * numLoci), Wildcard=0:(2*numLoci))))
  }else{
    result <- list(m = matrix(res$M, numLoci+1, numLoci + 1,byrow = TRUE, 
                              dimnames = list(match = 0:numLoci, partial = 0:numLoci)))
  }
  
  call <- list(loci=numLoci,single=single,collapse=collapse,vector=vector,wildcard=c(wildcard,wildcard.effect))
  if(length(res$row1)>0){
    result <- c(result,list(hits=data.frame(id1=id[res$row1],id2=id[res$row2],match=res$matches,partial=res$partial)))
    if(wildcard) result$hits <- cbind(result$hits,Fmatch=res$fmatches,Fpartial=res$fpartial)
    names(result$hits)[1:2] <- paste(nid,1:2,sep="")
    call <- c(call,list(hit=hit))
  }
  else result$hits <- NULL
  if(collapse) result$m <- structure(dbCollapse(result$m),.Names=0:(2*numLoci))
  else if(vector) result$m <- structure(t(result$m)[up.tri(result$m)],
                                        .Names=t(outer(dimnames(result$m)[[1]],dimnames(result$m)[[2]],paste,sep="/"))[up.tri(result$m)])
  if(!wildcard){
    if(nrow(x0)>0) result$excludedProfiles <- x0
    else result$excludedProfiles <- "none"
  }
  attributes(result)$call <- call
  attributes(result)$class <- "dbcompare"
  result
}

print.dbcompare <- function(x,...){
  if(is.matrix(x$m)){
    if(attributes(x)$call$wildcard[2]) x$m[lower.tri(x$m)] <- NA
    else x$m[!up.tri(x$m)] <- NA
  }
  if(length(attributes(x)$names)==1){
    if(is.matrix(x$m)) print.table(x$m,...)
    else print.default(x$m,...)
  }
  else{
    if(is.matrix(x$m)){
      cat("Summary matrix\n")
      print.table(x$m,...)
    }
    else{
      cat("Summary vector\n")
      print.default(x$m,...)
    }
    if(!is.null(x$hits)){
      cat(paste("\nProfiles with at least",attributes(x)$call$hit,"matching loci\n"))
      x$hits <- x$hits[with(x$hits,order(match,partial,decreasing=TRUE)),]
      if("Fmatch" %in% names(x$hits)) x$hits <- x$hits[with(x$hits,order(match,partial,match,partial,decreasing=TRUE)),]
      rownames(x$hits) <- 1:nrow(x$hits)
      print.data.frame(x$hits,...)
    }
  }
}

