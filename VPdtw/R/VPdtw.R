VPdtw <- function(reference,query,penalty=0,maxshift=50,Reference.type=c("random","median","mean","trimmed")) {

  ## We assume Sakoe Chiba DTW to allow for faster computation times
  if(!is.numeric(maxshift)) {
    return("Please specify maxshift as an integer value\n")
  }

  ## Figure out what kind of alignment we are doing -

  ## a reference vector to a query vector?
  ## no reference and a matrix of query vectors?
  ## something else = no implementation here
  
  if(is.null(reference) & !is.matrix(query)) {
    cat("Please specify a reference when passing a non-matrix query\n")
    return("Exiting....\n")
  }
  
  if(is.null(reference) & is.matrix(query)) {

    ## If no reference is specified, we choose one randomly, using
    ## median, mean or trimmed mean of the query matrix depending on
    ## the value of type as specified by the user.
    
    type <- match.arg(Reference.type,c("random","median","mean","trimmed"))
    ss <- sample(1:ncol(query),1)
    reference <- switch(type,
                        random = query[,ss],
                        median = apply(query,1,median,na.rm=TRUE),
                        mean = apply(query,1,mean,na.rm=TRUE),
                        trimmed = apply(query,1,mean,na.rm=TRUE,trim=0.1))
    reference <- na.omit(reference)

    ## if penalty is a number, then this means a constant penalty vector, create that vector
    if(length(penalty)==1) penalty <- rep(penalty,length(reference))

    ## Check penalty vector length
    if(length(penalty)<length(reference)) {
      warning("Penalty vector should be at least of length ",length(reference)," but it has length ",length(penalty),"\n")
      return("Exiting...\n")
    }

    ## information used in summary at end - what kind of reference do
    ## we have - in this block of code it IS null, include info on
    ## what kind was generated.
    
    information <- paste("Reference is NULL.")
    information <- paste(information,switch(type,
                                            random = paste("Query column #",ss,"is chosen at random."),
                                            median = "Query median used.",
                                            mean = "Query mean used.",
                                            trimmed = "Query trimmed mean used."))
    information <- paste(information,"\n")
  } else {
    information <- paste("Reference supplied by user.\n")
  }

  ## This function DoAlignment is a wrapped that takes a vector query,
  ## vector reference and vector penalty and (eventually) does a Sakoe
  ## Chiba DTW with width maxshift. This function calls
  ## signalMatchABand which calls the C++ code itself. This can be
  ## "applied" to matrices of queries. Another version DoAlignmentP is
  ## for apply to matrix of penalties.
  
  DoAlignment <- function(query,reference,penalty,maxshift) {
    ## Drop NAs
    reference <- na.omit(reference)
    query <- na.omit(query)
    result <- signalMatchABand(reference=reference,
                               query=query,
                               lambda=penalty,
                               maxshift=maxshift)
    result ## matrix four columns 
  }

  ## Ok we are ready to do our alignment.

  ## Scenario 1: all vectors:
  if(is.vector(reference) & is.vector(query) & is.vector(penalty)) {
    information <- "Reference is supplied by the user.\n"
    information <- paste(information,"Query vector is of length ",length(query),".\n",sep="")

    information <- paste(information,"Single Penalty vector supplied by user.\n",sep="")
    information <- paste(information,"Max allowed shift is ",maxshift,".\n",sep="")
    reference <- na.omit(reference)

    if(length(penalty)==1) penalty <- rep(penalty,length(reference))
    
    result <- DoAlignment(query,reference,penalty,maxshift)

    ## format output a little better:
    
    output <- vector("list",6)
    names(output) <- c("xVals","reference","query","penalty","warpedQuery","shift")
    
    output$xVals <- result[,"xVals"]
    
    output$reference <- result[,"reference"]

    output$query <- query
    output$penalty <- penalty

    output$warpedQuery <- result[,"warped query"]

    output$shift <- result[,"shift"]
    
    class(output) <- "VPdtw"

    ## Summary Statistics: 

    cost1 <- function(x) {
      ret <- c(sum(abs(x$warpedQuery - x$reference),na.rm=TRUE) +
               sum(x$penalty[x$xVals[which(diff(x$shift)==1)+1]],na.rm=TRUE) +
               2*sum(x$penalty[x$xVals[which(diff(x$shift)==-1)+1]],na.rm=TRUE),
               sum(!is.na(x$warpedQuery * x$reference),na.rm=TRUE),
               max(abs(x$shift),na.rm=TRUE),
               sum(diff(x$shift)==0,na.rm=TRUE)+1,
               sum(diff(x$shift)==1,na.rm=TRUE),
               sum(diff(x$shift)==-1,na.rm=TRUE))
      names(ret) <- c("Cost", "Overlap", "Max Obs. Shift","# Diag Moves", "# Expanded","# Dropped")
      ret
    }
    output$summary <- cost1(output)
              
    output$information <- information
    
    return(invisible(output))
  }

  ## Scenario 2: Matrix of queries, penalty vector
  
  if(is.matrix(query) & is.vector(penalty)) {
    information <- paste(information,"Query matrix is made up of ",ncol(query)," samples of length ",nrow(query),".\n",sep="")
    information <- paste(information,"Single Penalty vector supplied by user.\n",sep="")
    information <- paste(information,"Max allowed shift is ",maxshift,".\n",sep="")
    reference <- na.omit(reference)
    ## Align all of query vectors to reference using the same penalty
    queryL <- as.list(as.data.frame(query)) ## to ensure I get a list in the next round, can't apply to matrix
    result <- lapply(queryL,DoAlignment,reference=reference,penalty=penalty,maxshift=maxshift) ## produces a list

    ## result is now a list of length ncol(penalty) each part is a matrix

    xlim <- NULL
    for(ii in 1:length(result)) xlim <- c(xlim,range(result[[ii]][,1]))
    xlim <- range(xlim)

    xVals <- seq(xlim[1],xlim[2],by=1)

    output <- vector("list",6)
    names(output) <- c("xVals","reference","query","penalty","warpedQuery","shift")
    
    output$xVals <- xVals
    
    str <- which(xVals==1)
    end <- which(xVals==length(reference))
    output$reference <- rep(NA,length(xVals))
    output$reference[seq(str,end,by=1)] <- reference

    output$query <- query
    output$penalty <- penalty
    
    output$warpedQuery <- matrix(NA,length(xVals),ncol(query))
    colnames(output$warpedQuery) <- paste("warped query",1:ncol(query))

    output$shift <- matrix(NA,length(xVals),ncol(query))
    colnames(output$shift) <- paste("shift",1:ncol(query))
    
    for(ii in 1:ncol(query)) {
      colName <- paste("warped query",ii)
      str <- which(xVals==result[[ii]][1,1])
      end <- which(xVals==result[[ii]][nrow(result[[ii]]),1])
      output$warpedQuery[seq(str,end,by=1),colName] <- result[[ii]][,3]
      colName <- paste("shift",ii)
      output$shift[seq(str,end,by=1),colName] <- result[[ii]][,4]
    }

    ## matplot(output$xVals,output$warpedQuery,type="n",lty=1,col=c(2,3))
    ## lines(output$xVals,output$reference,lwd=2,col=1)
    ## matplot(output$xVals,output$warpedQuery,type="l",lty=1,col=c(2,3),add=TRUE)
    
    class(output) <- "VPdtw"

    ## Summary Statistics for each query separately
    
    cost2 <- function(x,ii) {
      ret <- c(sum(abs(x$warpedQuery[,ii] - x$reference),na.rm=TRUE) +
               sum(x$penalty[x$xVals[which(diff(x$shift[,ii])==1)+1]],na.rm=TRUE) +
               2*sum(x$penalty[x$xVals[which(diff(x$shift[,ii])==-1)+1]],na.rm=TRUE),
               sum(!is.na(x$warpedQuery[,ii] * x$reference),na.rm=TRUE),
               max(abs(x$shift[,ii]),na.rm=TRUE),
               sum(diff(x$shift[,ii])==0,na.rm=TRUE)+1,
               sum(diff(x$shift[,ii])==1,na.rm=TRUE),
                 sum(diff(x$shift[,ii])==-1,na.rm=TRUE))
      names(ret) <- c("Cost", "Overlap", "Max Obs Shift","# Diag Moves", "# Expanded","# Dropped")
      ret
    }
    
    output$summary <- NULL
    for(ii in 1:ncol(output$warpedQuery)) output$summary <- rbind(output$summary,cost2(output,ii))
    rownames(output$summary) <- paste("Query #",1:ncol(output$warpedQuery),":",sep="")
    output$information <- information
        
    return(invisible(output))
  }

  ## For doing alignment for many different penalties
  DoAlignmentP <- function(penalty,query,reference,maxshift) {
    ## Drop NAs
    reference <- na.omit(reference)
    query <- na.omit(query)
    result <- signalMatchABand(reference=reference,
                               query=query,
                               lambda=penalty,
                               maxshift=maxshift)
    result
  }

  ## Scenario 3: vector query and penalty matrix
  if(is.vector(query) & is.matrix(penalty)) {

    information <- paste(information,"Query vector of length ",length(query),".\n",sep="")
    
    information <- paste(information,"Penalty matrix made up of ",ncol(penalty)," penalties supplied by user.\n",sep="")

    information <- paste(information,"Max allowed shift is ",maxshift,".\n",sep="")
    reference <- na.omit(reference)
    ## Align query to reference using each of the penalties separately
    penaltyL <- as.list(as.data.frame(penalty)) ## to ensure I get a list in the next round, can't apply to matrix
    result <- lapply(penaltyL,DoAlignmentP,reference=reference,query=query,maxshift=maxshift) ## produces a list
    ## result is now a list of length ncol(penalty) each part is a matrix

    xlim <- NULL
    for(ii in 1:length(result)) xlim <- c(xlim,range(result[[ii]][,1]))
    xlim <- range(xlim)

    xVals <- seq(xlim[1],xlim[2],by=1)

    output <- vector("list",6)
    names(output) <- c("xVals","reference","query","penalty","warpedQuery","shift")
    
    output$xVals <- xVals
    
    str <- which(xVals==1)
    end <- which(xVals==length(reference))
    output$reference <- rep(NA,length(xVals))
    output$reference[seq(str,end,by=1)] <- reference

    output$query <- query
    output$penalty <- penalty
    
    output$warpedQuery <- matrix(NA,length(xVals),ncol(penalty))
    colnames(output$warpedQuery) <- paste("warped query penalty",1:ncol(penalty))

    output$shift <- matrix(NA,length(xVals),ncol(penalty))
    colnames(output$shift) <- paste("shift penalty",1:ncol(penalty))
    
    for(ii in 1:ncol(penalty)) {
      colName <- paste("warped query penalty",ii)
      str <- which(xVals==result[[ii]][1,1])
      end <- which(xVals==result[[ii]][nrow(result[[ii]]),1])
      output$warpedQuery[seq(str,end,by=1),colName] <- result[[ii]][,3]
      colName <- paste("shift penalty",ii)
      output$shift[seq(str,end,by=1),colName] <- result[[ii]][,4]
    }

    ## matplot(output$xVals,output$warpedQuery,type="n",lty=1,col=c(2,3))
    ## lines(output$xVals,output$reference,lwd=2,col=1)
    ## matplot(output$xVals,output$warpedQuery,type="l",lty=1,col=c(2,3),add=TRUE)
    class(output) <- "VPdtw"

    ## Summary Statistics for each query separately
    
    cost3 <- function(x,ii) {
      ret <- c(sum(abs(x$warpedQuery[,ii] - x$reference),na.rm=TRUE) +
               sum(x$penalty[,ii][x$xVals[which(diff(x$shift[,ii])==1)+1]],na.rm=TRUE) +
               2*sum(x$penalty[,ii][x$xVals[which(diff(x$shift[,ii])==-1)+1]],na.rm=TRUE),
               sum(!is.na(x$warpedQuery[,ii] * x$reference),na.rm=TRUE),
               max(abs(x$shift[,ii]),na.rm=TRUE),
               sum(diff(x$shift[,ii])==0,na.rm=TRUE)+1,
               sum(diff(x$shift[,ii])==1,na.rm=TRUE),
               sum(diff(x$shift[,ii])==-1,na.rm=TRUE))
      names(ret) <- c("Cost", "Overlap", "Max Obs Shift", "# Diag Moves", "# Expanded","# Dropped")
      ret
    }
    
    output$summary <- NULL
    for(ii in 1:ncol(output$warpedQuery)) output$summary <- rbind(output$summary,cost3(output,ii))
    rownames(output$summary) <- paste("Penalty #",1:ncol(output$warpedQuery),":",sep="")
    output$information <- information
    
    return(invisible(output))
  }

  if(is.matrix(query) & is.matrix(penalty)) {
    cat("Multiple queries and multiple penalties not yet implemented\n")
    cat("Please create loops and call VPdtw as required\n")
    return("Exiting....\n")
  }

  ## finished
}

print.VPdtw <- function(x,...) {
  cat(x$information)
  cat("\n")
  print(signif(x$summary,5))
}
