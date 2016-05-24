brt4df <- function(df,
                   varName,
                   max.order=1,
                   colNames,
                   auto=TRUE,
                   normalise=function(x) as.numeric(scale(log(x)))) {

  ## check that df is a data.frame
  if (!inherits(df,"data.frame"))
    stop("Argument df should be a data frame")

  ## check that varName is a df variable name
  if (!(varName %in% names(df)))
    stop(paste(varName,"should be a variable name of data frame df"))

  ## make sure that max.order makes sense
  max.order <- as.integer(max.order)
  if (max.order < 1) max.order <- as.integer(1)
  
  n <- dim(df)[1]
  if (auto) {
    eventIdx <- (1:n)[df$event==1]
  } else {
    eventIdx <- (1:n)[df[[varName]]==0]
    eventIdx <- eventIdx[!is.na(eventIdx)][-1]
  }
  
  lastOne <- eventIdx[length(eventIdx)]
  st <- df[,varName]
  rm(df)
  result <- matrix(integer(n*max.order),
                   nrow=n,
                   ncol=max.order)
  
  for (ordIdx in 1:max.order) {

    result[1:eventIdx[ordIdx],ordIdx] <- NA
    
    if (lastOne != n) {
      if (auto) offset <- st[lastOne]
      else offset <- st[lastOne-1]+1
      result[(lastOne+1):n,ordIdx] <- offset+seq(n-lastOne)
    }
    
    for (i in seq(along=eventIdx)[-(1:ordIdx)]) {
      mySeq <- (eventIdx[i-1]+1):eventIdx[i]
      if (auto) {
        if (ordIdx == 1) offset <- st[eventIdx[i-1]]
        else offset <- result[eventIdx[i-1],ordIdx-1]
      } else {
        if (ordIdx == 1) offset <- st[eventIdx[i-1]-1] + 1
        else offset <- result[eventIdx[i-1]-1,ordIdx-1] + 1
      }
      result[mySeq,ordIdx] <-  offset + seq(along=mySeq)
    } ## End of the loop on i

  } ## End of for loop on ordIdx
  
  if (inherits(normalise,"function")) {
    result <- apply(result,2,normalise)
    result <- cbind(normalise(st),result)
    if (missing(colNames)) colNames <- c(paste("r",varName,sep=""),
                                         paste("r",varName,".",(1:max.order)+1,sep=""))
    colnames(result) <- colNames
  } else {
    if (missing(colNames)) colNames <- paste(varName,".",(1:max.order)+1,sep="")
    colnames(result) <- colNames
  }
  
  as.data.frame(result)

}
