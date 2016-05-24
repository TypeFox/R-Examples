head2<-function(x, head = 3, tail=1, dotrows=1)
{
  if(!(is.matrix(x) || is.data.frame(x))){stop("A data frame or matrix is required")}
  # check number of rows in matrix or dataframe
  if(nrow(x)< head+dotrows+tail){x}
  else
  {
    ## fix -- tail adds rownames to matrix using addrownums, but not head 
    h1<-head(x,head + dotrows)
    if(is.null(rownames(h1))) { rownames(h1) <- paste("[", 1:(head+dotrows), ",]", sep = "") }
    x <- format(rbind(h1, tail(x,tail)))
    if(dotrows>0)
    {
      x[(head + 1):(head + dotrows),] <- "."
      for(i in 1:dotrows){rownames(x)[head+i]<-paste(".", substring("         ", 1, i-1))}
    }
    ## remove quotes from matrix elements?
    noquote(x)
  }
}

