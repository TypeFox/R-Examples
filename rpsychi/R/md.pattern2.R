md.pattern2 <- function (x)
{
    check1 <- sum(is.na(x))==0
    if(check1){
      x <- rbind(x, NA)
    }
    
    
    n <- nrow(x)
    p <- ncol(x)
    r <- 1 * is.na(x)
    rev.r <- 1 * !is.na(x)
    nmis <- colSums(r)

    pat.r <- NA
    for(i in 1:ncol(r)){
          pat.r <- paste(pat.r, rev.r[,i])
          }

    pat.r <- factor(pat.r)
    n.pat <- nlevels(pat.r)
    pat.table <- table(pat.r)

    out.pat <-  unlist(strsplit(substring(names(pat.table), 4), " "))
    out.pat <-  matrix(as.numeric(out.pat), ncol=ncol(x), nrow=n.pat, byrow=TRUE)
    row.pat <-  rowSums(out.pat == 0)
    temp <- cbind(out.pat, as.integer(pat.table), row.pat )
    temp <- temp[order(temp[,ncol(temp)-1],decreasing = TRUE),]
    
    ## output
    output <- matrix(ncol=ncol(x)+1, nrow=n.pat+1)
    colnames(output) <- c(names(sort(nmis, decreasing=TRUE)), "NA")
    rownames(output) <- c(temp[,ncol(temp)-1], "Sum")
    output["Sum",] <- c(sort(nmis, decreasing=TRUE), sum(nmis))
    output[,"NA"]  <- c(temp[,ncol(temp)], sum(nmis))
    temp <- temp[,-c(ncol(temp)-1,ncol(temp))]
    temp <- temp[,order(nmis, decreasing=TRUE)]
    output[1:n.pat,1:ncol(x)] <- temp
    
    if(check1){
        output <- output[-2, ]
        output[2,] <- 0
    }    
    
    return(output)
}
