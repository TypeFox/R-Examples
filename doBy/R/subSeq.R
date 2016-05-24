subSeq <- function (x, item = NULL) {
    rrr <- rle(x)
    len <- rrr$lengths
    val <- rrr$values

    first <- last <- rep.int(NA, length(val))
    first[1] <- 1
    last [1] <- len[1]
    if (length(val)>1){
	for (kk in 2:length(val)){
            first[kk] <- last[kk-1]+1
            last [kk] <- last[kk-1]+len[kk]
	}
    }
    midp <- floor(first + len/2)

    ans <- cbind(first=first, last=last, slength=len, midpoint=midp)

    if (!is.null(item)) {
        iii <- val==item
        ans <- as.data.frame.matrix(ans[iii,,drop=FALSE], stringsAsFactors=FALSE)
        ans$value <- val[iii]
    } else {
        ans <- as.data.frame.matrix(ans, stringsAsFactors=FALSE)
        ans$value <- val
    }
    ans
}



## Obsolete version; way too slow...

## subSeq <- function(x, item=NULL){
##   first  <-  c(1,which(c(FALSE,
##                          mapply(function(a,b){!identical(a,b)},
##                                 x[-1],x[-length(x)])
##                          )
##                        )
##                )
##   last   <-  c( first[-1] - 1, length(x))
##   rlen   <-  (last-first)+1
##   value  <-  x[first]
##   midp   <-  floor(first + rlen/2)
##   ans    <-  data.frame(value=value, first=first, last=last, slength=rlen,
##                         midpoint=midp, stringsAsFactors=FALSE)

##   if (!is.null(item)){
##     if (item %in% value){
##       iii <- unlist(lapply(as.character(ans$value), identical, as.character(item)))
##       ans <- ans[iii,]
##     } else {
##       ans <- NULL
##     }
##   }
##   return(ans)
## }


# Example:
# x <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 3)
# subseq(x)
# subseq(x,item=1)

# x <- c(T,T,F,F,F,T,F,T,T,T,T)
# subseq(x)
# subseq(x,item=TRUE)





