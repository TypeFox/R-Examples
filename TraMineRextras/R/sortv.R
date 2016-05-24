## computes a sorting key from successive states
## author: Gilbert Ritschard

sorti <- function(seqdata, start = "end", sort.index = "TRUE"){
    seqdata[seqdata==attr(seqdata,"void")] <- attr(seqdata,"nr")

	if (start %in% c("beg", "end")) {
    	end <- if (start=="end") { max(seqlength(seqdata)) } else { 1 }
    	beg <- if (start=="end") { 1 } else { max(seqlength(seqdata)) }
    }
    else{
        warning('start should be one of "beg" or "end" ')
        return(k <- 1:max(seqlength(seqdata)))
    }

    k <- do.call(order, as.data.frame(seqdata)[,end:beg])
    if (sort.index){
        return(k)
    }
    else{
        return(order(k))
    }
}


sortv <- function(seqdata, start = "end"){
    return(sorti(seqdata, start=start, sort.index = "FALSE"))
    }
##    asize <- length(alphabet(seqdata))
##    dframe <- as.data.frame(seqdata)
##    if ("%" %in% levels(dframe[,max(beg,end)]))
##    k <- as.numeric(dframe[,beg])
##    for (i in (beg:end)){
##        k <- k + asize^i*as.numeric(dframe[,i])
##    }
##    return(k)

