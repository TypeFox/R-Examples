##' @export
jumptimes <- function(time, status=TRUE, 
                      id,cause,
                      sample,sample.all=TRUE,
                      strata=NULL,num=NULL, ...) {
    if (missing(id)) {
        time <- if (missing(cause)) time[status>0]
                else time[status==cause]
    } else {
        ww <- na.omit(fast.reshape(cbind(time=time,status=status),id=id,num=num))
        statusvar <- grep("status",colnames(ww))
        timevar <- grep("time",colnames(ww))
        if (missing(cause)) {
            idx <- which(rowSums(ww[,statusvar]>0)==length(statusvar))
        } else {
            idx <- which(apply(as.matrix(ww[,statusvar]),1,function(x) all(x==cause)))
        }
            time <- na.omit(do.call(pmax,as.list(ww[idx,timevar])))
    }

    if (!missing(sample) && sample<length(time)) {
        time <- sort(unique(time))
        t0 <- seq(min(time),max(time),length.out=min(sample,length(time)))
        ii <- fast.approx(time,t0)
        ii <- unique(ii)
        time0 <- time[ii]
        if (length(ii)<sample && sample.all) {
            remain <- setdiff(seq(length(time)),ii)
            time0 <- c(time0,time[base::sample(remain,sample-length(ii))])
        }
        time <- time0
    } 
    return(sort(time))
}

## with(prt, jumptimes(time,status==2,id=id,sample=10))
