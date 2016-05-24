lower.record.times <-
function(sqnc,k) {
    if (sum(is.na(sqnc))>0) {
        stop("There are missing values in data. The function is not designated for such a data!")
    } else {
        if (length(k)==1 && k>=1 && k<=length(sqnc) && is.numeric(k) && floor(k)==k) {
            records<-lower.record.values(sqnc,k)
            dl<-length(records)
            L<-numeric(dl)
            L[1]<-k
            if (dl>1) {
                for (i in 2:(dl)) {
                    L[i]<-L[i-1]+which(sqnc[(L[i-1]+1):length(sqnc)]<records[i-1])[1]
                }
            }
            return(L)
        } else {
            stop("k must be an integer between 1 and length(sqnc)")
        }
    }
}

