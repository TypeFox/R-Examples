lower.record.values <-
function(sqnc,k) {
    if (sum(is.na(sqnc))>0) {
        stop("There are missing values in data. The function is not designated for such a data!")
    } else {
        if (length(k)==1 && k>=1 && k<=length(sqnc) && is.numeric(k) && floor(k)==k) {
            sqnc<--unique(sqnc)
            dl<-length(sqnc)
            wndw<-sort(sqnc[1:k])
            records<-wndw[1]
            if (k<dl) {
                for (i in (k+1):dl) {
                    if (sqnc[i]>wndw[1]) {
                        if (i<=dl-k+2) {
                            w<-which(sqnc[i]>wndw)
                        } else {
                            w<-which(sqnc[i]>wndw[1:(dl-i+2)])
                        }
                        where<-w[length(w)]
                        if (where>1) {
                            wndw[1:(where-1)]<-wndw[2:where]
                        }
                        wndw[where]<-sqnc[i]
                        records<-c(records,wndw[1])
                    }
                }
            }
            return(-records)
        } else {
            stop("k must be an integer between 1 and length(sqnc)")
        }
    }
}

