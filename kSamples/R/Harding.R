Harding <-
function (nn) 
{
nn <- sort(as.integer(nn))
nvec <- rev(cumsum(rev(nn)))
k <- length(nn)
L1 <- sum(nn[1:(k-1)]*nvec[2:k])+1
freq <- double(L1)

out <- .C("Harding0", k=as.integer(k), L1=as.integer(L1),
		nn=as.integer(nn), nvec=as.integer(nvec),
		freq=as.double(freq) )

freq <- out$freq

if(is.nan(sum(freq)) ||
        abs(sum(freq)-1) > .0000001) freq <- NaN

freq
}

