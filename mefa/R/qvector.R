`qvector` <-
function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, ...)
{
nrun <- length(probs)-1
qa <- quantile(x, probs, na.rm = na.rm, ...)
out <- rep(probs[2], length(x))
for (i in 2:nrun)
    out[which(x > qa[i] & x <= qa[(i+1)])] <- probs[(i+1)]
return(out)
}

