"bgrPoint" <-
function(sample)
# Calculate the bgr statistic given a sample concatenated over chains 
{
    numChains <- getNumChains()
    sampleSize <- length(sample)
    lenChain <- sampleSize %/% numChains
    if (is.R())
      dq <- quantile(sample, c(0.1, 0.9), names = FALSE)
    else
      dq <- quantile(sample, c(0.1, 0.9))
    d.delta <- dq[2] - dq[1]
    n.delta <- 0
    for (i in 1:numChains) {
        if (is.R())
          nq <- quantile(sample[((i - 1) * lenChain + 1) : (i * lenChain)], c(0.1, 0.9), names = FALSE)
        else
          nq <- quantile(sample[((i - 1) * lenChain + 1) : (i * lenChain)], c(0.1, 0.9))
        n.delta <- n.delta + nq[2] - nq[1]
    }
    n.delta <- n.delta / numChains
    bgr.stat <- d.delta / n.delta
    return(c(n.delta, d.delta, bgr.stat))
}
