w.hist <-
function(x, mu, breaks, cuts, ylim = NULL, freq = NULL, lwd = NULL) 
{
Gray = paste("gray", round(seq(from=10, to=100, len=length(cuts)-1)), sep="")
hist(x, col=Gray[1], xlim=range(breaks), ylim=ylim, breaks = breaks, freq = freq, lwd=lwd )
i=2
while(i<=length(cuts))
{
X=x
X[(X*(mu>=cuts[i]))==0]=NA
hist(X, col=Gray[i], xlim=range(breaks), ylim=ylim, breaks=breaks, freq=freq, lwd=lwd, add=TRUE)
i=i+1
}
}
