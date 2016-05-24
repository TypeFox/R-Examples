tscale <-
function(p)
{
if (length(p) != 3) {print("tscale must have input vector of length 3")}
p <- p/sum(p)
for (i in 1:3)
{
p[i] <- max(0,p[i])
p[i] <- min(1,p[i])
}
p <- p/sum(p)
p
}
