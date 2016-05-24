nk <- function(N, m)
{
Q <- SupportWR(N, m, ID = FALSE)
I <- matrix(0, choose(N+m-1, m), N)
for (i in 1:m) {
for (j in 1:choose(N+m-1, m)) {
for (k in 1:N) {
if (Q[j, i] == k)
I[j, k] <- sum(as.double(Q[j,]==k))
}
}
}
I
}