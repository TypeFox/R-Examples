dist2full <-
function(dis)
{
#### dis is a vector of distances between n points
### returns a symmetric matrix of size n x n
     n <- attr(dis, "Size")
     full <- matrix(0, n, n)
     full[lower.tri(full)] <- dis
     full + t(full)
}
