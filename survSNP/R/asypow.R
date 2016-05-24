asypow <-
function(n,theta,a,b,lambda0,q,p,alpha,z,exactvar)
     .Call("asypowRcpp", n, theta, a, b, lambda0, q, p, alpha, z, exactvar)

