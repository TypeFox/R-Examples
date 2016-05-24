P_mass <-
function(N,m,n,x) exp(lchoose(N,x)+lchoose(N-x,m-x)+lchoose(N-m,n-x)
                                -lchoose(N,m)-lchoose(N,n))
