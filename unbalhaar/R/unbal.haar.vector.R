unbal.haar.vector <-
function(a) {
n <- a[3] - a[1] + 1
m <- a[2] - a[1] + 1

return(c(rep(sqrt(1/m - 1/n), m), rep(-1/sqrt(n^2/m - n), n-m)))

}

