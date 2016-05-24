ChiSN <-
function(x){
k <- 1.5477; z <- x/k
(3*z^2 - 3*z^4 +z^6)*(abs(x) <= k) + 1*(abs(x) > k)}

