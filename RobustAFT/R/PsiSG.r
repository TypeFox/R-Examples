PsiSG <-
function(x){
k <- 1.5477; z <- x/k
(6*z-12*z^3+6*z^5)/k*(abs(x) <= k) + 0*(abs(x) > k)}

