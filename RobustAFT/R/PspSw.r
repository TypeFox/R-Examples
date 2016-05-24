PspSw <-
function(x){
k <- 1.718; k2  <- k*k; k4 <- k2*k2; k6 <- k2*k4
x2  <- x*x; x4 <- x2*x2
(6/k2-36*x2/k4+30*x4/k6)*(abs(x)<k)}

