hard.thresh.bu <-
function(buh.bu, sigma = 1) {
n <- dim(buh.bu$detail)[2]
th <- sigma * sqrt(2 * log(n))
buh.bu$detail[3,] <- buh.bu$detail[3,] * (abs(buh.bu$detail[3,]) > th)
return(buh.bu)
}

