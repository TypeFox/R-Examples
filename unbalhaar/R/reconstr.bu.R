reconstr.bu <-
function(buh.bu) {
y <- buh.bu
n <- dim(y$detail)[2] + 1

x <- y$smooth
detail <- y$detail

for (i in 1:(n-1)) {

ind <- round(detail[1,n-i])
x[(ind+2):n] <- x[(ind+1):(n-1)]
new.left <- detail[2,n-i] * x[ind] + sqrt(1 - detail[2,n-i]^2) * detail[3,n-i]
new.right <- sqrt(1 - detail[2,n-i]^2) * x[ind] - detail[2,n-i] * detail[3,n-i]
x[ind] <- new.left
x[ind+1] <- new.right

}
return(x)

}
