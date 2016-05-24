"qcochran" <-
function(p,n,k) {

f <- qf((1-p)/k,(n-1)*(k-1),n-1);
c <- 1/(1+(k-1)*f)

return(c)

}

