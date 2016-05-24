"optim.fun3" <-
function(x, data, param)
{
# Numerical maximum likelihood.
L1 <- x[1]
L2 <- x[2]
L3 <- x[3]
L4 <- x[4]
data <- sort(data)
check1 <- gl.check.lambda.alt(L1, L2, L3, L4, param)

if(check1==TRUE){
check2 <- qgl(0, L1, L2, L3, L4, param) <= min(data)
check3 <- qgl(1, L1, L2, L3, L4, param) >= max(data)
    }

    else {check2<-FALSE
       check3<-FALSE}

if(as.logical(check1 * check2 * check3)) {
response <-  - sum(log(dgl(data, L1, L2, L3, L4, param)
))
}
else if(as.logical(check1 * check2 * check3) == 0) {
response <- NA
}
return(response)
}

