"optim.fun2" <-
function(x, data, param, breaks, yy)
{
l1 <- x[1]
l2 <- x[2]
l3 <- x[3]
l4 <- x[4]
data <- sort(data)
check1 <- gl.check.lambda.alt(l1, l2, l3, l4, param)
if(check1==TRUE){
check2 <- qgl(0, l1, l2, l3, l4, param) <= min(data)
check3 <- qgl(1, l1, l2, l3, l4, param) >= max(data)
    }

    else {check2<-FALSE
       check3<-FALSE}
if(as.logical(check1 * check2 * check3)) {
# More efficient to put these outside the function so they are not evaluated everytime optim is invoked!
xx <- diff(pgl(breaks, l1, l2, l3, l4, param = param))
# The response is deliberately weighted so that the ones with highest frequencies need to have higher priority to be minimized. 
response <- sum(yy * (xx - yy)^2)
}
else response <- NA
return(response)
}

