"optim.fun4" <-
function(x, data1, data2, param1, param2)
{
# Numerical maximum likelihood.
L1 <- x[1]
L2 <- x[2]
L3 <- x[3]
L4 <- x[4]
M1 <- x[5]
M2 <- x[6]
M3 <- x[7]
M4 <- x[8]
p1 <- x[9]
data1 <- sort(data1)
check1 <- gl.check.lambda.alt(L1, L2, L3, L4, param1)
if(check1==TRUE){
check2 <- qgl(0, L1, L2, L3, L4, param1) <= min(data1)
check3 <- qgl(1, L1, L2, L3, L4, param1) >= max(data1)
    }
    else {
check2<-FALSE
check3<-FALSE}
data2 <- sort(data2)
check4 <- gl.check.lambda.alt(M1, M2, M3, M4, param2)

if(check4==TRUE)
{
check5 <- qgl(0, M1, M2, M3, M4, param2) <= min(data2)
check6 <- qgl(1, M1, M2, M3, M4, param2) >= max(data2)
}
else {
check5<-FALSE
check6<-FALSE
}


check7 <- p1 > 0
check8 <- p1 < 1

if(as.logical(check1 * check2 * check3 * check4 * check5 * check6 * check7 * check8)){
response<- -sum(log(p1)+log(dgl(data1,L1,L2,L3,L4,param1)))-sum(log(1-p1)+log(dgl(data2,M1,M2,M3,M4,param2)))
}


else if(as.logical(check1 * check2 * check3 * check4 * check5 * check6 * check7 * check8) == 0) {
response <- NA
}
return(response)
}

