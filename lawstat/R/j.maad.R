`j.maad` <-
function (x) 
{

DNAME=deparse(substitute(x))

## Strip NAs
x<-na.omit(x)


### Robust Standard Deviation J
J<-sqrt(pi/2)*mean(abs(x-median(x))) 

# return(J)

paste("MAAD estimated J =", J, "for data", DNAME)

}

