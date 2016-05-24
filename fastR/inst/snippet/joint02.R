require(cubature)
f <- function(x) { 6 * x[1] * x[2]^2 }
adaptIntegrate(f,c(0,0),c(1,1))
g <- function(x) { 
    if (x[1] > x[2]) {return(0)}   # set value to 0 if X > Y
    return(f(x))                   # else return joint pdf
    }
adaptIntegrate(g,c(0,0),c(1,1),tol=0.01)  # get less accuracy
