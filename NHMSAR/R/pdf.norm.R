pdf.norm <-
function (x, mean, sigma) 
{   
    d = dim(x)[1] 
    #d = dim(sigma)[1]
    T = dim(x)[2]
    if (is.null(d) | is.na(d)) {d = 1 ; T = length(x)}
    tmp = matrix(t(x-mean),T,d)
    mahal <- apply((tmp%*%solve(sigma))*tmp,1,sum)
    numer <- exp(-0.5*mahal)
    denom = (2*pi)^(d/2)*sqrt(abs(det(sigma)))

    numer/denom
}
