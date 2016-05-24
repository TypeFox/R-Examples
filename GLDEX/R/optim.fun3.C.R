optim.fun3.C <-
function (lambdas, data, param,inverse.eps = 1e-08, 
max.iterations = 500)
{
    
    u<-(-1)

    result <- .C("optim_fun3", param, as.double(lambdas[1]), as.double(lambdas[2]),
as.double(lambdas[3]), as.double(lambdas[4]),
    as.double(inverse.eps), as.integer(max.iterations), as.double(data),
    as.double(u), as.integer(length(data)),
    as.double(.Machine$double.eps))

    return(result[[9]])

}
