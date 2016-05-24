optim.fun3.C.m <-
function (lambdas, data, param,inverse.eps = 1e-08, 
max.iterations = 500)
{
    
    u<-rep(0,length=nrow(lambdas))
    npar<-nrow(lambdas)

    result <- .C("optim_fun3_v", param, as.double(lambdas[,1]), as.double(lambdas[,2]),
as.double(lambdas[,3]), as.double(lambdas[,4]),
    as.double(inverse.eps), as.integer(max.iterations), as.double(data),
    as.double(u), as.integer(length(data)),
    as.double(.Machine$double.eps),as.integer(npar))

    return(result[[9]])

}
