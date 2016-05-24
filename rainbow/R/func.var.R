func.var = function (functions) 
{
    nrow = dim(functions)[1]
    (nrow - 1) * apply(functions, 2, var) / nrow
}
