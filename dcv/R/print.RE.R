print.RE <-
function(x, ...){
cat('    Reduction of Error(RE)=', round(x[[1]],3),' \n 
    Mean squared erro(MSE)=', round(x[[2]],3),' \n
    Root mean squared error(RMSE)=', round(x[[3]],3),'\n')
}

