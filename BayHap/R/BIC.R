`BIC` <-
function(res){

     bic<--2*res[[1]]$versem_max+(res[[1]]$n_coeff*log(res[[1]]$numindiv))
     bic

}

