cOpt_gsym2_kernel <-
function(t, parameters)
{
   cte = parameters[1]
   w0b = parameters[2]
   x0b = parameters[3:length(parameters)]
   ROCt = 1-cdf_gaussianKernel(x0b, 1-t, w0b)
   res <- (1/cte)*(1-t)-ROCt
   return(res)
}
