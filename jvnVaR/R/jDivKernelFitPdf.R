jDivKernelFitPdf <-
function(s,h,x){
f <- 1/h * jNormPdf((x-s)/h) * ((s-x)/h)
return(mean(f))
}
