jGradKernelFitPdf2 <-
function(s,h,x){
f <- 1/h * (-1/h+((x-s)/h)^2)*jNormPdf((x-s)/h)
return(mean(f)^2)
}
