ptgpd <-
function(x,y,z, model="log", mar1=c(0,1,0.1), mar2=c(0,1,0.1), mar3=c(0,1,0.1), dep=1.5, A1=0, A2=0, B1=3, B2=3, checkconv=TRUE, ...)
{
models            = c("log","psilog","neglog","psineglog")
if(!(model %in% models)) stop(paste("'",model,"' is not a valid model.",sep="")) else{
if(model=="log")        ptgpd  = ptgpd_log      (x,y,z,mar1=mar1,mar2=mar2,mar3=mar3,dep=dep)
if(model=="psilog")     ptgpd  = ptgpd_psilog   (x,y,z,mar1=mar1,mar2=mar2,mar3=mar3,dep=dep,A1=A1,A2=A2,B1=B1,B2=B2,checkconv=checkconv)
if(model=="neglog")     ptgpd  = ptgpd_neglog   (x,y,z,mar1=mar1,mar2=mar2,mar3=mar3,dep=dep)
if(model=="psineglog")  ptgpd  = ptgpd_psineglog(x,y,z,mar1=mar1,mar2=mar2,mar3=mar3,dep=dep,A1=A1,A2=A2,B1=B1,B2=B2,checkconv=checkconv)
ptgpd
}}
