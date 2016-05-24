pbgpd_taj <-
function(x,y,mar1=c(0,1,0.1),mar2=c(0,1,0.1),a=2,b=1.5,...)
{
mu           = function(x,y) ((1/x)^(2*a)+2*(1+b)*(1/x/y)^(a)+(1/y)^(2*a))^(1/2/a)

param        = as.numeric(c(mar1,mar2,a,b))
mux          = param[1]; muy    = param[4]
sigx         = param[2]; sigy   = param[5]
gamx         = param[3]; gamy   = param[6]
a            = param[7]; b      = param[8]

Hxy               = NULL
error             = FALSE 
if(sigx<0 | sigy<0 | a<1 | b<=-1 | (b>(2*a-2))) error=TRUE
if(!error){

Hxy          = NA
tx           = (1+gamx*(x-mux)/sigx)^(1/gamx)
ty           = (1+gamy*(y-muy)/sigy)^(1/gamy)
tx0          = (1+gamx*(-mux)/sigx)^(1/gamx)
ty0          = (1+gamy*(-muy)/sigy)^(1/gamy)
c0           = -mu(tx0,ty0)
Hxy          = 1/c0*(mu(tx,ty)-mu(pmin(tx,tx0),pmin(ty,ty0)))
}else stop("invalid parameter(s)")
Hxy
}
