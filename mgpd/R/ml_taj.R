ml_taj <-
function(param,dat,mlmax=1e+15,fixed=FALSE,...)
{

loglik       = mlmax
hxy          = NA
x            = dat[,1]
y            = dat[,2]
error        = FALSE

mux          = param[1]; muy    = param[4]
sigx         = param[2]; sigy   = param[5]
gamx         = param[3]; gamy   = param[6]
a            = param[7]; b      = param[8]
mu           = function(x,y) ((1/x)^(2*a)+2*(1+b)*(1/x/y)^(a)+(1/y)^(2*a))^(1/2/a)
dxdymu       = function(x,y) (1/4)*((1/x)^(2*a)+(2*(1+b))*(1/(x*y))^a+(1/y)^(2*a))^
(1/(2*a))*(-(2*(1+b))*(1/(x*y))^a*a/y-2*(1/y)^(2*a)*a/y)*(-2*(1/x)^(2*a)*a/x-(2*(1+b))*
(1/(x*y))^a*a/x)/(a^2*((1/x)^(2*a)+(2*(1+b))*(1/(x*y))^a+(1/y)^(2*a))^2)+((1/x)^(2*a)+
(2*(1+b))*(1/(x*y))^a+(1/y)^(2*a))^(1/(2*a))*a*(1+b)*(1/(x*y))^a/(y*x*((1/x)^(2*a)+
(2*(1+b))*(1/(x*y))^a+(1/y)^(2*a)))-(1/2)*((1/x)^(2*a)+(2*(1+b))*(1/(x*y))^a+(1/y)^
(2*a))^(1/(2*a))*(-2*(1/x)^(2*a)*a/x-(2*(1+b))*(1/(x*y))^a*a/x)*(-(2*(1+b))*(1/(x*y))^
a*a/y-2*(1/y)^(2*a)*a/y)/(a*((1/x)^(2*a)+(2*(1+b))*(1/(x*y))^a+(1/y)^(2*a))^2)

if(sigx<0 | sigy<0 | a<1 | b<=-1 | (b>(2*a-2))) error = TRUE
if(fixed==TRUE) {mux=0}
if(error) loglik = mlmax
if(!error)
{
tx           = (1+gamx*(x-mux)/sigx)^(1/gamx)
ty           = (1+gamy*(y-muy)/sigy)^(1/gamy)
tx0          = (1+gamx*(-mux)/sigx)^(1/gamx)
ty0          = (1+gamy*(-muy)/sigy)^(1/gamy)
dtx          = (1/sigx)*pmax((1+gamx*(x-mux)/sigx),0)^(1/gamx-1)
dty          = (1/sigy)*pmax((1+gamy*(y-muy)/sigy),0)^(1/gamy-1)
c0           = -mu(tx0,ty0)
hxy          = 1/c0*dxdymu(tx,ty)*dtx*dty
hxy          = as.numeric(hxy*(1-((x<0)*(y<0))))
loglik       = -sum(log(hxy))
}
if(min(1+gamx*(x-mux)/sigx)<0) loglik=mlmax
if(min(1+gamy*(y-muy)/sigy)<0) loglik=mlmax
loglik
}
