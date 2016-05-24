ml_psilog <-
function(param,dat,mlmax=1e+15,fixed=FALSE,...)
{

loglik        = mlmax
lik           = NULL
x             = dat[,1]
y             = dat[,2]

if(fixed)     param[1]=0

lik           = try(dbgpd(x, y, model = "psilog",mar1 = param[1:3],mar2 = param[4:6],dep  = param[7],asy=param[8], p=param[9]),
silent=TRUE)

if(!is.null(lik)){
if(is.null(attr(lik,"class"))){
    loglik    = -sum(log(lik))
    if(min(1+param[3]*(x-param[1])/param[2])<0) loglik=mlmax
    if(min(1+param[6]*(y-param[4])/param[5])<0) loglik=mlmax
    }}

loglik
}
