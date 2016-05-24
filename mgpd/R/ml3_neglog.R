ml3_neglog <-
function(param,dat,mlmax=1e+15,fixed=FALSE,...)
{
#param=c(0,1,0.1,0,1,0.1,0,1,0.1,1.5)
#dat=potdata[,c(1:3)]
#fixed=FALSE
#mlmax=1e+15

loglik        = mlmax
lik           = NULL
x             = dat[,1]
y             = dat[,2]
z             = dat[,3]

if(fixed)     param[1]=0

lik           = try(dtgpd_neglog(x, y, z, mar1 = param[1:3],mar2 = param[4:6],mar3 = param[7:9],dep  = param[10]))
#points(lik,col=2)
#print(summary(lik-dens))
#print(-sum(log(lik)))

if(!is.null(lik)){
    loglik    = -sum(log(lik))
    if(min(1+param[3]*(x-param[1])/param[2])<0) loglik=mlmax
    if(min(1+param[6]*(y-param[4])/param[5])<0) loglik=mlmax
    if(min(1+param[9]*(z-param[7])/param[8])<0) loglik=mlmax}
loglik
}
