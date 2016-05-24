fbgpd <-
function(initpar,dat,model="log",fixed=FALSE,control=list(maxit=50000),psi=3,...)
{
est         = optim(initpar,bgpd_maxlik,dat=dat,model=model,control=control,fixed=fixed,psi=psi)
if(fixed==TRUE) est$par[1]=0
est
}
