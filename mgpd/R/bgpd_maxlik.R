bgpd_maxlik <-
function(param,dat,model="log",fixed=FALSE,psi=3,...)
{
#print(param)
mlmax=1e+15
models      = c("log","psilog","philog",
"neglog","psineglog","phineglog",
"bilog","negbilog","ct","taj")#,"smith","mix")
if(!(model %in% models)) stop(paste("'",model,"' is not a valid model.",sep="")) else{
if(model=="log")        ml  = ml_log(        param=param,dat=dat,mlmax=mlmax,fixed=fixed)
#if(model=="smith")      ml  = ml_smith(      param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="psilog")     ml  = ml_psilog(     param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="philog")     ml  = ml_philog(     param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="neglog")     ml  = ml_neglog(     param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="psineglog")  ml  = ml_psineglog(  param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="phineglog")  ml  = ml_phineglog(  param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="bilog")      ml  = ml_bilog(      param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="negbilog")   ml  = ml_negbilog(   param=param,dat=dat,mlmax=mlmax,fixed=fixed)
#if(model=="mix")        ml  = ml_mix(        param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="ct")         ml  = ml_ct(         param=param,dat=dat,mlmax=mlmax,fixed=fixed)
if(model=="taj")        ml  = ml_taj(        param=param,dat=dat,mlmax=mlmax,fixed=fixed)
ml
}
}
