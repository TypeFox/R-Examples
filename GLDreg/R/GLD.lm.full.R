GLD.lm.full<-
function (formula, data, param, maxit = 20000, fun, method = "Nelder-Mead", 
range = c(0.01, 0.99),n.simu=1000, summary.plot=TRUE, init=NULL) 
{
fit<-GLD.lm(formula, data, param, maxit, fun, method, diagnostics=FALSE, range,
init) 

formula.info<-unlist(strsplit(deparse(formula)," ~ "))

if(length(formula.info)==2 & formula.info[[2]]=="."){

data<-na.omit(data)

fit.simu.f<-fun.simu.gld.lm.alt(n.simu,
fit$y+rgl(nrow(fit$Fitted),
fit$'Estimated parameters'[match(c("L1","L2","L3","L4"),
names(fit$'Estimated parameters'))],param=fit$param)~.,
fit,data=data[,-c(match(formula.info[1],dimnames(data)[[2]]))],param=fit$param,
fun=fun,init=init)

}

if(length(formula.info)!=2 | formula.info[[2]]!="."){


fit.simu.f<-fun.simu.gld.lm.alt(n.simu,
update(fit$formula,fit$y+rgl(nrow(fit$Fitted),
fit$'Estimated parameters'[match(c("L1","L2","L3","L4"),
names(fit$'Estimated parameters'))],param=fit$param)~.),
fit,data=na.omit(subset(data,select=all.vars(fit$formula))),
param=fit$param,fun=fun,init=init)

}


fit.bc.f<-fun.simu.bias.correct.alt(fit.simu.f,fit)

result<-list(fit,fit.simu.f,fit.bc.f)
if(summary.plot==TRUE){
summaryGraphics.gld.lm(result)}
return(result)
}

 