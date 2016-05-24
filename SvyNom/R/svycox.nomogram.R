svycox.nomogram <-
function(.design,.model,.data,pred.at,fun.lab){
design.call=.design$call
svy.cox.fit=svycoxph(.model,x=TRUE,design=.design)
pred.lp.cox=predict(svy.cox.fit)
pred.survey.cox=predict(svy.cox.fit,type="curve",newdata=.data)
.rhs=.model[[3]]
f.form=paste("pred.lp.cox~",paste(all.vars(.model)[-(1:2)],collapse="+"))
.f=ols(as.formula(f.form),sigma=1,x=TRUE,y=TRUE,data=.data)
.ss3<-c(0.05,0.2,0.4,0.6,0.7,0.8,0.9,0.95,0.99)
.ss3.label<-100*.ss3
time.at=pred.survey.cox[[1]]$time[which(pred.survey.cox[[1]]$time>pred.at)[1]-1]
.baseline=exp(log(pred.survey.cox[[1]]$surv[names(time.at)])/exp(svy.cox.fit$linear.predictors[1]))

.tempfun=function(x) .baseline[[1]]^exp(x)

.nom=nomogram(.f, fun=.tempfun, funlabel=fun.lab,fun.at=.ss3,lp=T, vnames="labels")
return(list(nomog=.nom,design=.design,svy.cox=svy.cox.fit,preds=pred.survey.cox,
		pred.at=pred.at))
}

