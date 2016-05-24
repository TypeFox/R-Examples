svycox.validate <-
function(.boot.index,.nom,.data){

internal.validate.func=function(.boot.vec,.nom2,.data2){
boot.data=.data2[.boot.vec,]
design.boot=svydesign(ids=~1,strata=as.formula(paste("~",names(.nom2$design$strata))),
							probs=as.formula(paste("~",names(.nom2$design$prob))),
							fpc=as.formula(paste("~",colnames(.nom2$design$fpc$popsize))),				
						data=boot.data)
boot.fit<-svycoxph(formula(.nom2$svy.cox), x=TRUE,design=design.boot)
lp.boot=boot.fit$x%*%as.matrix(boot.fit$coefficients)-
		mean(boot.fit$x%*%as.matrix(boot.fit$coefficients))
lp.test=(.nom2$svy.cox)$x%*%as.matrix(boot.fit$coefficients)-
		mean((.nom2$svy.cox)$x%*%as.matrix(boot.fit$coefficients))
cindex.train<-1-rcorr.cens(lp.boot,Surv(boot.data$survival,boot.data$surv_cens))[[1]]
cindex.test<-1-rcorr.cens(lp.test,Surv(.data2$survival,.data2$surv_cens))[[1]]	
return(cindex.train-cindex.test)
}

val.res=apply(.boot.index,2,internal.validate.func,.nom2=.nom,.data2=.data)
print(mean(val.res))
return(val.res)

}

