model.sel.distbn <-
function(delta, x, min.branch, mod.id){
	num.par <- c(1,2,2,2)
	
	eaic <- waic <- laic <- vaic <- exp.out <- weib.out <- lnorm.out <- vrat.out <- NA
	
	if(mod.id[1]){
		exp.out <- censor.exp.x(delta,x,min.branch)
		eaic <- (-2*(exp.out$LL - num.par[1]))+(2*num.par[1]*(num.par[1]+1))/(length(x)-num.par[1]-1)
		}
	if(mod.id[2]){
		weib.out <- censor.weib.x(delta,x,min.branch)
		waic <- -2*(weib.out$LL - num.par[2])+(2*num.par[2]*(num.par[2]+1))/(length(x)-num.par[2]-1)
		}
	if(mod.id[3]){
		lnorm.out <- censor.lnorm.x(delta,x,min.branch)
		laic <- -2*(lnorm.out$LL - num.par[3])+(2*num.par[3]*(num.par[3]+1))/(length(x)-num.par[3]-1)
		}
	if(mod.id[4]){
		vrat.out <- censor.vrat.x(delta,x,min.branch)
		vaic <- -2*(vrat.out$LL - num.par[4])+(2*num.par[4]*(num.par[4]+1))/(length(x)-num.par[4]-1)
		}
	out<-list(exp.out,weib.out,lnorm.out,vrat.out)
	AICs <- c(eaic,waic,laic,vaic)
	AICs <- replace(AICs,mod.id==FALSE,NA)
	best.model.indx <- which(AICs==min(AICs,na.rm=TRUE))
	
	a<-data.frame(out[[best.model.indx]],AICs[best.model.indx],best.model.indx,num.par[best.model.indx])
	names(a) <- c("P1","P2","LL","AICc","Mod","n")
	return(a)
	}

