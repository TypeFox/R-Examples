
surv.km.base = function(tl, dl, tt, weight.perturb=NULL)	{
	if(is.null(weight.perturb)){	weight.perturb = rep(1, length(tl))}
	S.KM = survfit(Surv(tl,dl)~1, weights = weight.perturb)	
	S.tt.KM = approx(S.KM$time,S.KM$surv,tt)$y
	return(S.tt.KM)
}

surv.km = function(tl, dl, tt, var = FALSE, conf.int = FALSE, weight.perturb=NULL, perturb.vector=FALSE)	{
	S.KM.est = surv.km.base(tl=tl, dl=dl, tt=tt)
	if(var | conf.int)	{
		if(is.null(weight.perturb)){	
			weight.perturb = matrix(rexp(500*length(tl), rate=1), ncol = 500)
		}
		S.KM.p = apply(weight.perturb, 2, surv.km.base, tl=tl, dl=dl, tt = tt)
		if(conf.int)	{
			conf.l.normal.S = S.KM.est - 1.96*sd(S.KM.p)
			conf.u.normal.S = S.KM.est + 1.96*sd(S.KM.p)
			conf.l.quantile.S = quantile(S.KM.p, 0.025)
			conf.u.quantile.S = quantile(S.KM.p, 0.975)
		}
	}
	if(!(var) & !(conf.int) )	{
		return(list("S.estimate" = S.KM.est))
	}
	if(var & !(conf.int) & !perturb.vector)	{
		return(list("S.estimate" = S.KM.est, "S.var" = var(S.KM.p)))
	}
	if(conf.int & !perturb.vector)	{
		return(list("S.estimate" = S.KM.est, "S.var" = var(S.KM.p), "conf.int.normal.S" = as.vector(c(conf.l.normal.S, conf.u.normal.S)), "conf.int.quantile.S" = as.vector(as.vector(c(conf.l.quantile.S, conf.u.quantile.S)))))
	}	
	if(var & !(conf.int) & perturb.vector)	{
		return(list("S.estimate" = S.KM.est, "S.var" = var(S.KM.p), "perturb.vector" = S.KM.p))
	}
	if(conf.int & perturb.vector)	{
		return(list("S.estimate" = S.KM.est, "S.var" = var(S.KM.p), "conf.int.normal.S" = as.vector(c(conf.l.normal.S, conf.u.normal.S)), "conf.int.quantile.S" = as.vector(as.vector(c(conf.l.quantile.S, conf.u.quantile.S))),"perturb.vector" = S.KM.p))
	}	

}

delta.km = function(tl, dl, treat, tt, var = FALSE, conf.int = FALSE, weight.perturb=NULL)	{
	if(sum(unique(treat) %in% c(0,1)) != 2){
		treat.f = as.factor(treat)
		unique.treat = unique(treat.f)[1]
		treat = 1*(treat.f==unique.treat)
		print(paste("Treatment variable was not 0/1; converted to 0/1 by 1=", unique.treat,", 0=", unique(treat.f)[2], sep =""))
	}
	if(is.null(weight.perturb)){
		weightperturb.1 = NULL
		weightperturb.0 = NULL
	}
	if(!is.null(weight.perturb)){
		weightperturb.1 = weight.perturb[treat==1,]
		weightperturb.0 = weight.perturb[treat==0,]
	}
	surv.1 = surv.km(tl=tl[treat==1], dl = dl[treat==1], tt=tt, var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.1) 
	surv.0 = surv.km(tl=tl[treat==0], dl = dl[treat==0], tt=tt, var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.0)
	delta = surv.1$S.estimate - surv.0$S.estimate
	if(var | conf.int){
		delta.p =  surv.1$perturb.vector - surv.0$perturb.vector
		delta.var = var(delta.p)
	}
	if(conf.int)	{
		conf.l.normal.d = delta - 1.96*sd(delta.p)
		conf.u.normal.d = delta + 1.96*sd(delta.p)
		conf.l.quantile.d = quantile(delta.p, 0.025)
		conf.u.quantile.d = quantile(delta.p, 0.975)
	}
	if(var | conf.int){
		pval = pnorm(delta/sqrt(delta.var), lower.tail = (delta/sqrt(delta.var) < 0))*2
	}
	if(!var & !conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate, "delta.estimate" = delta))}
	if(var & !conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate,"delta.estimate" = delta, "S.var.1" = surv.1$S.var, "S.var.0" = surv.0$S.var, "delta.var" = delta.var, "p.value" = pval))}
	if(conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate, "delta.estimate" = delta, "S.var.1" = surv.1$S.var, "S.var.0" = surv.0$S.var, "delta.var" = delta.var,  "conf.int.normal.S.1" = surv.1$conf.int.normal.S,  "conf.int.normal.S.0" = surv.0$conf.int.normal.S, "conf.int.normal.delta" = c(conf.l.normal.d, conf.u.normal.d), "conf.int.quantile.S.1" = surv.1$conf.int.quantile.S,  "conf.int.quantile.S.0" = surv.0$conf.int.quantile.S,  "conf.int.quantile.delta" = as.vector(c(conf.l.quantile.d, conf.u.quantile.d)), "p.value" = pval))}
}

surv.iptw.km.base = function(tl, dl, tt, ps.weights,weight.perturb=NULL, perturb.ps = FALSE)	{
    if(is.null(weight.perturb)) {
    	S.KM = survfit(Surv(tl,dl)~1, weights = ps.weights)}
    if(!is.null(weight.perturb) & perturb.ps == FALSE) {
    	S.KM = survfit(Surv(tl,dl)~1, weights = ps.weights*weight.perturb)}
    if(!is.null(weight.perturb) & perturb.ps == TRUE) {
    	S.KM = survfit(Surv(tl,dl)~1, weights = weight.perturb)}
	S.tt.KM = approx(S.KM$time,S.KM$surv,tt)$y
	return(S.tt.KM)
}

surv.iptw.km = function(tl, dl, tt, var = FALSE, conf.int = FALSE, ps.weights,weight.perturb=NULL, perturb.ps=FALSE, perturb.vector=FALSE)	{
	S.KM.est = surv.iptw.km.base(tl=tl, dl=dl, tt=tt, ps.weights = ps.weights)
	if(var | conf.int)	{
		if(is.null(weight.perturb)){	
			weight.perturb = matrix(rexp(500*length(tl), rate=1), ncol = 500)
			S.KM.p = apply(weight.perturb, 2, surv.iptw.km.base, tl=tl, dl=dl, tt = tt,ps.weights = ps.weights, perturb.ps = FALSE)
		}
		if(!is.null(weight.perturb)){	
			S.KM.p = apply(weight.perturb, 2, surv.iptw.km.base, tl=tl, dl=dl, tt = tt,ps.weights = ps.weights, perturb.ps = perturb.ps)
		}		
		if(conf.int)	{
			conf.l.normal.S = S.KM.est - 1.96*sd(S.KM.p)
			conf.u.normal.S = S.KM.est + 1.96*sd(S.KM.p)
			conf.l.quantile.S = quantile(S.KM.p, 0.025)
			conf.u.quantile.S = quantile(S.KM.p, 0.975)
		}
	}
	if(!(var) & !(conf.int))	{
		return(list("S.estimate" = S.KM.est))
	}
	if(var & !(conf.int) & !perturb.vector)	{
		return(list("S.estimate" = S.KM.est, "S.var" = var(S.KM.p)))
	}
	if(conf.int & !perturb.vector)	{
		return(list("S.estimate" = S.KM.est, "S.var" = var(S.KM.p), "conf.int.normal.S" = as.vector(c(conf.l.normal.S, conf.u.normal.S)), "conf.int.quantile.S" = as.vector(as.vector(c(conf.l.quantile.S, conf.u.quantile.S))), "perturb.vector" = S.KM.p))
	}	
	if(var & !(conf.int) & perturb.vector)	{
		return(list("S.estimate" = S.KM.est, "S.var" = var(S.KM.p)))
	}
	if(conf.int & perturb.vector)	{
		return(list("S.estimate" = S.KM.est, "S.var" = var(S.KM.p), "conf.int.normal.S" = as.vector(c(conf.l.normal.S, conf.u.normal.S)), "conf.int.quantile.S" = as.vector(as.vector(c(conf.l.quantile.S, conf.u.quantile.S))), "perturb.vector" = S.KM.p))
	}
}


delta.iptw.km = function(tl, dl, treat, tt, var = FALSE, conf.int = FALSE,  ps.weights=NULL, weight.perturb=NULL, perturb.ps=FALSE,cov.for.ps = NULL)	{
	if(sum(unique(treat) %in% c(0,1)) != 2){
		treat.f = as.factor(treat)
		unique.treat = unique(treat.f)[1]
		treat = 1*(treat.f==unique.treat)
		print(paste("Treatment variable was not 0/1; converted to 0/1 by 1=", unique.treat,", 0=", unique(treat.f)[2], sep =""))
	}
	if(!var & !conf.int)	{
		weightperturb.1 = NULL
		weightperturb.0 = NULL
	}
	if(var | conf.int)	{
		if(is.null(weight.perturb)){
			weight.perturb = matrix(rexp(500*length(tl), rate=1), ncol = 500)
			weightperturb.1 = weight.perturb[treat==1,]
			weightperturb.0 = weight.perturb[treat==0,]
		}
		if(!is.null(weight.perturb)){
		weightperturb.1 = weight.perturb[treat==1,]
		weightperturb.0 = weight.perturb[treat==0,]
		}
	}	
	if(is.null(ps.weights))	{
		if(is.null(cov.for.ps))	{stop("must either supply propensity score weights or supply covariate information for the function to construct weights")}
		ps.weights = ps.wgt.fun(treat=treat, cov.for.ps = cov.for.ps)
		if(var | conf.int){
			ps.weights.perturb = apply(weight.perturb, 2, ps.wgt.fun, treat=treat, cov.for.ps = cov.for.ps)
			weightperturb.1 = weight.perturb[treat==1,]*ps.weights.perturb[treat==1,]
			weightperturb.0 = weight.perturb[treat==0,]*ps.weights.perturb[treat==0,]
			perturb.ps = TRUE
		}
		if(!var & !conf.int)	{
			perturb.ps=FALSE
		}	
	}
	surv.1 = surv.iptw.km(tl=tl[treat==1], dl = dl[treat==1], tt=tt, var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.1, ps.weights=ps.weights[treat==1], perturb.ps=perturb.ps) 
	surv.0 = surv.iptw.km(tl=tl[treat==0], dl = dl[treat==0], tt=tt, var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.0, ps.weights=ps.weights[treat==0], perturb.ps=perturb.ps) 
	delta = surv.1$S.estimate - surv.0$S.estimate
	if(var | conf.int){
		delta.p =  surv.1$perturb.vector - surv.0$perturb.vector
		delta.var = var(delta.p)
	}
	if(conf.int)	{
		conf.l.normal.d = delta - 1.96*sd(delta.p)
		conf.u.normal.d = delta + 1.96*sd(delta.p)
		conf.l.quantile.d = quantile(delta.p, 0.025)
		conf.u.quantile.d = quantile(delta.p, 0.975)
	}
	if(var | conf.int){
		pval = pnorm(delta/sqrt(delta.var), lower.tail = (delta/sqrt(delta.var) < 0))*2
	}
	if(!var & !conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate, "delta.estimate" = delta))}
	if(var & !conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate,"delta.estimate" = delta, "S.var.1" = surv.1$S.var, "S.var.0" = surv.0$S.var, "delta.var" = delta.var, "p.value" = pval))}
	if(conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate, "delta.estimate" = delta, "S.var.1" = surv.1$S.var, "S.var.0" = surv.0$S.var, "delta.var" = delta.var,  "conf.int.normal.S.1" = surv.1$conf.int.normal.S,  "conf.int.normal.S.0" = surv.0$conf.int.normal.S, "conf.int.normal.delta" = c(conf.l.normal.d, conf.u.normal.d), "conf.int.quantile.S.1" = surv.1$conf.int.quantile.S,  "conf.int.quantile.S.0" = surv.0$conf.int.quantile.S,  "conf.int.quantile.delta" = as.vector(c(conf.l.quantile.d, conf.u.quantile.d)), "p.value" = pval))}
}

surv.land.rct.base = function(tl, dl, tt, landmark, short=NULL, z.cov=NULL, weight.perturb=NULL, bw=NULL){
	if(is.null(weight.perturb)) {
    	weight.perturb = rep(1,length(tl))}
    if(is.null(z.cov) & is.null(short)) {stop("must supply either baseline covariate information or intermediate event information, otherwise you could just use regular Kaplan Meier.")}
	if(!is.null(z.cov) & is.null(short)) {
		betahat.cox.t <- coxph(Surv(tl, dl)~z.cov,  weights = weight.perturb)  
		score.t.i = c(as.matrix(z.cov)%*%betahat.cox.t$coef)
		Surv.score.t.i = Est.KM.FUN.weighted(xi=tl, di=dl, si=score.t.i, myt=tt, weight.perturb = weight.perturb, bw=bw)
		S.t.rct = sum(Surv.score.t.i*weight.perturb)/sum(weight.perturb)
	}
	if(!is.null(short)){
		d = dim(short)[2]
		if(d %% 2 > 0) {stop("the number of columns for the intermediate event information matrix needs to be a multiple of 2.")}
		short.new = matrix(nrow = length(tl[tl>landmark]), ncol = d)
		num.s = (d %/% 2)
		for(q in 1:num.s){
			short.new[,(2*q-1)] = 1*(short[tl>landmark,(2*q-1) ] <= landmark)
			short.new[,(2*q)] = pmin(short[tl>landmark,(2*q-1) ],landmark)	
	}
	}	
	if(is.null(z.cov) & !is.null(short)) {
		S.t0.propose = surv.km.base(tl=tl, dl=dl, tt=landmark, weight.perturb=weight.perturb)
		new.mat = short.new
		betahat.cox.t <- coxph(Surv(tl[tl>landmark], dl[tl>landmark])~new.mat, weights = weight.perturb[tl>landmark])  
		score.t.i = c(as.matrix(new.mat)%*%betahat.cox.t$coef)
		Surv.score.t.i = Est.KM.FUN.weighted(xi=tl[tl>landmark], di=dl[tl>landmark], si=score.t.i, myt=tt, weight.perturb = weight.perturb[tl>landmark], bw=bw)
		S.t.rct = (sum(Surv.score.t.i*weight.perturb[tl>landmark])/sum(weight.perturb[tl>landmark]))*S.t0.propose
	}
	if(!is.null(z.cov) & !is.null(short)) {
		betahat.cox.t0 <- coxph(Surv(tl, dl)~z.cov, weights = weight.perturb)  
		score.t0.i = c(as.matrix(z.cov)%*%betahat.cox.t0$coef)
		Surv.score.t0.i = Est.KM.FUN.weighted(xi=tl, di=dl,si=score.t0.i,myt=landmark, weight.perturb = weight.perturb, bw=bw)
		S.t0.propose = sum(Surv.score.t0.i*weight.perturb)/sum(weight.perturb)	
		new.mat = cbind(short.new, as.matrix(z.cov)[tl>landmark,])
		betahat.cox.t <- coxph(Surv(tl[tl>landmark], dl[tl>landmark])~new.mat, weights = weight.perturb[tl>landmark])  
		score.t.i = c(as.matrix(new.mat)%*%betahat.cox.t$coef)
		Surv.score.t.i = Est.KM.FUN.weighted(xi=tl[tl>landmark], di=dl[tl>landmark], si=score.t.i, myt=tt, weight.perturb = weight.perturb[tl>landmark], bw=bw)
		S.t.rct = (sum(Surv.score.t.i*weight.perturb[tl>landmark])/sum(weight.perturb[tl>landmark]))*S.t0.propose
		}
	return(S.t.rct)
}

surv.land.rct = function(tl, dl, tt, landmark, short=NULL, z.cov=NULL, var = FALSE, conf.int = FALSE, weight.perturb=NULL, perturb.vector=FALSE, bw = NULL)	{
	S.est = surv.land.rct.base(tl=tl, dl=dl, tt=tt, landmark=landmark, short=short, z.cov=z.cov)
	if(var | conf.int)	{
		if(is.null(weight.perturb)){	
			weight.perturb = matrix(rexp(500*length(tl), rate=1), ncol = 500)
		}
		S.p = apply(weight.perturb, 2, surv.land.rct.base, tl=tl, dl=dl, tt = tt, landmark=landmark, short=short, z.cov=z.cov)
		if(conf.int)	{
			conf.l.normal.S = S.est - 1.96*sd(S.p)
			conf.u.normal.S = S.est + 1.96*sd(S.p)
			conf.l.quantile.S = quantile(S.p, 0.025)
			conf.u.quantile.S = quantile(S.p, 0.975)
		}
	}
	if(!(var) & !(conf.int))	{
		return(list("S.estimate" = S.est))
	}
	if(var & !(conf.int) & !perturb.vector)	{
		return(list("S.estimate" = S.est, "S.var" = var(S.p)))
	}
	if(conf.int & !perturb.vector)	{
		return(list("S.estimate" = S.est, "S.var" = var(S.p), "conf.int.normal.S" = as.vector(c(conf.l.normal.S, conf.u.normal.S)), "conf.int.quantile.S" = as.vector(as.vector(c(conf.l.quantile.S, conf.u.quantile.S)))))
	}	
	if(var & !(conf.int) & perturb.vector)	{
		return(list("S.estimate" = S.est, "S.var" = var(S.p), "perturb.vector" = S.p))
	}
	if(conf.int & perturb.vector)	{
		return(list("S.estimate" = S.est, "S.var" = var(S.p), "conf.int.normal.S" = as.vector(c(conf.l.normal.S, conf.u.normal.S)), "conf.int.quantile.S" = as.vector(as.vector(c(conf.l.quantile.S, conf.u.quantile.S))), "perturb.vector" = S.p))
	}	

}

delta.land.rct = function(tl, dl, treat, tt, landmark, short=NULL, z.cov=NULL, var = FALSE, conf.int = FALSE, weight.perturb=NULL, bw = NULL)	{
	if(sum(unique(treat) %in% c(0,1)) != 2){
		treat.f = as.factor(treat)
		unique.treat = unique(treat.f)[1]
		treat = 1*(treat.f==unique.treat)
		print(paste("Treatment variable was not 0/1; converted to 0/1 by 1=", unique.treat,", 0=", unique(treat.f)[2], sep =""))
	}
	if(is.null(weight.perturb)){
		weight.perturb = matrix(rexp(500*length(tl), rate=1), ncol = 500)
		weightperturb.1 = weight.perturb[treat==1,]
		weightperturb.0 = weight.perturb[treat==0,]
	}
	if(!is.null(weight.perturb)){
		weightperturb.1 = weight.perturb[treat==1,]
		weightperturb.0 = weight.perturb[treat==0,]
	}
	if(is.null(z.cov) & is.null(short)) {stop("must supply either baseline covariate information or intermediate event information, otherwise you could just use regular Kaplan Meier.")}
	if(!is.null(short)){
		d = dim(short)[2]
		if(d %% 2 > 0) {stop("the number of columns for the intermediate event information matrix needs to be a multiple of 2.")}
	}	
	if(!is.null(z.cov) & !is.null(short)) {
		surv.1 = surv.land.rct(tl=tl[treat==1], dl = dl[treat==1], tt=tt, landmark = landmark, short = short[treat==1,], z.cov = as.matrix(z.cov)[treat==1,], var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.1) 
		surv.0 = surv.land.rct(tl=tl[treat==0], dl = dl[treat==0], tt=tt, landmark = landmark, short = short[treat==0,], z.cov = as.matrix(z.cov)[treat==0,], var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.0) 
	}
	if(is.null(z.cov)) {
		surv.1 = surv.land.rct(tl=tl[treat==1], dl = dl[treat==1], tt=tt, landmark = landmark, short = short[treat==1,], z.cov = NULL, var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.1) 
		surv.0 = surv.land.rct(tl=tl[treat==0], dl = dl[treat==0], tt=tt, landmark = landmark, short = short[treat==0,], z.cov = NULL, var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.0) 
	}
	if(is.null(short)) {
		surv.1 = surv.land.rct(tl=tl[treat==1], dl = dl[treat==1], tt=tt, landmark = landmark, short = NULL, z.cov = as.matrix(z.cov)[treat==1,], var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.1) 
		surv.0 = surv.land.rct(tl=tl[treat==0], dl = dl[treat==0], tt=tt, landmark = landmark, short =  NULL, z.cov = as.matrix(z.cov)[treat==0,], var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.0) 
	}
	delta = surv.1$S.estimate - surv.0$S.estimate
	if(var | conf.int){
		delta.p =  surv.1$perturb.vector - surv.0$perturb.vector
		delta.var = var(delta.p)
	}
	if(conf.int)	{
		conf.l.normal.d = delta - 1.96*sd(delta.p)
		conf.u.normal.d = delta + 1.96*sd(delta.p)
		conf.l.quantile.d = quantile(delta.p, 0.025)
		conf.u.quantile.d = quantile(delta.p, 0.975)
	}
	if(var | conf.int){
		pval = pnorm(delta/sqrt(delta.var), lower.tail = (delta/sqrt(delta.var) < 0))*2
	}
	if(!var & !conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate, "delta.estimate" = delta))}
	if(var & !conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate,"delta.estimate" = delta, "S.var.1" = surv.1$S.var, "S.var.0" = surv.0$S.var, "delta.var" = delta.var, "p.value" = pval))}
	if(conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate, "delta.estimate" = delta, "S.var.1" = surv.1$S.var, "S.var.0" = surv.0$S.var, "delta.var" = delta.var,  "conf.int.normal.S.1" = surv.1$conf.int.normal.S,  "conf.int.normal.S.0" = surv.0$conf.int.normal.S, "conf.int.normal.delta" = c(conf.l.normal.d, conf.u.normal.d), "conf.int.quantile.S.1" = surv.1$conf.int.quantile.S,  "conf.int.quantile.S.0" = surv.0$conf.int.quantile.S,  "conf.int.quantile.delta" = as.vector(c(conf.l.quantile.d, conf.u.quantile.d)), "p.value" = pval))}
}


surv.land.obs.base = function(tl, dl, tt, landmark, short=NULL, z.cov=NULL, ps.weights, weight.perturb=NULL, perturb.ps=FALSE, bw=NULL){
	if(is.null(weight.perturb)) {
    	new.weights = ps.weights}
    if(!is.null(weight.perturb) & perturb.ps == FALSE) {
    	new.weights = ps.weights*weight.perturb}
	if(!is.null(weight.perturb) & perturb.ps == TRUE) {
    	new.weights = weight.perturb}
    if(is.null(z.cov) & is.null(short)) {stop("must supply either baseline covariate information or intermediate event information, otherwise you could just use regular Kaplan Meier.")}
	if(!is.null(z.cov) & is.null(short)) {
		betahat.cox.t <- coxph(Surv(tl, dl)~z.cov,  weights = new.weights)  
		score.t.i = c(as.matrix(z.cov)%*%betahat.cox.t$coef)
		Surv.score.t.i = Est.KM.FUN.weighted(xi=tl, di=dl, si=score.t.i, myt=tt, weight.perturb = new.weights, bw=bw)
		S.t.obs = sum(Surv.score.t.i*new.weights)/sum(new.weights)
	}
	if(!is.null(short)){
		d = dim(short)[2]
		if(d %% 2 > 0) {stop("the number of columns for the intermediate event information matrix needs to be a multiple of 2.")}
		short.new = matrix(nrow = length(tl[tl>landmark]), ncol = d)
		num.s = (d %/% 2)
		for(q in 1:num.s){
			short.new[,(2*q-1)] = 1*(short[tl>landmark,(2*q-1) ] <= landmark)
			short.new[,(2*q)] = pmin(short[tl>landmark,(2*q-1) ],landmark)	
	}
	}	
	if(is.null(z.cov) & !is.null(short)) {
		S.t0.propose = surv.km.base(tl=tl, dl=dl, tt=landmark, weight.perturb=new.weights)
		new.mat = short.new
		betahat.cox.t <- coxph(Surv(tl[tl>landmark], dl[tl>landmark])~new.mat, weights = new.weights[tl>landmark])  
		score.t.i = c(as.matrix(new.mat)%*%betahat.cox.t$coef)
		Surv.score.t.i = Est.KM.FUN.weighted(xi=tl[tl>landmark], di=dl[tl>landmark], si=score.t.i, myt=tt, weight.perturb = new.weights[tl>landmark], bw=bw)
		S.t.obs = (sum(Surv.score.t.i*new.weights[tl>landmark])/sum(new.weights[tl>landmark]))*S.t0.propose
	}
	if(!is.null(z.cov) & !is.null(short)) {
		betahat.cox.t0 <- coxph(Surv(tl, dl)~z.cov, weights = new.weights)  
		score.t0.i = c(as.matrix(z.cov)%*%betahat.cox.t0$coef)
		Surv.score.t0.i = Est.KM.FUN.weighted(xi=tl, di=dl,si=score.t0.i,myt=landmark, weight.perturb = new.weights, bw=bw)
		S.t0.propose = sum(Surv.score.t0.i*new.weights)/sum(new.weights)	
		new.mat = cbind(short.new, as.matrix(z.cov)[tl>landmark,])
		betahat.cox.t <- coxph(Surv(tl[tl>landmark], dl[tl>landmark])~new.mat, weights = new.weights[tl>landmark])  
		score.t.i = c(as.matrix(new.mat)%*%betahat.cox.t$coef)
		Surv.score.t.i = Est.KM.FUN.weighted(xi=tl[tl>landmark], di=dl[tl>landmark], si=score.t.i, myt=tt, weight.perturb = new.weights[tl>landmark], bw=bw)
		S.t.obs = (sum(Surv.score.t.i*new.weights[tl>landmark])/sum(new.weights[tl>landmark]))*S.t0.propose
	}
	return(S.t.obs)
}

surv.land.obs = function(tl, dl, tt, landmark, short=NULL, z.cov=NULL, var = FALSE, conf.int = FALSE, ps.weights, weight.perturb=NULL, perturb.ps=FALSE, perturb.vector=FALSE,bw=NULL)	{
	S.est = surv.land.obs.base(tl=tl, dl=dl, tt=tt, landmark=landmark, short=short, z.cov=z.cov, ps.weights = ps.weights)
	if(var | conf.int)	{
		if(is.null(weight.perturb)){	
			weight.perturb = matrix(rexp(500*length(tl), rate=1), ncol = 500)
			S.p = apply(weight.perturb, 2, surv.land.obs.base, tl=tl, dl=dl, tt = tt, landmark=landmark, short=short, z.cov=z.cov,ps.weights = ps.weights, perturb.ps = FALSE)
		}
		if(!is.null(weight.perturb)){	
			S.p = apply(weight.perturb, 2, surv.land.obs.base, tl=tl, dl=dl, tt = tt, landmark=landmark, short=short, z.cov=z.cov,ps.weights = ps.weights, perturb.ps = perturb.ps)
		}
		if(conf.int)	{
			conf.l.normal.S = S.est - 1.96*sd(S.p)
			conf.u.normal.S = S.est + 1.96*sd(S.p)
			conf.l.quantile.S = quantile(S.p, 0.025)
			conf.u.quantile.S = quantile(S.p, 0.975)
		}
	}
	if(!(var) & !(conf.int))	{
		return(list("S.estimate" = S.est))
	}
	if(var & !(conf.int)& !perturb.vector)	{
		return(list("S.estimate" = S.est, "S.var" = var(S.p)))
	}
	if(conf.int & !perturb.vector)	{
		return(list("S.estimate" = S.est, "S.var" = var(S.p), "conf.int.normal.S" = as.vector(c(conf.l.normal.S, conf.u.normal.S)), "conf.int.quantile.S" = as.vector(as.vector(c(conf.l.quantile.S, conf.u.quantile.S)))))
	}	
	if(var & !(conf.int)& perturb.vector)	{
		return(list("S.estimate" = S.est, "S.var" = var(S.p), "perturb.vector" = S.p))
	}
	if(conf.int & perturb.vector)	{
		return(list("S.estimate" = S.est, "S.var" = var(S.p), "conf.int.normal.S" = as.vector(c(conf.l.normal.S, conf.u.normal.S)), "conf.int.quantile.S" = as.vector(as.vector(c(conf.l.quantile.S, conf.u.quantile.S))), "perturb.vector" = S.p))
	}	

}

delta.land.obs = function(tl, dl, treat, tt, landmark, short=NULL, z.cov=NULL, var = FALSE, conf.int = FALSE, ps.weights=NULL, weight.perturb=NULL, perturb.ps=FALSE,cov.for.ps = NULL,bw=NULL)	{
	if(sum(unique(treat) %in% c(0,1)) != 2){
		treat.f = as.factor(treat)
		unique.treat = unique(treat.f)[1]
		treat = 1*(treat.f==unique.treat)
		print(paste("Treatment variable was not 0/1; converted to 0/1 by 1=", unique.treat,", 0=", unique(treat.f)[2], sep =""))
	}
	if(!var & !conf.int)	{
		weightperturb.1 = NULL
		weightperturb.0 = NULL
	}
	if(var | conf.int)	{
		if(is.null(weight.perturb)){
			weight.perturb = matrix(rexp(500*length(tl), rate=1), ncol = 500)
			weightperturb.1 = weight.perturb[treat==1,]
			weightperturb.0 = weight.perturb[treat==0,]
		}
		if(!is.null(weight.perturb)){
		weightperturb.1 = weight.perturb[treat==1,]
		weightperturb.0 = weight.perturb[treat==0,]
		}
	}	
	if(is.null(z.cov) & is.null(short)) {stop("must supply either baseline covariate information or intermediate event information, otherwise you could just use regular Kaplan Meier.")}
	if(!is.null(short)){
		d = dim(short)[2]
		if(d %% 2 > 0) {stop("the number of columns for the intermediate event information matrix needs to be a multiple of 2.")}
	}	
	if(is.null(ps.weights))	{
		if(is.null(cov.for.ps))	{stop("must either supply propensity score weights or supply covariate information for the function to construct weights")}
		ps.weights = ps.wgt.fun(treat=treat, cov.for.ps = cov.for.ps)
		if(var | conf.int){
			ps.weights.perturb = apply(weight.perturb, 2, ps.wgt.fun, treat=treat, cov.for.ps = cov.for.ps)
			weightperturb.1 = weight.perturb[treat==1,]*ps.weights.perturb[treat==1,]
			weightperturb.0 = weight.perturb[treat==0,]*ps.weights.perturb[treat==0,]
			perturb.ps = TRUE
		}
		if(!var & !conf.int)	{
			perturb.ps=FALSE
		}	
	}	
	if(!is.null(z.cov) & !is.null(short)) {
		surv.1 = surv.land.obs(tl=tl[treat==1], dl = dl[treat==1], tt=tt, landmark = landmark, short = short[treat==1,], z.cov = as.matrix(z.cov)[treat==1,], var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.1, ps.weights=ps.weights[treat==1], perturb.ps=perturb.ps) 
		surv.0 = surv.land.obs(tl=tl[treat==0], dl = dl[treat==0], tt=tt, landmark = landmark, short = short[treat==0,], z.cov = as.matrix(z.cov)[treat==0,], var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.0, ps.weights=ps.weights[treat==0], perturb.ps=perturb.ps) 
	}
	if(is.null(z.cov)) {
		surv.1 = surv.land.obs(tl=tl[treat==1], dl = dl[treat==1], tt=tt, landmark = landmark, short = short[treat==1,], z.cov = NULL, var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.1, ps.weights=ps.weights[treat==1], perturb.ps=perturb.ps) 
		surv.0 = surv.land.obs(tl=tl[treat==0], dl = dl[treat==0], tt=tt, landmark = landmark, short = short[treat==0,], z.cov = NULL, var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.0, ps.weights=ps.weights[treat==0], perturb.ps=perturb.ps) 
	}
	if(is.null(short)) {
		surv.1 = surv.land.obs(tl=tl[treat==1], dl = dl[treat==1], tt=tt, landmark = landmark, short = NULL, z.cov = as.matrix(z.cov)[treat==1,], var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.1, ps.weights=ps.weights[treat==1], perturb.ps=perturb.ps) 
		surv.0 = surv.land.obs(tl=tl[treat==0], dl = dl[treat==0], tt=tt, landmark = landmark, short =  NULL, z.cov = as.matrix(z.cov)[treat==0,], var = var, conf.int = conf.int, perturb.vector = (var | conf.int), weight.perturb=weightperturb.0, ps.weights=ps.weights[treat==0], perturb.ps=perturb.ps) 
	}
	delta = surv.1$S.estimate - surv.0$S.estimate
	if(var | conf.int){
		delta.p =  surv.1$perturb.vector - surv.0$perturb.vector
		delta.var = var(delta.p)
	}
	if(conf.int)	{
		conf.l.normal.d = delta - 1.96*sd(delta.p)
		conf.u.normal.d = delta + 1.96*sd(delta.p)
		conf.l.quantile.d = quantile(delta.p, 0.025)
		conf.u.quantile.d = quantile(delta.p, 0.975)
	}
	if(var | conf.int){
		pval = pnorm(delta/sqrt(delta.var), lower.tail = (delta/sqrt(delta.var) < 0))*2
	}
	if(!var & !conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate, "delta.estimate" = delta))}
	if(var & !conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate,"delta.estimate" = delta, "S.var.1" = surv.1$S.var, "S.var.0" = surv.0$S.var, "delta.var" = delta.var, "p.value" = pval))}
	if(conf.int)	{return(list("S.estimate.1" = surv.1$S.estimate, "S.estimate.0" = surv.0$S.estimate, "delta.estimate" = delta, "S.var.1" = surv.1$S.var, "S.var.0" = surv.0$S.var, "delta.var" = delta.var,  "conf.int.normal.S.1" = surv.1$conf.int.normal.S,  "conf.int.normal.S.0" = surv.0$conf.int.normal.S, "conf.int.normal.delta" = c(conf.l.normal.d, conf.u.normal.d), "conf.int.quantile.S.1" = surv.1$conf.int.quantile.S,  "conf.int.quantile.S.0" = surv.0$conf.int.quantile.S,  "conf.int.quantile.delta" = as.vector(c(conf.l.quantile.d, conf.u.quantile.d)), "p.value" = pval))}
}


Est.KM.FUN.weighted <- function(xi,di,si,myt, weight.perturb = NULL, bw = NULL)
  {    
    if(is.null(bw))	{
    	bwini = bw.nrd(si)
    	n.s = length(si)
    	bw <- bwini/(n.s^0.11)
    }
    kerni.ss = Kern.FUN(si,si,bw)           
    tmpind = (xi<=myt)&(di==1); tj = xi[tmpind]; 
    if(is.null(weight.perturb)) {weight.perturb = rep(1, length(si))}
    kerni.1 = weight.perturb*kerni.ss
    pihamyt0.tj.ss = helper.si(tj, "<=", xi, Vi=kerni.1)   
    dLamhat.tj.ss = (kerni.1[tmpind,])/pihamyt0.tj.ss; 
	dLamhat.tj.ss[is.na(dLamhat.tj.ss)] = 0
    ret = apply(dLamhat.tj.ss,2,sum)
    Shat.ss  =exp(-ret)
    return(Shat.ss)
    }



 ps.wgt.fun <- function(treat, cov.for.ps, weight=NULL)
  {
    if(is.null(weight)) {weight = rep(1,length(treat))}
    prop.score = predict(glm(treat~cov.for.ps, family = "binomial", weights = weight), type = "response")
    Wi <- rep(0,length(treat))
    Wi[treat==1] <- 1/prop.score[treat==1]
    Wi[treat==0] <- 1/(1-prop.score[treat==0])
    return(Wi)
  }
  
Kern.FUN <- function(zz,zi,bw) 
  { 
    out = (VTM(zz,length(zi))- zi)/bw
	dnorm(out)/bw
           
  }
  
cumsum2 <- function(mydat)     
  {
    if(is.null(dim(mydat))) return(cumsum(mydat))
    else{
      out <- matrix(cumsum(mydat), nrow=nrow(mydat))
      out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
      return(out)
    }
  }

helper.si <- function(yy,FUN,Yi,Vi=NULL)   
  {  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
    if(is.null(Vi)){return(pos)}else{
      Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
      out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
      out[pos!=0,] <- Vi[pos,]
      if(is.null(dim(Vi))) out <- c(out)
      return(out) 
    }
  } 

      
  
  VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }
