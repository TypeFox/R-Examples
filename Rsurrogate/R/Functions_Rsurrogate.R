delta.estimate= function(yone,yzero, var = FALSE, conf.int = FALSE, weight = NULL, weight.perturb=NULL) {
	if(is.null(weight)){weight = rep(1,length(yone)+length(yzero))}
	delta = sum(weight[(1:length(yone))]*yone)/sum(weight[(1:length(yone))]) - sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*yzero)/(sum(weight[(length(yone)+1):(length(yone)+length(yzero))]))
	if(var | conf.int) { 
		if(is.null(weight.perturb)) {weight.perturb = matrix(rexp(500*(length(yone)+length(yzero)), rate=1), ncol = 500)}
		delta.p.vec = apply(weight.perturb, 2, function(x) sum(x[(1:length(yone))]*yone)/sum(x[(1:length(yone))]) - sum(x[(length(yone)+1):(length(yone)+length(yzero))]*yzero)/(sum(x[(length(yone)+1):(length(yone)+length(yzero))])) )
		}
	if(conf.int){
		conf.l.normal = delta - 1.96*sd(delta.p.vec)
		conf.u.normal = delta + 1.96*sd(delta.p.vec)
		conf.l.quantile = quantile(delta.p.vec, 0.025)
		conf.u.quantile = quantile(delta.p.vec, 0.975)
	}	
	if(!var & !conf.int) {return(list("delta" = delta))}
	if(var & !conf.int) {return(list("delta" = delta, "delta.var" = var(delta.p.vec)))}
	if(conf.int) {return(list("delta" = delta, "delta.var" = var(delta.p.vec), "conf.int.normal" = c(conf.l.normal, conf.u.normal), "conf.int.quantile" = as.vector(c(conf.l.quantile, conf.u.quantile))))}
}


delta.s.estimate = function(sone,szero,yone,yzero, weight.perturb = NULL, number,  type) {
	if(is.null(number)) {number = "single"}
	if(is.null(type)) {type = "robust"}
	if(substr(number,1,1) == "s") {number = "single"}
	if(substr(number,1,1) == "m") {number = "multiple"}
	if(substr(type,1,1) == "r") {type = "robust"}
	if(substr(type,1,1) == "m") {type = "model"}
	if(!(substr(type, 1,1) %in% c("r", "m"))) {print("Invalid type, default `robust' being used"); type = "robust"}
	if(!(substr(number, 1,1) %in% c("s", "m"))) {print("Invalid number, default `single' being used"); number = "single"}
	if(number == "multiple") {
		if(dim(as.matrix(sone))[2] <2) {print("Invalid number selection, there does not appear to be multiple surrogate markers; default `single' being used"); number = "single" }
	}
	if(number == "single" & dim(as.matrix(sone))[2] > 1) {print("Single setting being used but dimension of surrogate is greater than one; using first surrogate"); sone = as.matrix(sone)[,1]; szero = as.matrix(szero)[,1] }
	p.test = wilcox.test(x=yone, y=yzero, alternative = "two.sided")$p.value
	if(p.test > 0.05) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting")}
	if(number == "single" & type == "robust") {
	range.1 = range(sone)
	range.0 = range(szero)
	if((range.1[1] > range.0[1]) | (range.1[2] < range.0[2])) {
		print("Warning: observed supports to not appear equal, may need to consider a transformation or extrapolation")
	}
	}
	if(is.null(weight.perturb)) {weight = rep(1,length(yone)+length(yzero))}
	if(!is.null(weight.perturb)) {weight = weight.perturb}	
	if(type == "model" & number == "single") {
		reg.y = lm(c(yone,yzero)~c(sone,szero) + c(rep(1, length(sone)), rep(0,length(szero))) + c(rep(1, length(sone)), rep(0,length(szero)))*c(sone,szero), weights = weight)
		beta0 = reg.y$coeff[1]
		beta1 = reg.y$coeff[2]
		beta2 = reg.y$coeff[3]
		beta3 = reg.y$coeff[4]
		reg.s = lm(c(sone,szero) ~ c(rep(1, length(sone)), rep(0,length(szero))), weights = weight)
		alpha0 = reg.s$coeff[1]
		alpha1 = reg.s$coeff[2]
		m = beta0 + (beta1+beta3)*szero + beta2 
		delta.s = sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*(m-yzero))/sum(weight[(length(yone)+1):(length(yone)+length(yzero))])
	}
	if(type == "robust"  & number == "single") {
		h.select = bw.nrd(sone)*(length(sone)^(-0.10))
		mu.1.s0 = sapply(szero,pred.smooth,zz=sone, bw=h.select, y1=yone, weight = weight[(1:length(yone))])
		delta.s = sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*(mu.1.s0-yzero))/sum(weight[(length(yone)+1):(length(yone)+length(yzero))])	
	}
	if(type == "model"  & number == "multiple") {
		working = as.matrix(lm(yone~sone, weights = weight[(1:length(yone))] )$coeff)
		s.tilde.0 = cbind(rep(1,length(szero[,1])),szero)%*%working
		s.tilde.1 = cbind(rep(1,length(sone[,1])),sone)%*%working
		delta.s = sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*(s.tilde.0-yzero))/sum(weight[(length(yone)+1):(length(yone)+length(yzero))])
	}
	if(type == "robust" & number == "multiple") {
		working = as.matrix(lm(yone~sone, weights = weight[c(1:length(yone))])$coeff)
		s.tilde.0 = cbind(rep(1,length(szero[,1])),szero)%*%working
		s.tilde.1 = cbind(rep(1,length(sone[,1])),sone)%*%working
		mu.1.s0 = sapply(s.tilde.0,pred.smooth,zz=s.tilde.1, y1=yone, weight = weight[(1:length(yone))])
		delta.s = sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*(mu.1.s0-yzero))/sum(weight[(length(yone)+1):(length(yone)+length(yzero))])
	}
	return(delta.s)
}

R.s.estimate = function(sone,szero,yone,yzero, var = FALSE, conf.int = FALSE, weight.perturb = NULL,number, type) {
	if(is.null(number)) {number = "single"}
	if(is.null(type)) {type = "robust"}
	if(substr(number,1,1) == "s") {number = "single"}
	if(substr(number,1,1) == "m") {number = "multiple"}
	if(substr(type,1,1) == "r") {type = "robust"}
	if(substr(type,1,1) == "m") {type = "model"}
	if(substr(type,1,1) == "f") {type = "freedman"}
	if(!(substr(type, 1,1) %in% c("r", "m", "f"))) {print("Warning: Invalid type, default `robust' being used"); type = "robust"}
	if(!(substr(number, 1,1) %in% c("s", "m"))) {print("Warning: Invalid number, default `single' being used"); number = "single"}
	if(number == "multiple") {
		if(dim(as.matrix(sone))[2] <2) {print("Invalid number selection, there does not appear to be multiple surrogate markers; default `single' being used"); number = "single" }
	}
	if(number == "single" & dim(as.matrix(sone))[2] > 1) {print("Single setting being used but dimension of surrogate is greater than one; using first surrogate"); sone = as.matrix(sone)[,1]; szero = as.matrix(szero)[,1] }
	p.test = wilcox.test(x=yone, y=yzero, alternative = "two.sided")$p.value
	if(p.test > 0.05) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting")}
	if(mean(yone)-mean(yzero) < 0) {print("Warning: it looks like you need to switch the treatment groups")}
	if(number == "single" & type == "robust") {
	range.1 = range(sone)
	range.0 = range(szero)
	if((range.1[1] > range.0[1]) | (range.1[2] < range.0[2])) {
		print("Warning: observed supports to not appear equal, may need to consider a transformation or extrapolation")
	}
	}
	if(type != "freedman"){
		delta = delta.estimate(yone = yone, yzero=yzero)$delta	
		delta.s = delta.s.estimate(sone=sone, szero=szero, yone=yone, yzero=yzero, number = number, type = type)
		R.s = 1-delta.s/delta	
	}
	if(type == "freedman"  & number == "single") {
		reg.y = lm(c(yone,yzero)~c(sone, szero) + c(rep(1, length(sone)),rep(0,length(szero))))
		beta.treat.s = reg.y$coeff[3]
		reg.y = lm(c(yone,yzero)~c(rep(1, length(sone)),rep(0,length(szero))))
		beta.treat.nos = reg.y$coeff[2]
		R.s = 1-beta.treat.s/beta.treat.nos
	}
	if(type == "freedman"  & number == "multiple") {
		treat.ind = c(rep(1,length(yone)), rep(0, length(yzero)))
		beta.treat.s = summary(lm(c(yone, yzero)~cbind(treat.ind, rbind(sone, szero))))$coeff[2]
		beta.treat.nos = summary(lm(c(yone, yzero)~treat.ind))$coeff[2]
		R.s = 1-beta.treat.s/beta.treat.nos
	}
	if(var | conf.int){
		if(is.null(weight.perturb)) {
	weight.perturb = matrix(rexp(500*(length(yone)+length(yzero)), rate=1), ncol = 500)}
		if(type != "freedman"){
			delta.s.p.vec = apply(weight.perturb, 2, delta.s.estimate, sone=sone, szero=szero, yone = yone, yzero = yzero,number = number, type = type)
			delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.estimate, yone = yone, yzero=yzero, var= FALSE, conf.int = FALSE)))
			R.p = 1-delta.s.p.vec/delta.p.vec
		}
		if(type == "freedman"  & number == "single") {
			beta.treat.s.p<- apply( weight.perturb, 2 , function(x) (lm(c(yone,yzero)~c(sone, szero) + c(rep(1, length(sone)),rep(0,length(szero))), weights = x) )$coeff[3])
			beta.treat.nos.p = apply( weight.perturb, 2 , function(x) (lm(c(yone,yzero)~ c(rep(1, length(sone)),rep(0,length(szero))), weights = x) )$coeff[2])
			R.p = 1-beta.treat.s.p/beta.treat.nos.p
		}
		if(type == "freedman"  & number == "multiple") {
		beta.treat.s.p<- apply(weight.perturb, 2 , function(x) (lm(c(yone,yzero)~cbind(treat.ind, rbind(sone, szero)), weights = x))$coeff[2])
		beta.treat.nos.p = apply( weight.perturb, 2 , function(x) (lm(c(yone, yzero)~treat.ind,weights = x))$coeff[2])
			R.p = 1-beta.treat.s.p/beta.treat.nos.p
		}
	if(conf.int & type != "freedman")	{
		conf.l.normal.delta = delta - 1.96*sd(delta.p.vec)
		conf.u.normal.delta = delta + 1.96*sd(delta.p.vec)
		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.normal.delta.s = delta.s - 1.96*sd(delta.s.p.vec)
		conf.u.normal.delta.s = delta.s + 1.96*sd(delta.s.p.vec)
		conf.l.quantile.delta.s = quantile(delta.s.p.vec, 0.025)
		conf.u.quantile.delta.s = quantile(delta.s.p.vec, 0.975)
	}
	if(conf.int) {	
		conf.l.normal.R.s = R.s - 1.96*sd(R.p)
		conf.u.normal.R.s = R.s + 1.96*sd(R.p)
		conf.l.quantile.R.s = quantile(R.p, 0.025)
		conf.u.quantile.R.s = quantile(R.p, 0.975)
		if(type!= "freedman") {fieller.ci.calc = as.vector(fieller.ci(delta.s.p.vec, delta.p.vec, delta.s , delta))}
		if(type== "freedman") {fieller.ci.calc=as.vector(fieller.ci(beta.treat.s.p, beta.treat.nos.p, beta.treat.s, beta.treat.nos))}
	}	
	}
	if(type != "freedman") {
		if(!var & !conf.int) {return(list("delta" = delta, "delta.s" =delta.s, "R.s" = R.s))}
		if(var & !conf.int) {return(list("delta" = delta, "delta.s" =delta.s, "R.s" = R.s, "delta.var" = var(delta.p.vec), "delta.s.var" = var(delta.s.p.vec), "R.s.var" = var(R.p)))}
		if(conf.int) {return(list("delta" = delta, "delta.s" =delta.s, "R.s" = R.s, "delta.var" = var(delta.p.vec), "delta.s.var" = var(delta.s.p.vec), "R.s.var" = var(R.p), "conf.int.normal.delta" = c(conf.l.normal.delta, conf.u.normal.delta), "conf.int.quantile.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.normal.delta.s" = c(conf.l.normal.delta.s, conf.u.normal.delta.s), "conf.int.quantile.delta.s" = as.vector(c(conf.l.quantile.delta.s, conf.u.quantile.delta.s)),
"conf.int.normal.R.s" = c(conf.l.normal.R.s, conf.u.normal.R.s), "conf.int.quantile.R.s" = as.vector(c(conf.l.quantile.R.s, conf.u.quantile.R.s)), "conf.int.fieller.R.s" = fieller.ci.calc))}	
	}
	if(type == "freedman") {
		if(!var & !conf.int) {return(list("R.s" = as.vector(R.s)))}
		if(var & !conf.int) {return(list("R.s" = as.vector(R.s), "R.s.var" = var(R.p)))}
		if(conf.int) {return(list("R.s" = as.vector(R.s),  "R.s.var" = as.vector(var(R.p)), 
"conf.int.normal.R.s" = as.vector(c(conf.l.normal.R.s, conf.u.normal.R.s)), "conf.int.quantile.R.s" = as.vector(as.vector(c(conf.l.quantile.R.s, conf.u.quantile.R.s))), "conf.int.fieller.R.s" = fieller.ci.calc))}
	}
}

fieller.ci = function(perturb.delta.s, perturb.delta, delta.s, delta) {
	num = (perturb.delta.s-(delta.s/delta)*perturb.delta)^2
	sigma11 = var(perturb.delta.s)
	sigma22 = var(perturb.delta)
	sigma12 = cov(perturb.delta.s, perturb.delta)
	den = sigma11 - 2*(delta.s/delta)*sigma12 + (delta.s/delta)^2*sigma22
	c.alpha = quantile(num/den,probs=0.95)
	#quadratic
	a = delta^2-sigma22*c.alpha
	b = -2*delta.s*delta + c.alpha*sigma12*2
	c = delta.s^2 - c.alpha*sigma11
	root.1 = (-b+sqrt(b^2 - 4*a*c))/(2*a)
	root.2 = (-b-sqrt(b^2 - 4*a*c))/(2*a)
	return(c(1-root.1,1-root.2))
	
}

VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }
    
 
Kern.FUN <- function(zz,zi,bw) 
  { 
    out = (VTM(zz,length(zi))- zi)/bw
	dnorm(out)/bw
           
  }
  

pred.smooth <-function(zz,zi.one, bw=NULL,y1, weight = NULL) {
  	if(is.null(bw)) { bw = bw.nrd(zz)/((length(zz))^(0.10))}
  	if(is.null(weight)) {weight = rep(1,length(y1))}
  	return(sum(weight*Kern.FUN(zz,zi.one,bw=bw)*y1)/sum(weight*Kern.FUN(zz,zi.one,bw=bw)))
  }

