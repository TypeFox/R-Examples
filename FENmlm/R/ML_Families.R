
#Les fonctions des "family":
#	- grad.NL : prend en argument la jacobienne non lineaire et renvoie le gradient de la LL
#	- grad.L : prend en argument les donnees et renvoie le gradient
#	- hess : prend e argument les jacobiennes et renvoie la hessienne (il restera qd meme la partie NL de la hessienne a calculer)
#	- ll: la vraisemblance, prend les coefs et renvoie la vraisemblance
#	- dummy: prend en compte la variable factorielle et renvoie l'optimisation vav de ces variables (seulement pour le Poisson pour l'instant). Je pourrais aussi changer les dummies par des sparse matrices mais je ne sais pas si cela ira plus vite.

# TODO
# write functions to get the STDerror of the  intercept and the dummies
#Write the scores.


#	- poisson
#	- logit
#	- probit
#	- tobit
#	- negative binomial
#  - gaussian

#***************************#
#### ===== POISSON ===== ####
#***************************#

ml_poisson = function(){
	#Cette fonction renvoie une famille de fonctions
	ll = function(y,mu,env,...){
		sum(y*mu - exp(mu) - lfactorial(y))
	}

	#Derivee
	ll_dl = function(y,mu,...){
		c(y - exp(mu))
	}

	#Derivee seconde
	ll_d2 = function(y,mu,...){
		-c(exp(mu))
	}

	ll_TEST_score = function(y,mu,env,...){
		#used to compute the scores numerically
		#Not to be used by end-users
		c(y*mu - exp(mu) - lfactorial(y))
	}

	expected.predictor = function(mu,env,...){
		exp(mu)
	}

	closedFormDummies = function(dum,S,y,mu,env,...){
		#We send only the dummies (not the vector of dummies)
		sum_y_i = as.vector(S%*%y)
		mu_dum = as.vector(S%*%c(exp(mu)))
		log(sum_y_i) - log(mu_dum)
	}

	ll0 = function(cste,y){
		n = length(y)
		#la fonction minimise, on renvoie "-"
		- (cste * sum(y) - n * exp(cste) - sum(lfactorial(y)))
	}

	grad0 = function(cste,y){
		#la fonction minimise, on renvoie "-"
		- sum(y-exp(cste))
	}

	model0 = function(y,mu){
		# TO be REDEFINED !!
		x = log(sum(y)) - log(sum(exp(mu)))
		ll = sum(y*x - exp(x) - lfactorial(y))
		list(loglik=ll,constant=x)
	}

	return(list(ll=ll,expected.predictor=expected.predictor,ll0=ll0,grad0=grad0,ll_TEST_score=ll_TEST_score,ll_dl=ll_dl,ll_d2=ll_d2,model0=model0,closedFormDummies=closedFormDummies))
}

#************************#
#### ===== LOGIT ==== ####
#************************#


ml_logit = function(){

	ll = function(y,mu,env,...){
		sum(y*mu - log(1+exp(mu)))
	}

	#Derivee
	ll_dl = function(y,mu,...){
		c(y - exp(mu)/(1+exp(mu)))
	}

	#Derivee seconde
	ll_d2 = function(y,mu,...){
		- c(exp(mu)/(1+exp(mu))^2)
	}

	ll_TEST_score = function(y,mu,env,...){
		#used to compute the scores numerically
		#Not to be used by end-users
		c(y*mu - log(1+exp(mu)))
	}

	expected.predictor = function(mu,env,...){
		exp(mu)/(1+exp(mu))
	}

	initDummy = function(S,y,mu,env,...){
		#guess for the dummy:
		ni1 = as.vector(S%*%y)
		ni = rowSums(S)
		log(ni1) - log(ni-ni1) - as.vector(S%*%mu)/ni
	}

	guessDummy = function(sum_y,n_group,mu,...){
		#guess for the dummy:
		log(sum_y) - log(n_group-sum_y) - mu
	}

	dum_fx = function(x,sum_y,S,mu,dum,...){
		sum_y - as.vector(S%*%c(exp(mu+x[dum]) / (1+exp(mu+x[dum]))))
	}

	dum_dfx = function(x,S,mu,dum){
		-as.vector(S%*%c(exp(mu+x[dum]) / (1+exp(mu+x[dum]))^2))
	}

	ratio_fx_dfx = function(x,dum,S,y,mu,env,...){
		ni1 = as.vector(S%*%y)

		m = as.vector(S%*%c(exp(mu+x[dum]) / (1+exp(mu+x[dum]))))
		Value = c(ni1 - m)
		Derivee = -as.vector(S%*%c(exp(mu+x[dum]) / (1+exp(mu+x[dum]))^2))
		Value / Derivee
	}

	ll0 = function(cste,y){
		n = length(y)
		#la fonction minimise, on renvoie "-"
		- (cste*sum(y) - n*log(1+exp(cste)))
	}

	grad0 = function(cste,y){
		#la fonction minimise, on renvoie "-"
		- sum(y-exp(cste)/(1+exp(cste)))
	}

	model0 = function(y){
		n = length(y)
		x = log(sum(y)) - log(length(y) - sum(y))
		ll = x*sum(y) - n*log(1+exp(x))
		list(loglik=ll,constant=x)
	}

	return(list(ll=ll,expected.predictor=expected.predictor,ll0=ll0,grad0=grad0,ll_TEST_score=ll_TEST_score,model0=model0,ll_dl=ll_dl,ll_d2=ll_d2,ratio_fx_dfx=ratio_fx_dfx,initDummy=initDummy,guessDummy=guessDummy,dum_fx=dum_fx,dum_dfx=dum_dfx))
}


#*************************#
#### ===== NEGBIN ==== ####
#*************************#

ml_negbin = function(){
	#Cette fonction renvoie une famille de fonctions
	ll = function(y,mu,env,coef,...){
		theta = coef[".theta"]
		#theta = get(".theta",env)
		#theta = attr(mu,".theta")
		sum(lgamma(theta+y) - lgamma(theta) - lgamma(y+1) + theta*log(theta) + y*mu - (theta+y)*log(theta+exp(mu)))
	}

	#Derivee
	ll_dl = function(y,mu,coef,env,...){
		theta = coef[".theta"]
		#theta = get(".theta",env)
		#theta = attr(mu,".theta")
		c(y - (theta+y)*exp(mu)/(theta+exp(mu)))
	}

	#Derivee croisee
	ll_dx_dother = function(y,mu,coef,env,...){
		#Means the second derivative of the LL wrt to the linear part and theta
		theta = coef[".theta"]
		#theta = get(".theta",env)
		#theta = attr(mu,".theta")
		c(-exp(mu)*(exp(mu)-y)/(theta+exp(mu))^2)
	}

	#Derivee seconde:
	ll_d2 = function(y,mu,coef,env,...){
		theta = coef[".theta"]
		#theta = get(".theta",env)
		#theta = attr(mu,".theta")
		- theta * (theta+y) * c(exp(mu)/(theta+exp(mu))^2)
	}

	ll_TEST_score = function(y,mu,env,coef,...){
		#used to compute the scores numerically
		#Not to be used by end-users
		theta = coef[length(coef)]
		#theta = attr(mu,".theta")
		c(lgamma(theta+y) - lgamma(theta) - lgamma(y+1) + theta*log(theta) + y*mu - (theta+y)*log(theta+exp(mu)))
	}

	grad.theta = function(theta,y,mu,...){
		sum( psigamma(theta+y) - psigamma(theta) + log(theta) + 1 - log(theta+exp(mu)) - (theta+y)/(theta+exp(mu)) )
	}

	scores.theta = function(theta,y,mu){
		psigamma(theta+y) - psigamma(theta) + log(theta) + 1 - log(theta+exp(mu)) - (theta+y)/(theta+exp(mu))
	}

	hess.theta = function(theta,y,mu){
		sum( psigamma(theta+y,1) - psigamma(theta,1) + 1/theta - 1/(theta+exp(mu)) + (y-exp(mu))/(theta+exp(mu))^2 )
	}

	hess.thetaL = function(theta,jacob.mat,y,dxi_dbeta,dxi_dother,ll_d2,ll_dx_dother){
		res = crossprod((jacob.mat+dxi_dbeta),dxi_dother*ll_d2+ll_dx_dother)
		return(res)
		#DEPREC
		#if there's no dummies, then dxi_dbeta = 0
# 		ll_dxdt = exp(mu)*(y-exp(mu))/(theta+exp(mu))^2
#
# 		d_mu = jacob.mat + dxi_dbeta
# 		H = crossprod(d_mu,dxi_dbeta * ll_d2 + ll_dxdt)
#		#The old version was:
#		#H = crossprod(d_mu,ll_dxdt)
# 		return(as.matrix(H))
	}

	hess_theta_part = function(theta,y,mu,dxi_dother,ll_dx_dother,ll_d2){
		#La derivee vav de theta en prenant en compte les dummies
		d2ll_d2theta = hess.theta(theta,y,mu)
		res = sum(dxi_dother^2*ll_d2 + 2*dxi_dother*ll_dx_dother) + d2ll_d2theta
		return(res)
		#DEPREC -- when theta-conditionned
# 		#Here we add the part of the Hessian
# 		#That corresponds to the theta profiled
# 		e_mu = exp(mu)
# 		n = length(y)
# # print(n)
# 		#First we get the sum(d2ll_d2theta)
# 		d2l_d2t = sum( psigamma(theta+y,1) - 1/(theta+e_mu) + (y-e_mu)/(theta+e_mu)^2 ) - n*psigamma(theta,1) + n*1/theta
# # 		print(d2l_d2t)
#
# 		#Then we get the dtheta_dbeta
# 		d2ll_dxdt = e_mu*(y-e_mu)/(theta+e_mu)^2
#
# 		dt_db = - colSums(jacob.mat * d2ll_dxdt) / d2l_d2t
#
# 		#Finally:
# # 		print(tcrossprod(dt_db) * d2l_d2t)
# 		- tcrossprod(dt_db) * d2l_d2t

	}

	scores_theta_part = function(theta,jacob.mat,y,mu){
		#We add the part of the score due to the profiled theta
		e_mu = exp(mu)
		n = length(y)

		#First we get the sum(d2ll_d2theta)
		d2l_d2t = sum( psigamma(theta+y,1) - 1/(theta+e_mu) + (y-e_mu)/(theta+e_mu)^2 ) - n*psigamma(theta,1) + n*1/theta

		#Then we get the dtheta_dbeta
		d2ll_dxdt = e_mu*(y-e_mu)/(theta+e_mu)^2
		dt_db = - colSums(jacob.mat * d2ll_dxdt) / d2l_d2t

		#We also get dll_dtheta
		dl_dt = psigamma(theta+y) - psigamma(theta) + log(theta) + 1 - log(theta+e_mu) - (theta+y)/(theta+e_mu)

		#Finally:
		tcrossprod(dl_dt,dt_db)
	}

	expected.predictor = function(mu,env,...){
		exp(mu)
	}

	ratio_fx_dfx = function(x,dum,S,y,mu,env,coef,...){
		sum_yi = as.vector(S%*%y)
		theta = coef[".theta"]
 		#theta = get(".theta",env)
		m = as.vector(S%*%c((theta+y) * exp(mu+x[dum]) / (theta+exp(mu+x[dum]))))
		Value = c(sum_yi - m)
		Derivee = - as.vector(S%*%c(theta * (theta+y) * exp(mu+x[dum]) / (theta+exp(mu+x[dum]))^2))
		Value / Derivee
	}

	initDummy = function(S,y,mu,env,...){
		sum_yi = as.vector(S%*%y)
		ni = rowSums(S)
		log(sum_yi) - log(ni) - as.vector(S%*%mu)/ni
	}

	guessDummy = function(sum_y,n_group,mu,...){
		#guess for the dummy:
		log(sum_y) - log(n_group) - mu
	}

	dum_fx = function(x,sum_y,S,mu,dum,y,coef,env,...){
		theta = coef[".theta"]
		#theta = get(".theta",env)
		sum_y - as.vector(S%*%c((theta+y) * exp(mu+x[dum]) / (theta+exp(mu+x[dum]))))
	}

	dum_dfx = function(x,S,mu,dum,y,coef,env,...){
		theta = coef[".theta"]
		#theta = get(".theta",env)
		- as.vector(S%*%c(theta * (theta+y) * exp(mu+x[dum]) / (theta+exp(mu+x[dum]))^2))
	}

	ll0 = function(coef,y){
		n = length(y)
		#La fonction minimise, on renvoie "-"
		theta = coef[2]
		cste = coef[1]
		ll = sum(lgamma(theta+y) - lgamma(theta) - lgamma(y+1) + theta*log(theta) + y*cste - (theta+y)*log(theta+exp(cste)))
		-ll
	}

	grad0 = function(coef,y){
		#La fonction minimise, on renvoie "-"
		theta = coef[2]
		cste = coef[1]

		grad.cste = sum(y - (theta+y)*exp(cste)/(theta+exp(cste)))
		grad.theta = sum(psigamma(theta+y) - psigamma(theta) + log(theta) + 1 - log(theta+exp(cste)) - (theta+y)/(theta+exp(cste)))

		- c(grad.cste,grad.theta)
	}

	init_theta = function(y,mu){
		#We know that when theta is close to 0, the foc tends to +infinity
		#So starting from a positive value is good because we are then 100%
		#sure that we will find a solution
		e_mu = exp(mu)
		n = length(y)

		theta = 1
		while( (Value <- sum(psigamma(theta+y) - log(1+e_mu/theta) - (y-e_mu)/(theta+e_mu)) - n*psigamma(theta)) < 0) theta = theta/10
		theta
	}

	ratio_fx_dfx_theta = function(theta,y,mu){
		e_mu = exp(mu)
		n = length(y)

		#fonction pour optimiser theta
		#Il faut juste empecher la fonction d'aller trop a gauche
		Value = sum(psigamma(theta+y) - log(1+e_mu/theta) - (y-e_mu)/(theta+e_mu)) - n*psigamma(theta)
		Derivee = sum(psigamma(theta+y,1) - 1/(theta+e_mu) + (y-e_mu)/(theta+e_mu)^2) - n*psigamma(theta,1) + n*1/theta
		if( (theta - Value / Derivee) <= 0){
			return(theta/2)
		} else{
			return(Value / Derivee)
		}
	}

	get_theta = function(theta,y,mu){
		#Attention a la precision de calcul! c'est important
		#quand on calcule les dummies
		show=FALSE
		if(show) cat("\n\nENTRANCE theta = ",theta,"\n")
		#This function gets the optimal value of theta wrt mu
		e_mu = exp(mu)
		n = length(y)

		# Condition de Premier Ordre:
		cpo = function(theta) sum(psigamma(theta+y) - log(1+e_mu/theta) - (y-e_mu)/(theta+e_mu)) - n*psigamma(theta)
		# Condition de Second Ordre
		cso = function(theta) sum(psigamma(theta+y,1) - 1/(theta+e_mu) + (y-e_mu)/(theta+e_mu)^2) - n*psigamma(theta,1) + n*1/theta

		# Step 1: get the boundaries + set the starting value
		Value = sum(psigamma(theta+y) - log(1+e_mu/theta) - (y-e_mu)/(theta+e_mu)) - n*psigamma(theta)
		if(Value<0){
			borne_sup = theta
			# Brute force to get a positive theta
			theta = theta/2
			while( cpo(theta) < 0 ) theta = theta/2
			borne_inf = theta
		} else{
			borne_inf = theta
			# Brute force to get a negative theta
			theta = 2*theta
			while( cpo(theta) > 0 ) theta = 2*theta
			borne_sup = theta
		}

		x1 = (borne_sup + borne_inf)/2

		# Setp 2: optimization

		precision = 1e-6 ; itermax = 100 ; value.precision = 1e-9
		iter = 0 ; ok = TRUE

		while(ok){
			# 1st step: initialisation des bornes
			Value = cpo(x1)
			if(Value>0) borne_inf = x1
			else borne_sup = x1

			if(show){
				cat("\nStart iter",iter,":")
				cat(sprintf("\nx1: %.10f\nValues: %e\n",x1,Value))
				cat(sprintf("borne_inf: %.10f\nborne_sup: %.10f\n",borne_inf,borne_sup))
			}

			# 2nd step: NR iteration
			x0 = x1
			Derivee = cso(x1)
			x1 = x0 - Value / Derivee

			if(show) cat(sprintf("X new (NR): %.10f\n",x1))

			#3rd step: dichotomy (if necessary)
			# Update of the value if it goes out of the boundaries
			if(x1>=borne_sup | x1<=borne_inf) x1 = (borne_inf + borne_sup)/2

			if(show) cat(sprintf("X new (dicho): %.10f\n",x1))

			if((iter<-iter+1)==itermax) stop("The algorithm getting the dummies diverged.\nThe error is unknown.")
			if(anyNA(x1)) stop("The algorithm getting the dummies diverged.\nThe error is unknown.")
			if(max(abs(x0-x1))<precision | abs(Value)<value.precision) ok = FALSE
		}
		x1
	}

	theta0 = function(y){
		# Quand il y a seulement la constante!
		y_bar = mean(y)
		n = length(y)

		ratio_fx_dfx_theta = function(x){
			Value = sum(psigamma(x+y) - (x+y)/(x+y_bar)) + n*( -psigamma(x) + log(x) + 1 - log(x+y_bar))
			Derivee = sum(psigamma(x+y,1)) + n*( -psigamma(x,1) + 1/x - 1/(x+y_bar))
			Value / Derivee
		}
		#il faut trouver la bonne initialisation de theta!
		init=555
		theta0 = qNR(init,ratio_fx_dfx_theta)

	}

	model0 = function(y){
		x = log(mean(y))
		theta = theta0(y)
		#get theta!
		list(loglik=ll,constant=x)
	}

	return(list(ll=ll,expected.predictor=expected.predictor,ll0=ll0,grad0=grad0,ll_TEST_score=ll_TEST_score,hess.thetaL=hess.thetaL,hess.theta=hess.theta,grad.theta=grad.theta,scores.theta=scores.theta,ll_dl=ll_dl,ll_d2=ll_d2,ratio_fx_dfx=ratio_fx_dfx,initDummy=initDummy,guessDummy=guessDummy,dum_fx=dum_fx,dum_dfx=dum_dfx,ratio_fx_dfx_theta=ratio_fx_dfx_theta,init_theta=init_theta,hess_theta_part=hess_theta_part,scores_theta_part=scores_theta_part,get_theta=get_theta,ll_dx_dother=ll_dx_dother))
}


#*************************#
#### ===== PROBIT ==== ####
#*************************#


ml_probit = function(){

	IHR = function(x){
		#Actually it's the inverse hazard rate
		#I have to create such a function to compute
		#properly the hazard rate that is not numerically
		#computable for large values of x
		#yet it tends to x
		v = dnorm(x)/pnorm(-x)
		qui = which(is.na(v) | !is.finite(v))
		v[qui] = x[qui]
		v
	}

	ll = function(y,mu,env,...){
		mu = mu * (2*y-1)
		sum(pnorm(mu,log.p = TRUE))
	}

	#Derivee
	ll_dl = function(y,mu,...){
		y_sign = 2*y-1
		#c(y_sign * dnorm(mu)/pnorm(mu*y_sign))
		c(y_sign * IHR(-mu*y_sign))
	}

	#Derivee seconde
	ll_d2 = function(y,mu,...){
		y_sign = 2*y-1
		#dnorm(mu)/pnorm(mu*y_sign)*(-(y_sign*mu) - y_sign*dnorm(mu)/pnorm(y_sign*mu))
		IHR(-mu*y_sign)*(-(y_sign*mu) - y_sign*IHR(-y_sign*mu))
	}

	ll_TEST_score = function(y,mu,env,...){
		#used to compute the scores numerically
		#Not to be used by end-users
		mu = mu * (2*y-1)
		c(-pnorm(mu,log.p = TRUE))
	}

	expected.predictor = function(mu,env,...){
		#mu = mu * (2*y-1)#=>no, not the LL, the exp.pred.
		pnorm(mu)
	}

	initDummy = function(S,y,mu,env,...){
		#guess for the dummy:
		#we assume mu is at the average level
		ni1 = as.vector(S%*%y)
		ni = rowSums(S)
		qnorm(ni1/ni)-as.vector(S%*%mu)/ni
	}

	guessDummy = function(sum_y,n_group,mu,...){
		#Dummy when mu is constant within the group
		qnorm(sum_y/n_group)-mu
	}

	dum_fx = function(x,sum_y,S,mu,dum,y,coef,env,...){
		mu = mu + x[dum]
		as.vector(S%*%c(IHR(-mu) * (y/pnorm(-mu) - IHR(mu))))
	}

	dum_dfx = function(x,sum_y,S,mu,dum,y,coef,env,...){
		stop("not yet implemented")
		-as.vector(S%*%c(exp(mu+x[dum]) / (1+exp(mu+x[dum]))^2))
	}

	ratio_fx_dfx = function(x,dum,S,y,mu,env,...){
		stop("not yet implemented")
		mu = mu + x[dum]
		Value = as.vector(S%*%c(IHR(-mu) * (y/pnorm(-mu) - IHR(mu))))

		Derivee = "to be written"
		Value / Derivee
	}

	ll0 = function(cste,y){
		y_sign = 2*y-1
		#la fonction minimise, on renvoie "-"
		- sum(pnorm(y_sign*cste,log.p = TRUE))
	}

	grad0 = function(cste,y){
		#la fonction minimise, on renvoie "-"
		y_sign = 2*y-1
		-c(y_sign * dnorm(cste)/pnorm(cste*y_sign))
	}

	model0 = function(y){
		stop("not yet implemented")
		n = length(y)
		x = log(sum(y)) - log(length(y) - sum(y))
		ll = x*sum(y) - n*log(1+exp(x))
		list(loglik=ll,constant=x)
	}

	return(list(ll=ll,expected.predictor=expected.predictor,ll0=ll0,grad0=grad0,ll_TEST_score=ll_TEST_score,model0=model0,ll_dl=ll_dl,ll_d2=ll_d2,ratio_fx_dfx=ratio_fx_dfx,initDummy=initDummy,guessDummy=guessDummy,dum_fx=dum_fx,dum_dfx=dum_dfx))
}

#***************************#
#### ===== GAUSSIAN ==== ####
#***************************#


ml_gaussian = function(){

	ll = function(y,mu,env,...){
		sigma = sqrt(mean((y-mu)^2))
		n = length(y)
		-1/2/sigma^2*sum((y-mu)^2) - n*log(sigma) - n*log(2*pi)/2
	}

	#Derivee
	ll_dl = function(y,mu,env,...){
		sigma = sqrt(mean((y-mu)^2))
		(y-mu)/sigma^2
	}

	#Derivee seconde
	ll_d2 = function(y,mu,env,...){
		sigma = sqrt(mean((y-mu)^2))
		rep(-1/sigma^2,length(y))
	}

	ll_TEST_score = function(y,mu,env,...){
		#used to compute the scores numerically
		#Not to be used by end-users
		sigma = sqrt(mean((y-mu)^2))
		c(-1/2/sigma^2*(y-mu)^2 - log(sigma) - log(2*pi)/2)
	}

	expected.predictor = function(mu,env,...){
		mu
	}

	closedFormDummies = function(dum,S,y,mu,env,...){
		#We send only the dummies (not the vector of dummies)
		as.vector(S%*%(y-mu)) / as.vector(rowSums(S))
	}

	ll0 = function(cste,y){
		#la fonction minimise, on renvoie "-"
		n = length(y)
		sigma = sqrt(sum((y-cste)^2)/n)
		-(-1/2/sigma^2*sum((y-cste)^2) - n*log(sigma) - n*log(2*pi)/2)
	}

	grad0 = function(cste,y){
		#la fonction minimise, on renvoie "-"
		n = length(y)
		sigma = sqrt(sum((y-cste)^2)/n)
		-sum((y-cste)/sigma^2)
	}

	model0 = function(y){
		n = length(y)
		x = mean(y)
		sigma = sqrt(sum((y-x)^2)/n)
		ll = -1/2/sigma^2*sum((y-x)^2) - n*log(sigma) - n*log(2*pi)/2
		list(loglik=ll,constant=x)
	}

	return(list(ll=ll,expected.predictor=expected.predictor,ll0=ll0,grad0=grad0,ll_TEST_score=ll_TEST_score,model0=model0,ll_dl=ll_dl,ll_d2=ll_d2,closedFormDummies=closedFormDummies))
}

#************************#
#### ===== TOBIT ==== ####
#************************#

# Not done yet! => I have to be careful and revise the
# derivative of the .sigma wrt the dummies

ml_tobit = function(){
	# we separate each case with U: uncensored and C: censored

	IHR = function(x){
		#Actually it's the inverse hazard rate
		#I have to create such a function to compute
		#properly the hazard rate that is not numerically
		#computable for large values of x
		#yet it tends to x
		v = dnorm(x)/pnorm(-x)
		qui = which(is.na(v) | !is.finite(v))
		v[qui] = x[qui]
		v
	}

	ll = function(y,mu,env,coef,...){
		#sigma = get(".sigma",env)
		sigma = coef[".sigma"]
		y_U = y[y>0]
		mu_U = mu[y>0]
		y_C = y[y==0]
		mu_C = mu[y==0]
		n_U = length(y_U)

		sum(pnorm(-mu_C/sigma,log.p = TRUE)) - 1/2/sigma^2*sum((y_U-mu_U)^2) - n_U*log(sigma) - n_U/2*log(2*pi)
	}

	#Derivee
	ll_dl = function(y,mu,env,coef,...){
		#sigma = get(".sigma",env)
		sigma = coef[".sigma"]
		y_U = y[y>0]
		mu_U = mu[y>0]
		y_C = y[y==0]
		mu_C = mu[y==0]

		ll_dl = rep(NA,length(y))
		ll_dl[y==0] = -1/sigma*IHR(mu_C/sigma)
		ll_dl[y>0] = 1/sigma^2*(y_U-mu_U)

		ll_dl
	}

	#Derivee seconde
	ll_d2 = function(y,mu,env,coef,...){
		#sigma = get(".sigma",env)
		sigma = coef[".sigma"]
		y_U = y[y>0]
		mu_U = mu[y>0]
		y_C = y[y==0]
		mu_C = mu[y==0]

		mu_C = mu_C/sigma
		coef_C = IHR(mu_C)

		ll_d2 = rep(NA,length(y))
		ll_d2[y==0] = 1/sigma^2*coef_C*( mu_C - coef_C )
		ll_d2[y>0] = - 1/sigma^2

		ll_d2
	}

	ll_TEST_score = function(y,mu,env,coef,...){
		#used to compute the scores numerically
		#Not to be used by end-users
		#sigma = get(".sigma",env)
		sigma = coef[".sigma"]
		y_U = y[y>0]
		mu_U = mu[y>0]
		y_C = y[y==0]
		mu_C = mu[y==0]

		n_U = length(y_U)

		ll = rep(NA,length(y))
		ll[y==0] = pnorm(-mu_C/sigma,log.p = TRUE)
		ll[y>0] = 1/2/sigma^2*(y_U-mu_U)^2 - n_U*log(sigma) - n_U/2*log(2*pi)
		ll
	}

	grad.sigma = function(sigma,y,mu,...){

		y_U = y[y>0]
		mu_U = mu[y>0]
		y_C = y[y==0]
		mu_C = mu[y==0]

		n_U = length(y_U)

		grad.sigma = 0
		grad.sigma = grad.sigma + 1/sigma^2*sum(mu_C*dnorm(mu_C/sigma)/pnorm(-mu_C/sigma))
		grad.sigma = grad.sigma + 1/sigma^3*sum((y_U-mu_U)^2) - n_U*1/sigma

		grad.sigma
	}

	scores.sigma = function(sigma,y,mu){
		y_U = y[y>0]
		mu_U = mu[y>0]
		y_C = y[y==0]
		mu_C = mu[y==0]

		score.sigma = rep(NA,length(y))
		score.sigma[y==0] = mu_C/sigma^2*IHR(mu_C/sigma)
		score.sigma[y>0] = 1/sigma^3*(y_U-mu_U)^2 - 1/sigma

		score.sigma
	}

	hess.sigma = function(sigma,y,mu){
		y_U = y[y>0]
		mu_U = mu[y>0]
		y_C = y[y==0]
		mu_C = mu[y==0]

		coef_C = IHR(mu_C/sigma)

		hess.sigma = 0
		hess.sigma = hess.sigma + sum( mu_C/sigma^3*coef_C*(-2+mu_C^2/sigma^2-mu_C/sigma*coef_C) )
		hess.sigma = hess.sigma + sum( -3/sigma^4*(y_U-mu_U)^2 + 1/sigma^2 )

		hess.sigma
	}

	hess.sigmaL = function(sigma,jacob.mat,y,mu,dxi_dbeta,ll_d2){
		#if there's no dummies, then dxi_dbeta = 0
		y_U = y[y>0]
		mu_U = mu[y>0]
		y_C = y[y==0]
		mu_C = mu[y==0]

		mu_s_C = mu_C/sigma
		coef_C = IHR(mu_s_C)/sigma

		ll_dxdt = crossprod(jacob.mat[y==0,],coef_C/sigma*(1-mu_s_C^2+mu_C*coef_C)) - 2*crossprod(jacob.mat[y>0,],(y_U-mu_U)/sigma^3)

		d_mu = jacob.mat + dxi_dbeta
		H = crossprod(d_mu,dxi_dbeta * ll_d2 + ll_dxdt)

		return(as.matrix(H))
	}

	expected.predictor = function(mu,env,coef,...){
		#sigma = get(".sigma",env)
		sigma = coef[".sigma"]
		sigma * dnorm(mu/sigma) + pnorm(mu/sigma) * mu
	}

	initDummy = function(S,y,mu,env,...){
		stop("not yet implemented")
		#guess for the dummy:
		#we assume mu is at the average level
		ni1 = as.vector(S%*%y)
		ni = rowSums(S)
		qnorm(ni1/ni)-as.vector(S%*%mu)/ni
	}

	guessDummy = function(sum_y,n_group,mu,...){
		stop("not yet implemented")
		#Dummy when mu is constant within the group
		qnorm(sum_y/n_group)-mu
	}

	dum_fx = function(y,x,sum_y,S,mu,dum,env,...){
		#On calcule la condition de premier ordre des dummies
		# ATTENTION PAS IMPLEMENTE CORRECTEMENT
		stop("not yet implemented")
		sigma = get(".sigma",env)
		id_U = which(y>0)
		id_C = which(y==0)

		mu = mu + x[dum]

		y_U = y[id_U]
		mu_U = mu[id_U]
		y_C = y[id_C]
		mu_C = mu[id_C]

		mu_s_C = mu_C/sigma
		coef_C = IHR(mu_s_C)

		cpo = rep(NA,length(y))
		cpo[id_C] = - coef_C
		cpo[id_U] = (y_U - mu_U)/sigma^2

		as.vector(S%*%cpo)
	}

	dum_dfx = function(x,S,mu,dum,y,coef,env,...){
		#On calcule la derivee des dummies
		sigma = get(".sigma",env)
		id_U = which(y>0)
		id_C = which(y==0)

		mu = mu + x[dum]

		y_C = y[id_C]
		mu_C = mu[id_C]

		mu_s_C = mu_C/sigma
		coef_C = IHR(mu_s_C)

		cso = rep(NA,length(y))
		cso[id_C] = coef_C/sigma^2 * (mu_s_C-coef_C)
		cso[id_U] = -1/sigma^2

		as.vector(S%*%cso)
	}

	ratio_fx_dfx = function(x,dum,S,y,mu,env,...){
		#On calcule la cond. de premier ordre, celle de second
		#ordre, et alors le pas optimal
		sigma = get(".sigma",env)
		id_U = which(y>0)
		id_C = which(y==0)

		mu = mu + x[dum]

		y_U = y[id_U]
		mu_U = mu[id_U]
		y_C = y[id_C]
		mu_C = mu[id_C]

		mu_s_C = mu_C/sigma
		coef_C = IHR(mu_s_C)

		cpo = rep(NA,length(y))
		cpo[id_C] = - coef_C
		cpo[id_U] = (y_U - mu_U)/sigma^2

		Value = as.vector(S%*%cpo)

		cso = rep(NA,length(y))
		cso[id_C] = coef_C/sigma^2 * (mu_s_C-coef_C)
		cso[id_U] = -1/sigma^2

		Derivee = as.vector(S%*%cso)

		Value / Derivee
	}

	ll0 = function(coef,y){
		#La fonction minimise, on renvoie "-"
		sigma = coef[2]
		cste = coef[1]

		y_U = y[y>0]
		n_U = length(y_U)
		n_C = sum(y==0)

		ll = n_C*pnorm(-cste/sigma,log.p = TRUE) - 1/2/sigma^2*sum((y_U-cste)^2) - n_U*log(sigma) - n_U/2*log(2*pi)

		-ll
	}

	grad0 = function(coef,y){
		#La fonction minimise, on renvoie "-"
		sigma = coef[2]
		cste = coef[1]

		y_U = y[y>0]
		n_U = length(y_U)
		n_C = sum(y==0)

		grad.cste = 0
		grad.cste = grad.cste + -n_C/sigma*IHR(cste/sigma)
		grad.cste = grad.cste + 1/sigma^2*sum(y_U-cste)

		grad.sigma = 0
		grad.sigma = grad.sigma + n_C*cste/sigma^2*IHR(cste/sigma)
		grad.sigma = grad.sigma + 1/sigma^3*sum((y_U-cste)^2) - n_U*1/sigma

		- c(grad.cste,grad.sigma)
	}

	model0 = function(y){
		stop("not yet implemented")
		n = length(y)
		x = log(sum(y)) - log(length(y) - sum(y))
		ll = x*sum(y) - n*log(1+exp(x))
		list(loglik=ll,constant=x)
	}

	return(list(ll=ll,expected.predictor=expected.predictor,ll0=ll0,grad0=grad0,ll_TEST_score=ll_TEST_score,model0=model0,ll_dl=ll_dl,ll_d2=ll_d2,ratio_fx_dfx=ratio_fx_dfx,initDummy=initDummy,guessDummy=guessDummy,dum_fx=dum_fx,dum_dfx=dum_dfx))
}

