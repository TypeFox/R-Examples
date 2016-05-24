
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coxph_mpl=function(formula,data,subset,na.action,control,...){
	#  get and organise information 
	# (same tests as in coxph(), thanks to the survival package)
	mc = match.call(expand.dots = FALSE)
	m  = match(c("formula","data","subset","na.action"),names(mc),0)
	mc = mc[c(1,m)]	
	if (m[1]==0){stop("A formula argument is required")}
	data.name = if(m[2]!=0){mc[m[2]][[1]]}else{"-"}
	mc[[1]] = as.name("model.frame")
	mc$formula = if(missing(data)) terms(formula)
                     else              terms(formula, data=data)
	mf = eval(mc,parent.frame())
	if (any(is.na(mf))) stop("Missing observations in the model variables")
	if (nrow(mf) ==0) stop("No (non-missing) observations")
	mt = attr(mf,"terms")
	# extract response
	y    = model.extract(mf, "response")
	type = attr(y, "type")
	if(!inherits(y, "Surv")){stop("Response must be a survival object")}
	if(type!="right"&&type!="counting"){
		stop(paste("Cox model doesn't support \"", type, "\" survival data",sep = ""))
		}
	t_i         = y[,1L]
	observed    = y[,2L]==1L
	n           = length(t_i)
	n.obs       = sum(observed)	
	# control arguments
	extraArgs <- list(...)
	if (length(extraArgs)) {
	controlargs <- names(formals(coxph_mpl.control)) 
	m <- pmatch(names(extraArgs), controlargs, nomatch=0L)
	if (any(m==0L))
	    stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m==0L]),
		 domain = NA, call. = FALSE)
	}	
	if (missing(control)) control <- coxph_mpl.control(n.obs, ...)
	# ties
	t_i.obs  = t_i[observed]	
	ties     = duplicated(t_i.obs)
	if(any(ties)){
		if(control$ties=="epsilon"){
			if(length(control$seed)>0){
				old <- .Random.seed
				on.exit({.Random.seed <<- old})
				set.seed(control$seed)
			}
			t_i.obs[ties] = t_i.obs[ties]+runif(sum(ties),-1e-11,1e-11)
			t_i[observed] = t_i.obs
		}else{	
			t_i.obs = t_i.obs[!ties]
			n.obs   = length(t_i.obs)
		}
	}
	# X and centered X matrix
	X           = model.matrix(mt, mf, contrasts)
	X           = X[,!apply(X, 2, function(x) all(x==x[1])), drop=FALSE]
	if(ncol(X)==0){
		X   = matrix(0,n,1)
		noX = TRUE
	}else{  noX = FALSE}
	p           = ncol(X)	
	mean_j      = apply(X, 2, mean)	
	XC          = X - rep(mean_j, each=n)	
	# knot sequence and psi matrices
	knots  = knots_mpl(control, t_i.obs, range(t_i))
	m      = knots$m
	psi    = basis_mpl(t_i,knots,control$basis,control$order,which=1)
	PSI    = basis_mpl(t_i,knots,control$basis,control$order,which=2)
	R      = penalty_mpl(control,knots)
	R_star = rbind(matrix(0,p,p+m),cbind(matrix(0,m,p),R))	
	# God save the Queen
	lambda     = control$smooth
	Beta       = rep(0,p)
	Theta      = rep(1,knots$m)
	correction = 1
	full.iter  = 0
	for(iter in 1:control$max.iter[1]){
		fit <- .Call("coxph_mpl",
				status2 = as.integer(observed), X2=XC, meanX2 = mean_j, R2 = R,
				psi2  = psi, PSI2  = PSI,
				Beta0 = Beta, Theta0 = Theta/correction, 
				lambda2 = as.double(lambda), kappa2 = control$kappa,
				convVal2 = control$tol, minTheta2 = control$epsilon, maxiter2 = control$max.iter[2])			
		## sigma
		correction = fit$ploglik[2]
		M_theta_m1 = Theta = fit$coef$Theta
		Beta       = fit$coef$Beta
		H          = fit$matrices$H
		HRinv      = matrix(0,p+m,p+m)	
		sigma2_old = 1/(2*lambda)
		full.iter  = full.iter+fit$iter
		#
		pos   = c(if(noX){FALSE}else{rep(TRUE,p)},M_theta_m1>control$min.theta)
		HRinv[pos,pos] = solve(H[pos,pos]+(1/sigma2_old)*R_star[pos,pos])
		# update
		sigma2     = c(t(M_theta_m1)%*%R%*%M_theta_m1/
			      (m-(sum(diag(HRinv%*%R_star))/sigma2_old)))	
		lambda_old = ifelse(iter>1,lambda,-1)
		lambda     = 1/(2*sigma2)
		# check for convergence
		cat(iter,"\tlambda =",lambda,"\titer =",fit$iter,"\n")	
		if((abs(lambda-lambda_old)/abs(lambda_old))<(control$tol)){break}
		}
	### God save the King		
	if(control$max.iter[1]>1) control$smooth = lambda 
	M_theta_m1 = fit$coef$Theta
	M_1    = fit$matrices$M1
	M_2    = fit$matrices$M2
	H      = fit$matrices$H
	Q      = fit$matrices$Q
	p      = length(fit$coef$Beta)
	m      = length(fit$coef$Theta)
	pos    = c(if(noX){FALSE}else{rep(TRUE,p)},M_theta_m1>control$min.theta)
	Minv_1 = Minv_2 = Hinv = matrix(0,p+m,p+m) 
	# correction factor and correction matrix (jacobian from the delta method)
	M_corr_mpm = cbind(matrix(fit$ploglik[2]*(rep(fit$coef$Theta,p)*rep(-fit$coef$Beta,each=m)),ncol=p),
	                   diag(rep(fit$ploglik[2],m)))
	# M1QM1 and M1HM1
	temp = try(solve(M_1[pos,pos]),silent=TRUE) 
	if(class(temp)!="try-error"){
		Minv_1[pos,pos] = temp
		cov_NuNu_M1QM1  = Minv_1%*%Q%*%Minv_1
		cov_NuNu_M1HM1  = Minv_1%*%H%*%Minv_1
		cov_NuNu_M1QM1[-(1:p),-(1:p)] = M_corr_mpm%*%cov_NuNu_M1QM1%*%t(M_corr_mpm)
		cov_NuNu_M1HM1[-(1:p),-(1:p)] = M_corr_mpm%*%cov_NuNu_M1HM1%*%t(M_corr_mpm)		
		se.Eta_M1QM1    = suppressWarnings(sqrt(diag(cov_NuNu_M1QM1)))
		se.Eta_M1HM1    = suppressWarnings(sqrt(diag(cov_NuNu_M1HM1)))
	}else{
		cov_NuNu_M1QM1  = cov_NuNu_M1HM1 = matrix(NA,p+m,p+m)
		se.Eta_M1QM1    = se.Eta_M1HM1   = rep(NA,p+m)
	}
	# M2QM2 and M2HM2
	temp = try(solve(M_2[pos,pos]),silent=TRUE) 
	if(class(temp)!="try-error"){
		Minv_2[pos,pos] = temp		
		cov_NuNu_M2QM2  = Minv_2%*%Q%*%Minv_2
		cov_NuNu_M2HM2  = Minv_2%*%H%*%Minv_2
		cov_NuNu_M2QM2[-(1:p),-(1:p)] = M_corr_mpm%*%cov_NuNu_M2QM2%*%t(M_corr_mpm)
		cov_NuNu_M2HM2[-(1:p),-(1:p)] = M_corr_mpm%*%cov_NuNu_M2HM2%*%t(M_corr_mpm)				
		se.Eta_M2QM2    = suppressWarnings(sqrt(diag(cov_NuNu_M2QM2)))
		se.Eta_M2HM2    = suppressWarnings(sqrt(diag(cov_NuNu_M2HM2)))
	}else{
		cov_NuNu_M2QM2  = cov_NuNu_M2HM2 = matrix(NA,p+m,p+m)
		se.Eta_M2QM2    = se.Eta_M2HM2   = rep(NA,p+m)
	}
	# Hessian
	temp = try(solve(H[pos,pos]),silent=TRUE) 
	if(class(temp)!="try-error"){
		Hinv[pos,pos] = temp
		cov_NuNu_H    = Hinv
		cov_NuNu_H[-(1:p),-(1:p)] = M_corr_mpm%*%cov_NuNu_H%*%t(M_corr_mpm)						
		se.Eta_H      = suppressWarnings(sqrt(diag(cov_NuNu_H)))
	}else{
		cov_NuNu_H  = matrix(NA,p+m,p+m)
		se.Eta_H    = rep(NA,p+m)
	}
	## out
	mx.seNu.l5=as.data.frame(cbind(se.Eta_M1QM1,se.Eta_M1HM1,se.Eta_M2QM2,se.Eta_M2HM2,se.Eta_H))
	colnames(mx.seNu.l5)=c("M1QM1","M1HM1","M2QM2","M2HM2","H")
	rownames(mx.seNu.l5)[(p+1):(p+m)] = paste("Theta",1:m,sep="")		
	rownames(mx.seNu.l5)[(1:p)] = paste("Beta",1:p,sep="")
	fit$se    = list(Beta=mx.seNu.l5[1:p,],Theta=mx.seNu.l5[(p+1):(p+m),])
	fit$covar = list(M1QM1=cov_NuNu_M1QM1,M1HM1=cov_NuNu_M1HM1,
			 M2QM2=cov_NuNu_M2QM2,M2HM2=cov_NuNu_M2HM2,H=cov_NuNu_H)  	
	### output
	fit         = fit[is.na(match(names(fit),c("flag","matrices")))]
	fit$knots   = knots
	fit$control = control
	fit$call    = match.call()
	fit$dim     = list(n = n, n.obs = sum(observed), n.ties = sum(ties), p = p, m = knots$m)
	fit$data    = list(time = t_i, observed = observed, X = X, name = data.name)
	fit$iter    = c(iter,full.iter,fit$iter)
	fit$loglik  = fit$ploglik[3]
	class(fit)  = "coxph_mpl"
	fit
	}
	
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print.coxph_mpl=function(x, ...) {
	cat("\n")
	print(x$call)
	cat("\nLog-likelihood : ",x$ploglik[1],"\n",sep="")	
	cat("\nRegression parameters :\n")
	vect=c(x$coef$Beta)
	names(vect)=dimnames(x$data$X)[[2]]
	print(vect, ...)
	cat("\nBaseline hasard parameters : \n")
	print(x$coef$Theta, ...)
	cat("\n")
	}

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary.coxph_mpl=function(object,se="M2QM2",full=FALSE,...) {
	col.names = c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
	seB   = object$se$Beta[[se]]
	matxB = cbind(object$coef$Beta,seB,object$coef$Beta/seB,2*(1-pnorm(abs(object$coef$Beta/seB))))
	dimnames(matxB)=list(paste(" ",dimnames(object$data$X)[[2]],sep=""),col.names)
	pos = object$coef$Theta>object$control$min.theta
	if(full){
		seT   = object$se$Theta[[se]]
		matxT = cbind(object$coef$Theta,seT,object$coef$Theta/seT,2*(1-pnorm(abs(object$coef$Theta/seT))))
		dimnames(matxT)=list(format(1:object$dim$m,just="right"),col.names)
		matxT = matxT[pos,]
	}else{	matxT = object$coef$Theta[pos]
	        names(matxT) = seq(1,object$dim$m)[pos]}
	out = list(Beta = matxB, Theta = matxT, inf = list(iter = object$iter, call=object$call, data = object$data$name,
		control = object$control, dim = object$dim, trunc = trunc, full = full, ploglik = object$ploglik[1]))
	class(out) = "summary.coxph_mpl"
	out
}

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
coxph_mpl.control <- function(n.obs=NULL, basis = "uniform", smooth = NULL, max.iter=c(50,1e+05), tol=1e-7, 
		n.knots = NULL, n.events_basis = NULL, range.quant = c(0.075,.9),
		cover.sigma.quant = .25, cover.sigma.fixed=.25, min.theta = 1e-10,
		penalty = 2L, order = 3L, kappa = 1/.6, epsilon = 1e-50, ties = "unique", seed = NULL){
	basis       = basis.name_mpl(basis)
	max.iter     = c(ifelse(is.null(smooth),ifelse(max.iter[1]>0,as.integer(max.iter[1]),50),1L),
	                 ifelse(max.iter[2]>0,as.integer(max.iter[2]),1e+05))	
	tol          = ifelse(tol>0 & tol<1,tol,1e-7)	
	order        = ifelse(order>0 & order<6,as.integer(order),3L)	
	min.theta    = ifelse(min.theta>0 & min.theta<1e-3,min.theta,1e-10)	
	penalty      = penalty.order_mpl(penalty,basis,order)
	kappa        = ifelse(kappa>1, kappa, 1/.6)
	cover.sigma.quant  = ifelse(cover.sigma.quant>0 & cover.sigma.quant<0.4,cover.sigma.quant,0.3)
	cover.sigma.fixed  = ifelse(cover.sigma.fixed>0 & cover.sigma.fixed<0.4,cover.sigma.fixed,0.3)	
	if(all(range.quant<1) & all(range.quant>0) & length(range.quant)==2){
		range.quant = range.quant[order(range.quant)]
		}else{range.quant = c(0.075,.9)}
	if(is.null(n.knots)|sum(n.knots)<3|length(n.knots)!=2){
		n.knots    = if(basis!='uniform' & basis!='msplines'){c(0,20)}else{c(8,2)}
	}
	if(!is.null(n.events_basis)){
		n.events_basis = ifelse(n.events_basis<1|n.events_basis>floor(n.obs/2),
			max(round(3.5*log(n.obs)-7.5),1L),round(n.events_basis))
		}else{n.events_basis = max(round(3.5*log(n.obs)-7.5),1L)}	
	if(!is.null(smooth)){
		smooth = ifelse(smooth<0,0,smooth)
		}else{smooth=0}
	out = list(basis = basis, smooth = smooth, max.iter = max.iter, tol = tol,
		order = order, penalty = penalty, n.knots = n.knots, range.quant = range.quant,
		cover.sigma.quant = cover.sigma.quant, cover.sigma.fixed = cover.sigma.fixed,
		n.events_basis = as.integer(n.events_basis), min.theta = min.theta, ties = ties,
		seed = as.integer(seed), kappa = kappa, epsilon = epsilon)
	class(out) = "coxph_mpl.control"
	out
	}


	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
basis.name_mpl <- function(k){
	if(k == "discr"| k == "discretized" | k == "discretised" | k == "unif" | k == "uniform"){"uniform"
	}else{if(k == "m" | k == "msplines" | k == "mspline"){"msplines"
	}else{if(k == "gauss" | k == "gaussian"){"gaussian"
	}else{if(k == "epa" | k == "epanechikov"){"epanechikov"
	}else{stop("Unkown basis choice", call. = FALSE)}}}}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
penalty.order_mpl <- function(p,basis,order){
	p = as.integer(p)
	switch(basis,
	'uniform'  = ifelse(p>0 & p<3,p,2),
	'gaussian' = ifelse(p>0 & p<3,p,2),
	'msplines' = order-1,
	'epa'      = 2)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knots_mpl=function(control,events,range){
	n.events  = length(events)
	## uniform
	if(control$basis=="uniform"){
		m          = floor(n.events/control$n.events_basis) 
		n.i_u      = round(control$n.events_basis*c(rep(1L,m-1L),1L+n.events/control$n.events_basis-m))
		Alpha      = c(range[1],unlist(lapply(split(events[order(events)],rep(1L:m,n.i_u)),max)))
		if(Alpha[length(Alpha)]<range[2]){Alpha[length(Alpha)]=range[2]}		
		Alpha      = unique(Alpha)
		m          = length(Alpha)-1
		Delta      = Alpha[2L:(m+1L)]-Alpha[1L:m]
		list(m=m,Alpha=Alpha,Delta=Delta)
	## other
	}else{
		# quantile knots
		if(control$n.knots[1]>0){
			Alpha1 = quantile(events,prob=seq(control$range.quant[1],control$range.quant[2],length=control$n.knots[1]))
		}else{Alpha1 = NULL}
		control$n.knots[2] = max(control$n.knots[2]+2,2)
		Alpha2  = c(range[1],seq(ifelse(control$n.knots[1]==0,range[1],max(Alpha1)),
			    range[2],length=control$n.knots[2])[-1])
		Alpha   = unique(c(Alpha1,Alpha2))
		Alpha   = Alpha[order(Alpha)]
		n.Alpha = length(Alpha)
		## gaussian-basis
		if(control$basis=="gaussian"){
			Sigma = Delta = rep(0,n.Alpha)
			for(aw in 1:n.Alpha){
				if(aw>1 & aw<(n.Alpha-control$n.knots[2]+2)){
					while(sum(events>(Alpha[aw]-2*Sigma[aw])&events<(Alpha[aw]+2*Sigma[aw]))<(n.events*control$cover.sigma.quant)){
						Sigma[aw] = Sigma[aw] + 0.001}
				}else{Sigma[aw] = control$cover.sigma.fixed*(Alpha[n.Alpha]-Alpha[1])/3}
				Delta[aw]= pnorm((range[2]-Alpha[aw])/Sigma[aw])-
				           pnorm((range[1]-control$epsilon-Alpha[aw])/Sigma[aw])
			}
			list(m=n.Alpha, Alpha=Alpha, Sigma=Sigma, Delta=Delta)
		## m-splines and epanechikov
		}else{
			m = n.Alpha+control$order-2
			list(m=m, Alpha=Alpha, Delta=rep(1,m))
		}
	}}
	

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
basis_mpl = function(x,knots,basis,order,which=c(1,2)){
	which.matrix = rep(T,2)
	which.matrix[-which]=FALSE	
	n        = length(x)
	Alpha    = knots$Alpha
	Delta    = knots$Delta
	n.Alpha  = length(Alpha)
	m        = ifelse(basis=="msplines"|basis=="epanechikov",n.Alpha+order-2,knots$m)
	M_Psi_nm = M_psi_nm = matrix(0,n,m)
	##
	if(basis=="uniform"){
		u_i = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)
		for(i in 1:n){
			M_psi_nm[i,u_i[i]]   = 1
			M_Psi_nm[i,1:u_i[i]] = c(if(u_i[i]>1){Delta[1:(u_i[i]-1)]},
			                         x[i]-Alpha[u_i[i]])
			}
	##
	}else{
	if(basis=="gaussian"){
		Sigma = knots$Sigma
		for(u in 1:m){
			M_psi_nm[,u] =  dnorm((x-Alpha[u])/Sigma[u])/(Sigma[u]*Delta[u])
			M_Psi_nm[,u] = (pnorm((x-Alpha[u])/Sigma[u])-
				        pnorm((Alpha[1]-Alpha[u])/Sigma[u]))/Delta[u]
			}
	##
	}else{
		seq1n = 1:n
		Alpha_star   = as.numeric(c(rep(Alpha[1],order-1L),Alpha,rep(Alpha[n.Alpha],order-1L)))
		M_psi_nm     = M_Psi_nm = cbind(M_psi_nm,0)		
		if(which.matrix[1]){
			Alpha_star_x = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L
			if(basis=="msplines"){
				M_psi_nm[(Alpha_star_x-1L)*n+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
				if(order>1){
				for(ow in 2L:order){
					uw_x = Alpha_star_x-ow+1L
					for(pw in 0:(ow-1L)){
						pos_x = (uw_x+pw-1L)*n+seq1n
						M_psi_nm[pos_x]=(ow/((ow-1)*(Alpha_star[1:m+ow]-Alpha_star[1:m])))[uw_x+pw]*
							       ((x-Alpha_star[uw_x+pw])*M_psi_nm[pos_x]+
								(Alpha_star[uw_x+pw+ow]-x)*M_psi_nm[pos_x+n])	
						}
					}
				}
			# Epanechikov
			}else{
				uw_x = Alpha_star_x-order+1L
				for(pw in 0:(order-1L)){
					pos_x = (uw_x+pw-1L)*n+seq1n
					pos_1 = (uw_x+pw)==1
					pos_m = (uw_x+pw)==m
					pos_other = pos_1==FALSE & pos_m==FALSE	
					# 1<u<m
					M_psi_nm[pos_x[pos_other]]=(6*(x-Alpha_star[uw_x+pw])*(x-Alpha_star[uw_x+pw+order])/ 
						       ((Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_other]
					# case u=1
					M_psi_nm[pos_x[pos_1]]=(12*(x-Alpha_star[uw_x+pw+order])*(x-2*Alpha_star[uw_x+pw]+Alpha_star[uw_x+pw+order])/ 
						       ((2*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])^3))[pos_1]
					# case u=m
					M_psi_nm[pos_x[pos_m]]=(12*(x-Alpha_star[uw_x+pw])*(x+Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])/ 
						       ((2*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])^3))[pos_m]
					}
			}
			M_psi_nm = M_psi_nm[,1:m,drop=FALSE]
		}	
		if(which.matrix[2]){
			rank.x   = rank(x)
			x        = x[order(x)]	
			Alpha_x  = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)	
			# integral equals 1
			up_u     = cumsum(tabulate(Alpha_x,n.Alpha-1))
			for(uw in 1:(m-order+1)){M_Psi_nm[min(n,up_u[uw]+1):n,uw] = 1}	
			# other cases
			if(basis=="msplines"){
				Alpha_star2 = c(rep(Alpha[1],order),Alpha,rep(Alpha[n.Alpha],order))	
				factor_v    = c((Alpha_star2[(order+2):length(Alpha_star2)]-Alpha_star2[1:(length(Alpha_star2)-order-1)])/
					       (order+1),rep(0,order-1))
				M_psi2_nm   = cbind(basis_mpl(x,knots,basis=basis,order=order+1,which=1),matrix(0,n,order-1))
				pos_xo  = rep((Alpha_x-1L)*n,1)+seq1n
				pos_xo1 = rep(pos_xo,order)+rep(1:order,each=n)*n
				for(ow in 0:(order-1)){
					M_Psi_nm[pos_xo+ow*n] = apply(matrix(M_psi2_nm[pos_xo1+ow*n]*
						 factor_v[rep(Alpha_x,order)+rep((1:order)+ow,each=n)],ncol=order),1,sum)
					}
			}else{
				Alpha_star_x = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L
				uw_x = Alpha_star_x-order+1L
				for(pw in 0:(order-1L)){
					pos_x = (uw_x+pw-1L)*n+seq1n
					pos_1 = (uw_x+pw)==1
					pos_m = (uw_x+pw)==m
					pos_other = pos_1==FALSE & pos_m==FALSE	
					# 1<u<m
					M_Psi_nm[pos_x[pos_other]]=((x-Alpha_star[uw_x+pw])^2*(2*x+Alpha_star[uw_x+pw]-3*Alpha_star[uw_x+pw+order])/ 
						       ((Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_other]
					# case u=1
					M_Psi_nm[pos_x[pos_1]]=((x-Alpha_star[uw_x+pw])*(x^2-2*x*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw]^2+
								6*Alpha_star[uw_x+pw]*Alpha_star[uw_x+pw+order]-3*Alpha_star[uw_x+pw+order]^2)/ 
						       (2*(Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_1]
					# case u=m
					M_Psi_nm[pos_x[pos_m]]=((x-Alpha_star[uw_x+pw])^2*(x+2*Alpha_star[uw_x+pw]-3*Alpha_star[uw_x+pw+order])/ 
						       (2*(Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_m]
					}
			}
			M_Psi_nm = M_Psi_nm[rank.x,1:m,drop=FALSE]
	}}}
	# pdf and cdf
	if(all(which.matrix)){list(psi=M_psi_nm,Psi=M_Psi_nm)
	# pdf or cdf
	}else{if(which.matrix[1]){M_psi_nm}else{M_Psi_nm}}
	}

		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
penalty_mpl=function(control,knots){
	#
	penalty = control$penalty
	m       = knots$m
	## uniform
	if(control$basis=="uniform"){
		D = diag(m)*c(1,-2)[penalty]
		E = diag(m+1)[-(m+1),-1]*c(-1,1)[penalty]
		B = D+E+list(0,t(E))[[penalty]]
		B[c(m+1,(m-1)*m)[1:penalty]] = c(-1,2)[penalty]
		M_R_mm = t(B)%*%B		
	## gaussian
	}else{
		M_R_mm = matrix(0,m,m)
		if(control$basis=="gaussian"){
			# useful function
			int_rij_2.fun=function(x,mu_i,mu_j,sig_i,sig_j,t1,tn){
				K     = 4*(pnorm((t1-mu_i)/sig_i)-pnorm((tn-mu_i)/sig_i))*(pnorm((t1-mu_j)/sig_j)-pnorm((tn-mu_j)/sig_j))
				q1q2a = 4*dnorm(x,mu_i,sig_i)*dnorm(x,mu_j,sig_j)*sig_i*sig_j*2*pi*sqrt(sig_i^2+sig_j^2)*
					((mu_j-x)*sig_i^6*((x-mu_i)^2-sig_i^2)+
					  sig_i^4*sig_j^2*((x-4*mu_j+3*mu_i)*sig_i^2-(x-mu_i)^2*(3*x-4*mu_j+mu_i))-
					  sig_i^2*sig_j^4*(x-mu_j)^2*(3*x-4*mu_i+mu_j)+
					  sig_j^6*((x+3*mu_j-4*mu_i)*sig_i^2-(x-mu_j)^2*(x-mu_i))+
					  sig_j^8*(x-mu_i)
					)
				q1q2b = 2*dnorm(mu_j,mu_i,sqrt(sig_i^2+sig_j^2))*pi*sqrt(sig_i^2+sig_j^2)*sig_i^3*sig_j^3*
					(mu_j^4-4*mu_j^3*mu_i+mu_i^4+6*mu_j^2*(mu_i^2-sig_i^2-sig_j^2)-
					 6*mu_i^2*(sig_i^2+sig_j^2)+3*(sig_i^2+sig_j^2)^2-4*mu_i*mu_j*(mu_i^2-3*(sig_i^2+sig_j^2))
					)*
					(2*pnorm(((x-mu_j)*sig_i^2+(x-mu_i)*sig_j^2)/(sig_i*sig_j*sqrt(sig_i^2+sig_j^2)))-1)
				q3    = pi*sig_i^3*sig_j^3*(sig_i^2+sig_j^2)^(9/2)*K
				(q1q2a+q1q2b)/q3
				}
			int_rij_1.fun=function(x,mu_i,mu_j,sig_i,sig_j,t1,tn){
				K  = 4*(pnorm((t1-mu_i)/sig_i)-pnorm((tn-mu_i)/sig_i))*(pnorm((t1-mu_j)/sig_j)-pnorm((tn-mu_j)/sig_j))
				q2 =  dnorm((mu_i-mu_j)/sqrt(sig_i^2+sig_j^2))*2*pi*
				      sig_i*sig_j*(sig_i^2+sig_j^2-(mu_i-mu_j)^2)*
				      (2*pnorm(((x-mu_j)*sig_i^2+(x-mu_i)*sig_j^2)/(sig_i*sig_j*sqrt(sig_i^2+sig_j^2)))-1)
				q1 =   4*pi*dnorm((x-mu_i)/sig_i)*dnorm((x-mu_j)/sig_j)*
				       sqrt(sig_i^2+sig_j^2)*((mu_i-x)*sig_i^2+(mu_j-x)*sig_j^2)
				q3 =  pi*sig_i*sig_j*(sig_i^2+sig_j^2)^(5/2)*K
			       (q1+q2)/q3	
			}
			# fill
			for(i in 1:m){
				for(j in i:m){
					if(penalty==2){
						M_R_mm[i,j] = M_R_mm[j,i] = 
						int_rij_2.fun(knots$Alpha[m],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])-
						int_rij_2.fun(knots$Alpha[1],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])
					}else{  M_R_mm[i,j] = M_R_mm[j,i] = 
						int_rij_1.fun(knots$Alpha[m],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])-
						int_rij_1.fun(knots$Alpha[1],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])}
				}
			}
		## other
		}else{  Alpha      = knots$Alpha
			n.Alpha    = length(Alpha)
			order      = control$order
			Alpha_star = c(rep(Alpha[1],order-1L),Alpha,rep(Alpha[n.Alpha],order-1L))
			# msplines
			if(control$basis=="msplines"){
				seq1n        = 1L:(n.Alpha-1)
				n.Alpha_star = length(Alpha_star)
				Alpha_star_x = sapply(Alpha[-1],function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L				
				## 2derivative of each u,v pair
				M_d2f_mm = matrix(0,n.Alpha-1,n.Alpha+order-1L)
				M_d2f_mm[(Alpha_star_x-1L)*(n.Alpha-1)+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
				for(ow in 2L:order){
					pw   = 1L:ow 
					uw_x = Alpha_star_x-ow+1L
					for(pw in 0:(ow-1L)){
						M_d2f_mm[(uw_x+pw-1L)*(n.Alpha-1)+seq1n]=
							(ow/(Alpha_star[1:(n.Alpha+ow)+ow]-Alpha_star[1:(n.Alpha+ow)]))[uw_x+pw]*
							(M_d2f_mm[(uw_x+pw-1L)*(n.Alpha-1)+seq1n]-M_d2f_mm[(uw_x+pw)*(n.Alpha-1)+seq1n])
					}
				}
				M_d2f_mm = M_d2f_mm[,1:m,drop=FALSE]
				for(uw in 1:m){
					for(vw in uw:m){
						M_R_mm[uw,vw] = M_R_mm[vw,uw] = 
							sum((M_d2f_mm[,uw]*M_d2f_mm[,vw])*(Alpha[-1]-Alpha[-n.Alpha]))
					}
				}
			# epanechikov
			}else{	for(uw in 1:m){
					f_u = ifelse(uw==1|uw==m,4,1)
					for(vw in uw:m){
						if(Alpha_star[vw]<Alpha_star[uw+order]){
							f_v = ifelse(vw==1|vw==m,4,1)
							M_R_mm[uw,vw] = M_R_mm[vw,uw] = (144*(Alpha_star[uw+order]-Alpha_star[vw]))/
								(f_v*f_u*(Alpha_star[uw+order]-Alpha_star[uw])^3*(Alpha_star[vw+order]-Alpha_star[vw])^3)
						}
					}
				}
			}
	}}
	M_R_mm
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.coxph_mpl=function(x,se="M2QM2",ask=TRUE,which=1:4,upper.quantile=.95,...){
	which.plot=rep(TRUE,4)
	if(!is.null(which)){which.plot[-which]=FALSE}
	if(sum(which.plot)==1){ask=FALSE}	
	if(ask){oask <- devAskNewPage(TRUE)
		on.exit(devAskNewPage(oask))
		}
	control = x$control
	knots   = x$knots
	pos     = x$coef$Theta<x$control$min.theta
	# x axis points
	n.x       = 1000
	V_x_X     = seq(x$knots$Alpha[1],max(x$knots$Alpha),length=n.x)	
	colw      = terrain.colors(x$dim$m+1)
	prob      = upper.quantile
	quant     = quantile(x$data$time,prob=prob)
	M_psi_Xm  = basis_mpl(V_x_X,knots,control$basis,control$order,which=1)	
	# Bases
	if(which.plot[1]){
		plot(1,1,pch="",xlim=range(V_x_X),ylim=max(M_psi_Xm)*c(-.05,1),axes=FALSE,
		     xlab="Survival time",ylab=expression(psi[u]^{o}*(t)),
		     main=paste(if(control$basis=="uniform"){"Uniform"}else{
				if(control$basis=="gaussian"){"Gaussian"}else{
				if(control$basis=="msplines"){"M-spline"}else{
				if(control$basis=="epanechikov"){"Epanechikov"}}}},
				" bases used to approximate the baseline hazard\n",
				"(",x$dim$m," bases)",sep=""),...)
		abline(v=knots$Alpha,col=gray(.9))
		abline(v=quant,col=gray(.5),lty=2)
		for(u in 1:x$dim$m){lines(V_x_X,M_psi_Xm[,u],col=colw[u],lty=2*pos[u]+1)}
		axis(2,las=2,...)
		points(x$data$time,jitter(rep(-.05*max(M_psi_Xm),x$dim$n),10),col=(!x$data$observed)+1,pch=1)
		axis(1,...)
		if(any(pos)){
		legend("topright",legend=c("Observed","Censored","Knots",as.expression(bquote(widehat(F)^-1*(.(prob)))),
		       expression(hat(theta)[u]==0),expression(hat(theta)[u]>0)),
		       pch=c(1,1,NA,NA,NA,NA),lty=c(0,0,1,2,3,1),col=c("black","red",gray(.9),gray(.5),colw[1],colw[1]),ncol=1,bty="n",cex=.75)
		}else{
		legend("topright",legend=c("Observed","Censored","Knots",as.expression(bquote(widehat(F)^-1*(.(prob))))),
		       pch=c(1,1,NA,NA),lty=c(0,0,1,2),col=c("black","red",gray(.9),gray(.5)),ncol=1,bty="n",cex=.75)			
		}
	}
	# Baseline hazard, cumulated baseline hazard and survival functions
	if(any(which.plot[2:4])){
		cov_ThetaTheta = x$covar[[se]][-c(1:x$dim$p),-c(1:x$dim$p)]	
		xlim = c(x$knots$Alpha[1],quantile(x$data$time,prob=upper.quantile))
		plot_bh = function(j,V_x_X,Theta,covar,control,knots,pos,prob,xlim,...){
			M_Ppsi_Xm   = basis_mpl(V_x_X,knots,control$basis,control$order,which=as.numeric(j>1)+1)
			V_sd2.Hh0.X = diag(M_Ppsi_Xm%*%covar%*%t(M_Ppsi_Xm))
			pos.var     = V_sd2.Hh0.X>0
			V_Hh0_X     = c(M_Ppsi_Xm%*%matrix(Theta,ncol=1))[pos.var]
			V_x_X       = V_x_X[pos.var]
			V_sd.Hh0.X  = sqrt(V_sd2.Hh0.X[pos.var])
			upper       = V_Hh0_X+2*V_sd.Hh0.X
			lower       = V_Hh0_X-2*V_sd.Hh0.X
			lower[lower<0] = 0
			if(j==3){
				V_Hh0_X = exp(-V_Hh0_X)
				upper   = exp(-upper)
				lower   = exp(-lower)
			}
			plot(1,1,pch="",xlim=xlim,ylim=c(ifelse(j<3,0,min(upper[V_x_X<xlim[2]])),max(upper[V_x_X<xlim[2]])),axes=FALSE,
			     main=paste("Estimate of the",c(" baseline hazard"," cumulative baseline hazard"," baseline survival")[j]," function",sep=""),
			     xlab="Survival time",ylab=c(expression(h[0]*(t)),expression(H[0]*(t)),expression(S[0]*(t)))[j],...)
			rect(quant,0,max(V_x_X),max(upper)*1.5,col=gray(.9),border = NA)
			abline(v=quant,col=gray(.5),lty=2)			
			axis(2,las=2,...)
			axis(1,pos=ifelse(j<3,0,min(upper[V_x_X<xlim[2]])),...)
			x = c(V_x_X,V_x_X[length(V_x_X):1])
			y = c(upper,lower[length(V_x_X):1])
			confcol = paste(substr(colw[length(Theta)],1,7),70,sep="")
			polygon(x,y,col=confcol,border = "gray")
			lines(V_x_X,V_Hh0_X,lwd=1.1,col=colw[1])	
			legend(ifelse(j<3,"topleft","topright"),legend=c("Estimate","95% conf. interval",
			       as.expression(bquote(widehat(F)^-1*(.(prob))))),pch=c(NA,15,NA),lty=c(1,0,3),
			       col=c(colw[1],confcol,gray(.5)),ncol=1,bty="n",cex=.75)		
		}
		if(which.plot[2]){plot_bh(1,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim,...)}
		if(which.plot[3]){plot_bh(2,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim,...)}	
		if(which.plot[4]){plot_bh(3,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim,...)}				
	}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print.summary.coxph_mpl=function(x,se="M2QM2",...) {
	inf = x$inf
	cat("\n")
	print(inf$call)
	cat("\n-----\n\n")
	cat("Cox Proportional Hazards Model Fit Using MPL","\n\n\n")
	cat("Convergence               : ",
		if(inf$control$max.iter[1]==1){
		ifelse(inf$iter[3]<inf$control$max.iter[2],
		paste(c("Yes (",inf$iter[2]," iter.)"),collapse=""),"NO")		
		}else{
		ifelse(inf$iter[1]<inf$control$max.iter[1]&inf$iter[3]<inf$control$max.iter[2],
		paste(c("Yes (",inf$iter[1]," + ",inf$iter[2]," iter.)"),collapse=""),"NO")
		},"\n\n")
	cat("Data             : ",inf$data,"\n",sep="")	
	number=format(c(inf$dim$n,inf$dim$n.obs,inf$dim$n-inf$dim$n.obs))
	percentage=format(c(0,inf$dim$n.obs,inf$dim$n-inf$dim$n.obs)/inf$dim$n*100)
	cat("Number of obs.   : ",number[1],"\n",sep="")
	cat("Number of events : ",number[2]," (",percentage[2],"%)\n",sep="")
	cat("Number of cens.  : ",number[3]," (",percentage[3],"%)\n\n",sep="")
	cat("Regression parameters : ",deparse(inf$call[[2]]),"\n",sep="")
	printCoefmat(x$Beta, P.values=TRUE, has.Pvalue=TRUE,...)	
	cat("\nBaseline hasard parameters approximated using",
	    if(inf$control$basis=="uniform"){"a step function"}else{
	    if(inf$control$basis=="gaussian"){"Gaussian splines"}else{
	    if(inf$control$basis=="msplines"){"M-splines"}else{
	    "Epanechikov splines"}}},":\n")
	if(inf$control$basis=="uniform"){
	cat(paste(" (",inf$dim$m," equal events bins)\n",sep=""))
	}else{
	cat(paste(" (2 (min/max) + ",inf$control$n.knots[1]," quantile knots + ",inf$control$n.knots[2]," equally spaced knots",
	    if(inf$control$basis=="msplines"|inf$control$basis=="epanechikov"){
	paste(" +",inf$control$order,"(order) - 2")}," = ",inf$dim$m," parameters)",sep=""),"\n")}
	if(inf$full){
	printCoefmat(x$Theta, P.values=TRUE, has.Pvalue=TRUE,...)	
	}else{print(x$Theta,...)}
	cat("\n-----\n\n")
	}

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
coef.summary.coxph_mpl=function(object, parameters = "Beta", ...) {
	object[[parameters]]
	}
	
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coef.coxph_mpl=function(object, parameters = "Beta", ...) {
	out = object$coef[[parameters]]
	if(parameters == "Beta") names(out) = colnames(object$data$X)
	else names(out) = 1:object$dim$m
	out
	}
	
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
residuals.coxph_mpl=function(object,klein=FALSE,...) {
	control  = object$control
	M_Psi_Xm = basis_mpl(object$data$time,object$knots,control$basis,control$order,which=2)
	out = as.data.frame(matrix(NA,object$dim$n,4),...)
	colnames(out)  = c("time","observed","coxsnell","martingale")
	out$time       = object$data$time
	out$observed   = object$data$observed
	out$coxsnell   = exp(object$data$X%*%object$coef$Beta)*M_Psi_Xm%*%object$coef$Theta
	out$martingale = out$observed-out$coxsnell	
	if(klein){
		control$smooth=NULL
		fitcoxsnellres = try(coxph_mpl(Surv(out$coxsnell,object$data$observed)~1,control=control),silent=T)
		if(class(fitcoxsnellres)=="try-error"){klein=F
		}else{out$H0coxsnell = basis_mpl(out$coxsnell,fitcoxsnellres$knots,fitcoxsnellres$control$basis,
					fitcoxsnellres$control$order,which=2)%*%matrix(fitcoxsnellres$coef$Theta,ncol=1)
		}
	}	
	class(out) =c("residuals.coxph_mpl","data.frame")
	out
	}
	
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.residuals.coxph_mpl=function(x,ask=TRUE,which=1:3,upper.quantile=.95,...){
	prob = upper.quantile	
	which.plot=rep(TRUE,3)
	if(!is.null(which)){which.plot[-which]=FALSE}
	if(!any(colnames(x)=="H0coxsnell")){which.plot[3]=FALSE}
	if(sum(which.plot)==1){ask=FALSE}	
	if(ask){oask <- devAskNewPage(TRUE)
		on.exit(devAskNewPage(oask))
		}
	# martingale
	if(which.plot[1]){
		plot(1:nrow(x),x$martingale,col=(!x$observed)+1,ylab="",
		     xlab="Index",main="Martingale Residuals",axes=FALSE,
		     ylim=c(min(x$martingale),1+(1-min(x$martingale))*.075),...)
		mtext(expression(delta[i]-plain(e)^(x[i]^T*hat(beta))*widehat(H)[0](t[i])),2,padj=-2,...)
		abline(h=1,col="light gray")
		abline(h=0,col="light gray",lty=2)
		axis(1,...)
		axis(2,...)
		legend("top",ncol=2,legend=c("observed","censored"),col=c(1,2),pch=1,cex=.75)
		}
	# coxsnell
	if(which.plot[2]){
		plot(1:nrow(x),x$coxsnell,col=(!x$observed)+1,ylab="",
		     xlab="Index",main="Cox & Snell Residuals",axes=FALSE,
		     ylim=c(0,max(x$coxsnell)*1.075),...)
		mtext(expression(plain(e)^(x[i]^T*hat(beta))*widehat(H)[0](t[i])),2,padj=-2,...)
		abline(h=1,col="light gray",lty=2)
		axis(1,pos=0,...)
		axis(2,...)
		legend("top",ncol=2,legend=c("observed","censored"),col=c(1,2),pch=1,cex=.75)
		}
	# plot of Klein and Moeschberger
	if(which.plot[3]){
		plot(x$coxsnell[order(x$coxsnell)],x$H0coxsnell[order(x$H0coxsnell)],
		     ylab="",type='s',main="Cox & Snell Residuals",axes=FALSE,
		     xlab=expression(r[Ci] ==plain(e)^(x[i]^T*hat(beta))*widehat(H)[0](t[i])),
		     ylim=c(0,max(x$coxsnell)*1.075),...)
		mtext(expression(hat(H)[0](r[Ci])),2,padj=-2)
		abline(0,1,col="red",lty=2)
		abline(v=quantile(x$coxsnell,0.95),col="light gray",lty=2)
		abline(h=quantile(x$H0coxsnell,0.95),col="light gray",lty=2)		
		axis(1,pos=0,...)
		axis(2,pos=0,...)
		legend("top",legend=as.expression(bquote(widehat(F)^-1*(.(prob)))),col="light gray",
		       pch=NA,cex=.75,lty=2,bg="white")
		}
	}

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
predict.coxph_mpl=function(object,se="M2QM2",type="risk",i=NULL,time=NULL,upper.quantile=.95,...) {
	prob  = upper.quantile
	covar = object$covar[[se]]
	Beta  = object$coef$Beta
	Theta = object$coef$Theta	
	p     = object$dim$p
	m     = object$dim$m	
	if(length(i)>1){warning("only the first observation will be considered\n",call. = FALSE)}
	# time
	if(is.null(time)){
		n.x   = 1000
		V_x_X = seq(object$knots$Alpha[1],max(object$knots$Alpha),length=n.x)
	}else{	n.x   = length(time) 
		V_x_X = time
	}
	# x
	xTB = if(is.null(i)){apply(object$data$X,2,mean)}else{object$data$X[i[1],,drop=FALSE]}
	Mu  = c(exp(xTB%*%Beta))
	# risk
	out = data.frame(time = V_x_X, mid = NA, se = NA, low = NA, high = NA,...)
	if(type=="risk"){
		M_psi_Xm = basis_mpl(V_x_X,object$knots,object$control$basis,object$control$order,which=1)
		out$mid  = Mu*M_psi_Xm%*%Theta
		# correction factor
		M_corr_mpm = matrix(c(rep(M_psi_Xm%*%Theta*Mu,p)*rep(Beta,each=n.x),Mu*M_psi_Xm),ncol=p+m)
		out$se     = sqrt(diag(M_corr_mpm%*%object$covar[[se]]%*%t(M_corr_mpm)))
		out$low    = out$mid - 2*out$se; out$low[out$low<0] = 0
		out$high   = out$mid + 2*out$se
	# survival
	}else{
		M_Psi_Xm = basis_mpl(V_x_X,object$knots,object$control$basis,object$control$order,which=2)
		out$mid  = exp(-Mu*M_Psi_Xm%*%Theta)
		# correction factor
		M_corr_mpm = rep(-out$mid*Mu,p+m)*matrix(c(rep(M_Psi_Xm%*%Theta,p)*rep(Beta,each=n.x),M_Psi_Xm),ncol=p+m)
		out$se     = sqrt(diag(M_corr_mpm%*%object$covar[[se]]%*%t(M_corr_mpm)))
		out$low    = out$mid - 2*out$se; out$low[out$low<0] = 0
		out$high   = out$mid + 2*out$se; out$high[out$high>1] = 1	
	}
	# out
	attributes(out)$inf = list(i=i[1], upper.quantile=prob, upper.value = quantile(object$data$time,prob), 
		max = max(object$data$time), user.time = !is.null(time), m=m, risk = type=="risk")
	colnames(out)[2]=type
	class(out) =c("predict.coxph_mpl","data.frame")
	out
	}
	
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.predict.coxph_mpl=function(x,...){
	inf     = attr(x,"inf")
	colw    = terrain.colors(3)[1:2]
	pos.var = x$se>0
	main    = paste(ifelse(inf$risk,"Predicted instantaneous risk at time t",
	          "Predicted probability of survival after time t"),"\n",
	          ifelse(is.null(inf$i),"for 'average' covariates",paste("for observation",inf$i)),sep="")
	plot(1,1,pch="",axes=FALSE,xlab="Time",ylim=c(0,max(x$high)),main=main,
	     xlim=if(inf$user.time){c(0.5,nrow(x)+.5)}else{c(min(x$time),inf$upper.value)},ylab="",...)
	axis(2,...)
	abline(h=0,col="light gray")
	mtext(if(inf$risk){expression(hat(h)(t[i]))}else{expression(widehat(S)(t[i]))},2,padj=-2,...)
	if(inf$user.time){
		axis(1,at=1:nrow(x),labels=x$time,tick=TRUE,pos=0,...)
		for(tw in 1:nrow(x)){
			if(x$se[tw]>0){
				arrows(tw, y0=x$low[tw], y1=x$high[tw], angle=90, code=3, col=colw[2],...)
			}
		}
		points(1:nrow(x),x[,2],col=colw[1],...)
		legend("topleft",legend=c("Estimate","95% conf. interval"),pch=c(NA,15),lty=c(1,0),
		       col=c(colw[1],confcol),ncol=1,bty="n",cex=.75)				
	}else{
		axis(1,pos=0,...)
		rect(inf$upper.value,0,inf$max,max(x$high)*1.5,col=gray(.9),border = NA)
		if(inf$upper.quantile<1)abline(v=inf$upper.value,col=gray(.5),lty=2)
		xx = c(x$time[pos.var],x$time[pos.var][length(x$time[pos.var]):1])
		yy = c(x$high[pos.var],x$low[pos.var][length(x$time[pos.var]):1])
		confcol = paste(substr(colw[2],1,7),70,sep="")
		polygon(xx,yy,col=confcol,border = "gray",...)
		lines(x$time,x[,2],lwd=1.1,col=colw[1],...)	
		legend(ifelse(inf$risk,"topleft","topright"),legend=c("Estimate","95% conf. interval",
		       as.expression(bquote(widehat(F)^-1*(.(inf$upper.quantile))))),pch=c(NA,15,NA),lty=c(1,0,3),
		       col=c(colw[1],confcol,gray(.5)),ncol=1,bty="n",cex=.75)		
	}
	}


	
