
#On va essayer de generaliser notre fonction
#look at _rmcoll function from stata to remove multicollinear variables

#TODO:
#ajouter la possibilite de partir avec les coefs que l'on souhaite

#Mettre du multicore pour le calcul des jacobiennes
feNmlm <- function(fml,data,linear.fml,start,lower,upper,env,dummy,start.init,nl.gradient,linear.start=0,jacobian.method=c("simple","Richardson"),useHessian=TRUE,d.hessian,opt_method=c("nlminb","optim"),debug=FALSE,family=c("poisson","negbin","logit"),opt.control=list(),optim.method="BFGS",...){

	opt_method <- match.arg(opt_method)
	jacobian.method <- match.arg(jacobian.method)
	family = match.arg(family)

	# Some settings (too complicated to be tweaked by the user)
	# Nber of evaluations of the NL part to be kept in memory
	# Default keeps the two last evaluations
	NLsave=2
	# Not implemented yet (maybe never, cause it's not the point):
	rmcoll=FALSE

	famFuns = switch(family,
						  poisson = ml_poisson(),
						  negbin = ml_negbin(),
						  logit = ml_logit(),
						  tobit = ml_tobit(),
						  probit = ml_probit(),
						  gaussian = ml_gaussian())

	stopifnot(class(fml)=="formula")
	if(length(fml)!=3) stop("The formula must be two sided.\nEG: a~exp(b/x), or a~0 if there is no nonlinear part.")
	call = match.call()
	dataNames = names(data)
	#The LHS must contain only values in the DF
	namesLHS = all.vars(fml[[2]])
	if(!all(namesLHS%in%dataNames)) stop("Some elements on the LHS of the formula are not in the dataset:\n",paste0(namesLHS[!namesLHS%in%dataNames],collapse=", "))
	#Now the nonlinear part:
	allnames = all.vars(fml[[3]])
	nonlinear.params = allnames[!allnames %in% dataNames]
	nonlinear.varnames = allnames[allnames %in% dataNames]

	#The dependent variable: lhs==left_hand_side
	lhs = as.vector(eval(fml[[2]],data))

	#creation de l'environnement
	if(missing(env)) env <- new.env()
	else stopifnot(class(env)=="environment")

	#
	# First check
	#

	#NA are not allowed !!!
	if(anyNA(lhs)) stop("The left hand side of the fomula has NA values. Please provide data without NA.")
	if(family%in%c("poisson","negbin") & any(lhs<0)) stop("Negative values of the dependant variable \nare not allowed for the \"",family,"\" family.",sep="")
	if(family%in%c("logit") & !all(lhs==0 | lhs==1)) stop("The dependant variable has values different from 0 or 1.\nThis is not allowed with the \"logit\" family.")

	#
	# Controls and setting of the linear part:
	#

	isLinear = FALSE
	if(!missing(linear.fml)){
		isLinear <- TRUE
		if(class(linear.fml)!="formula" | length(linear.fml)!=2) stop("'linear.fml' must be a formula like, for ex., ~x1+x2-1")
		linear.varnames <- all.vars(linear.fml)
		if(!all(linear.varnames%in%dataNames)) stop(paste("In 'linear.fml', some variables are not in the data:\n",paste(linear.varnames[!linear.varnames%in%dataNames],collapse=', '),".",sep=""))
		if(!missing(dummy)){
			#if dummies are provided, we make sure there is an
			#intercept so that factors can be handled properly
			linear.fml = update(linear.fml,~.+1)
		}
		linear.mat = model.matrix(linear.fml,data)
		linear.params <- colnames(linear.mat)
		N_linear <- length(linear.params)
		if(anyNA(linear.mat)) stop("Evaluation of the linear part returns NA.\nAs NAs are not supported, please solve this problem first before running the algorithm.")
		if(!is.numeric(linear.start)) stop("'linear.start' must be numeric!")
	} 	else linear.params <- linear.start <- linear.varnames <- NULL
	params <- c(nonlinear.params,linear.params)
	lparams <- length(params)
	varnames <- c(nonlinear.varnames,linear.varnames)
	#attention les parametres non lineaires peuvent etre vides
	if(length(nonlinear.params)==0) isNL = FALSE
	else isNL = TRUE

	# Control for NAs
	if(anyNA(data[,varnames])){
		varWithNA = varnames[which(apply(data[,varnames],2,anyNA))]
		stop("Some variables used for this estimation contain NA. NA are not supported, please remove them first.\nFYI, the variables are:\n",paste0(varWithNA,collapse = "; "),call. = FALSE)
	}

	#
	# Handling dummies
	#

	isDummy = FALSE
	if(!missing(dummy)){
		# TODO:
		# - mettre un meilleur controle des classes pures (mettre un while, ici ne suffit pas)
		# - ajouter la possibilite de mettre une "reference" (ie oter une dummy)

		isDummy = TRUE
		if(class(dummy)!="character" | any(!dummy%in%names(data))) stop("Dummies must be a variable name!")
 		#qui = which(sapply(dummy,function(x) is.factor(class(data[[x]]))))
 		#if(length(qui)>0) stop("The variable(s) defining the cluster(s) must be a factor!\nIt concerns:",paste(dummy[qui],collapse=", "))

		g = length(dummy)
		dum_all = dum_names = S_all = list()
		obs2remove = c()
		for(i in 1:g){
			dum = data[[dummy[i]]]
			#in order to avoid "unclassed" values > real nber of classes: we re-factor the dummy
			dum = factor(dum)
			dum_names[[i]] = thisNames = levels(dum)
			dum = unclass(dum)
			#We create the matrix to compute the sums per class:

			n = length(dum)
			k = max(dum)
			S = Matrix::Matrix(0,k,n,sparse = TRUE)
			S[cbind(dum,1:n)] = 1

			dum_all[[i]] = dum
			S_all[[i]] = S

			#We delete "all zero" outcome
			sum_y_clust = as.vector(S%*%lhs)
			n_perClust = as.vector(rowSums(S))
			if(family%in%c("poisson","negbin")) qui = which(sum_y_clust==0)
			if(family=="logit") qui = which(sum_y_clust==0 | sum_y_clust==n_perClust)

			if(length(qui>0)){
				#We first delete the data:
				dummyOmises = thisNames[qui] #not currently used
				obs2remove = unique(c(obs2remove,which(dum%in%qui)))
			}
		}

		#We remove the problems
		if(length(obs2remove)>0){
			data = data[-obs2remove,]

			#We recreate the linear matrix and the LHS
			if(isLinear) linear.mat = model.matrix(linear.fml,data)
			lhs = eval(fml[[2]],data)

			#Then we recreate the dummies
			for(i in 1:g){
				dum = data[[dummy[i]]]
				dum = factor(dum)
				dum_names[[i]] = levels(dum)
				dum_all[[i]] = dum = unclass(dum)
				#We create the matrix to computes the sums per class:
				n = length(dum)
				k = max(dum)
				S = Matrix::Matrix(0,k,n,sparse = TRUE)
				S[cbind(dum,1:n)] = 1
				S_all[[i]] = S
			}

			#Then the warning message
			#TODO : make it fancier
			warning(length(dummyOmises)," clusters (",length(obs2remove)," observations) removed because of only ",ifelse(family=="logit","zero/one","zero")," outcomes.",call. = FALSE,immediate. = TRUE)
		}

		#If there is a linear intercept, we withdraw it
		#we add some more controls to avoid colinear variables

		#TAKES too long => not run // user should be aware of collinearity issues
		if(FALSE & isLinear){
			tol = 1e-8
			var2remove = c()
			for(i in 1:g){
				S = S_all[[i]]
				k = nrow(S)
				ptm_est = proc.time()
				middle = Matrix::Diagonal(k,1/rowSums(S))
				#Orthogonal projection on the dummies subspace
				# estimation = linear.mat-t(S)%*%middle%*%S%*%linear.mat
				estimation = linear.mat-crossprod(sqrt(middle)%*%S)%*%linear.mat
				var2remove = unique(c(var2remove,which(colSums(estimation**2)<tol)))
				cat("Estimation collinearity:",gt(ptm_est),"\n")
				print(var2remove)
			}
			if(length(var2remove)>0){
				name_withdrawn = colnames(linear.mat)[var2remove]
				if(ncol(linear.mat)==length(var2remove)) isLinear = FALSE
				else{
					linear.mat = linear.mat[,-var2remove,drop=FALSE]
					linear.params <- colnames(linear.mat)
					N_linear <- length(linear.params)
					params <- c(nonlinear.params,linear.params)
					lparams <- length(params)
					varnames <- c(nonlinear.varnames,linear.varnames)
				}
				#we don't bother to tell we've withdrawn the intercept
				name_withdrawn = name_withdrawn[name_withdrawn!="(Intercept)"]
				if(length(name_withdrawn)>=1) warning("The following variables were removed because of collinearity with the clusters:\n",paste0(name_withdrawn,collapse=" ; "),immediate. = TRUE,call. = FALSE)
			}
		}
		#We drop the intercept:
		if("(Intercept)"%in%colnames(linear.mat)){
			var2remove = which(colnames(linear.mat)=="(Intercept)")
			if(ncol(linear.mat)==length(var2remove)) isLinear = FALSE
			else{
				linear.mat = linear.mat[,-var2remove,drop=FALSE]
				linear.params <- colnames(linear.mat)
				N_linear <- length(linear.params)
				params <- c(nonlinear.params,linear.params)
				lparams <- length(params)
				varnames <- c(nonlinear.varnames,linear.varnames)
			}
		}
	} else if(FALSE & "(Intercept)"%in%linear.params){
		#We treat the intercept as a dummy:
		isDummy = TRUE
		n = nrow(data)
		dum = rep(1,n)
		dum_names = "(Intercept)"
		sum_y_clust = sum(lhs)
		n_perClust = n
		S = matrix(1,nrow=1,ncol = n)
		#we withdraw the intercept
		if(ncol(linear.mat)==1){
			isLinear = FALSE
			linear.params = NULL
			N_linear = 0
			params <- c(nonlinear.params)
			lparams <- length(params)
			varnames <- c(nonlinear.varnames)
		}
		else{
			qui = which(linear.params=="(Intercept)")
			linear.mat = linear.mat[,-qui,drop=FALSE]
			linear.params <- colnames(linear.mat)
			N_linear <- length(linear.params)
			params <- c(nonlinear.params,linear.params)
			lparams <- length(params)
			varnames <- c(nonlinear.varnames,linear.varnames)
		}
		#we set it as a list for compatibility
		dum_all = list(dum)
		S_all = list(S)
	}

	#
	# Collinearity
	#

	if(rmcoll & isLinear){
		#We withdraw perfectly collinear variables
		#no: too long. User should be aware beforehand
		l = lm(update(fml,linear.fml),data)
		m = coef(l)
		if(anyNA(m)){
			name_withdrawn = names(m)[is.na(m)]
			linear.mat = linear.mat[,!colnames(linear.mat)%in%name_withdrawn,drop=FALSE]
			linear.params <- colnames(linear.mat)
			N_linear <- length(linear.params)
			params <- c(nonlinear.params,linear.params)
			lparams <- length(params)
			varnames <- c(nonlinear.varnames,linear.varnames)
			stop("The following variables were removed because of perfect collinearity:\n",paste0(name_withdrawn,collapse=" ; "),immediate. = TRUE,call. = FALSE)
		}
		#To add: the possibility to delete collinear variables (when hessian cant be computed)
	}


	#
	# Checks for MONKEY TEST
	#

	if(lparams==0) stop("No parameter to be estimated.")
	if(!is.logical(useHessian)) stop("'useHessian' must be of type 'logical'!")

	#
	# Controls: The non linear part
	#

	if(isNL){
		if(missing(start.init)){
			if(missing(start)) stop("There must be starting values.")
			if(typeof(start)!="list") stop("start must be a list.")
			if(any(!names(start) %in% params)) stop(paste("Some parameters in 'start' are not in the formula:\n",paste(names(start)[!names(start) %in% params],collapse=", "),".",sep=""))
			if(any(!nonlinear.params %in% names(start))) stop(paste("Hey, some parameters have no starting values:\n",paste(nonlinear.params[!nonlinear.params%in%names(start)],collapse=", "),".",sep=""))
		}
		else{
			if(length(start.init)>1) stop("start.init musn't be a vector.")
			if(class(start.init)!="numeric") stop("start.init must be numeric!")
			if(!is.finite(start.init)) stop("Infinites values as starting values, you must be kidding me...")

			if(missing(start)){
				start <- list()
				start[nonlinear.params] <- start.init
			}
			else{
				if(typeof(start)!="list") stop("start must be a list.")
				if(any(!names(start) %in% params)) stop(paste("Some parameters in 'start' are not in the formula:\n",paste(names(start)[!names(start) %in% params],collapse=", "),".",sep=""))

				missing.params <- nonlinear.params[!nonlinear.params%in%names(start)]
				start[missing.params] <- start.init
			}
		}
	} else start <- list()

	#
	# Controls: The upper and lower limits
	#

	if(!missing(lower)){
		if(typeof(lower)!="list") stop("'lower' MUST be a list.")
		if(any(!names(lower)%in%params)){
			text <- paste("Hey, some parameters in 'lower' are not in the formula:\n",paste(names(lower)[!names(lower)%in%params],collapse=", "),".",sep="")
			stop(text)
		}
	}
	if(!missing(upper)){
		if(typeof(upper)!="list") stop("'upper' MUST be a list.")
		if(any(!names(upper)%in%params)){
			text <- paste("Hey, some parameters in 'upper' are not in the formula:\n",paste(names(upper)[!names(upper)%in%params],collapse=", "),".",sep="")
			stop(text)
		}
	}

	# Now setting upper and lower

	if(!missing(lower)){
		lower[params[!params%in%names(lower)]] <- -Inf
		lower <- unlist(lower[params])
	}	else {
		#TODO B
		lower <- rep(-Inf,lparams)
		names(lower) <- params
	}
	if(!missing(upper)){
		upper[params[!params%in%names(upper)]] <- Inf
		upper <- unlist(upper[params])
	}	else upper <- rep(Inf,lparams)

	lower <- c(lower)
	upper <- c(upper)

	#
	# Controls: user defined gradient
	#

	if(!missing(nl.gradient)){
		isGradient=TRUE
		if(class(nl.gradient)!="formula" | length(nl.gradient)==3) stop("'nl.gradient' must be a formula like, for ex., ~f0(a1,x1,a2,x2). f0 giving the gradient.")
	}
	else isGradient=FALSE

	if(!missing(d.hessian)){
		hessianArgs=list(d=d.hessian)
	} else hessianArgs = NULL
	assign("hessianArgs",hessianArgs,env)

	#
	# Sending data in the environment
	#

	#Initial checks are done
	nonlinear.params <- names(start) #=> in the order the user wants

	# control for the linear start => we can provide coefficients
	# from past estimations. Coefficients that are not provided are set
	# to 0
	if(length(linear.start)>1){
		what = linear.start[linear.params]
		what[is.na(what)] = 0
		linear.start = what
	}
	start[linear.params] <- linear.start
	params <- names(start)
	start <- unlist(start)
	start <- c(start)
	lparams <- length(params)
	names(start) <- params

	#For the negative binomial:
	if(family=="negbin"){
		params = c(params,".theta")
		start = c(start,1)
		names(start) = params
		lower = c(lower,1e-4)
		# 		theta = 1 # no matter the init, its handled in the function getting theta
		# 		assign(".theta",theta,env)
	} else if(family=="tobit"){
		params = c(params,".sigma")
		start = c(start,1)
		names(start) = params
		lower = c(lower,1e-4)
	}

	#On balance les donnees a utiliser dans un nouvel environnement
	for(i in varnames) assign(i, data[[i]],env)
	if(isLinear) assign("linear.mat",linear.mat,env)
	if(isGradient) assign(".call_gradient",nl.gradient[[2]],env)
	#The dummies
	assign("isDummy",isDummy,env)
	if(isDummy){
		assign(".dummy",dum_all,env)
		assign(".S",S_all,env)
		#the saved dummies
		assign(".savedDummy",rep(0,length(lhs)),env)
	}
	#other
	assign(".nl.call",fml[[3]],env)
	assign(".lhs",lhs,env)
	assign("isNL",isNL,env)
	assign("isLinear",isLinear,env)
	assign("isGradient",isGradient,env)
	assign("linear.params",linear.params,env)
	assign("nonlinear.params",nonlinear.params,env)
	assign("params",params,env)
	assign("nobs",length(lhs),env)
	assign("debug",debug,env)
	assign("jacobian.method",jacobian.method,env)
	assign(".famFuns",famFuns,env)
	assign(".family",family,env)
	assign("iter",0,env)
	#Pour gerer les valeurs de mu:
	assign(".coefMu",list(),env)
	assign(".valueMu",list(),env)
	assign(".wasUsed",TRUE,env)
	#Pour les valeurs de la Jacobienne non lineaire
	assign(".JC_nbSave",0,env)
	assign(".JC_nbMaxSave",1,env)
	assign(".JC_savedCoef",list(),env)
	assign(".JC_savedValue",list(),env)

	#On teste les valeurs initiales pour informer l'utilisateur
	for(var in nonlinear.params) assign(var,start[var],env)
	mu = eval(fml[[3]],envir= env)

	#On sauvegarde les valeurs de la partie non lineaire
	assign(".nbMaxSave",NLsave,env) #nombre maximal de valeurs a sauvegarder
	assign(".nbSave",1,env) #nombre de valeurs totales sauvegardees
	assign(".savedCoef",list(start[nonlinear.params]),env)
	assign(".savedValue",list(mu),env)
	if(isLinear) mu <- mu + c(linear.mat%*%unlist(start[linear.params]))

	if(length(mu)!=nrow(data)) stop("Wow, must be a big problem... length(lhs)!=length(eval(fml))")
	if(anyNA(mu)) stop("Hum, must be a problem, evaluating formula returns NA.\nMaybe only these starting values can't be computed, or maybe there's another BIGGER problem.")

	# ^_- the theta was here

	#Check of the user-defined gradient, if given
	if(isGradient){
		for(nom in nonlinear.params) assign(nom,start[nom],env)
		test <- eval(nl.gradient[[2]],envir=env)
		if(!class(test)%in%c("list","data.frame")) stop("The function called by 'nl.gradient' must return an object of type 'list' or 'data.frame'.")
		if(!all(nonlinear.params%in%names(test))) stop(paste("The gradient must return a value for each parameter. Some are missing:\n", paste(nonlinear.params[!nonlinear.params%in%names(test)],collapse=", "),".",sep=""))
		if(!all(names(test)%in%nonlinear.params)) warning(paste("Some values given by 'nl.gradient' are not in the parameters:\n", paste(names(test)[!names(test)%in%nonlinear.params],collapse=", "),".",sep=""))
		if(mean(sapply(test[nonlinear.params],length))!=length(lhs)) stop("Strange, the length of the vector returned by 'nl.gradient' does not match with the data.")
		#we save 1 gradient:
		jacob.mat = as.matrix(test[nonlinear.params])
		assign(".JC_nbSave",1,env)
		assign(".JC_savedCoef",list(start[nonlinear.params]),env)
		assign(".JC_savedValue",list(jacob.mat),env)
	}

	#Mise en place du calcul du gradient
	gradient = NULL
	if(opt_method=="nlminb" | optim.method%in%c("BFGS","CG","L-BFGS-B")) gradient = ll_glm_gradient
	hessian <- NULL
	if(useHessian) hessian <- ll_glm_hessian

	L = list(...)
	if(!is.null(L$give.params) && L$give.params) return(list(coef=start,env=env))

	#
	# Maximizing the likelihood
	#

	opt <- NULL
	if(opt_method=="nlminb"){
# 		try(opt <- nlminb(start=start,objective=ll_glm,env=env,lower=lower,upper=upper,gradient=gradient,hessian=hessian),silent=FALSE)
		opt <- stats::nlminb(start=start,objective=ll_glm,env=env,lower=lower,upper=upper,gradient=gradient,hessian=hessian,control=opt.control)
	} else if(opt_method=="optim"){
		try(opt <- stats::optim(par=start,fn=ll_glm,gr=gradient,env=env,hessian=FALSE,control=opt.control),silent=FALSE)
		opt$objective <- opt$value
		opt$message = opt$convergence
	}

	if(is.null(opt)){
		stop("Could not achieve maximization.")
	}

	convStatus = TRUE
	if(opt_method=="nlminb" && !opt$message %in% c("relative convergence (4)","both X-convergence and relative convergence (5)")){
		warning("The result is not reliable, the optimization did not converge.",call. = FALSE)
		convStatus = FALSE
	} else if(opt_method=="optim"){
		if(opt$convergence!=0){
			warning("The result is not reliable, the optimization did not converge.\nConvergence code",opt$convergence,". See optim help page for more info.",call. = FALSE)
			opt$message = paste0("No convergence. Status ",opt$convergence)
			convStatus = FALSE
		} else opt$message = "Convergence"
	}

	#
	# Computing all the parameters
	#

	coef <- opt$par

	#The Hessian
	hessian = ll_glm_hessian(coef,env=env)

	var <- NULL
	try(var <- solve(hessian),silent = TRUE)
	if(is.null(var)){
		#var <- MASS::ginv(hessian)
		if(!isDummy){
			Qr = qr(hessian)
			collvar = params[Qr$pivot[-(1:Qr$rank)]]
			warning("Could not achieve to get the covariance matrix. The information matrix was singular.\nTry to manually suppress some collinear variables.FYI the suspects are:\n",paste(collvar,collapse=", "),call. = FALSE)
		} else {
			var <- MASS::ginv(hessian)
			collvar = params[diag(var)==0]
			if(length(collvar)>0){
				warning("Could not achieve to get the covariance matrix. The information matrix was singular.\nTry to manually suppress some collinear variables.FYI the suspects (collinear with the dummies) are:\n",paste(collvar,collapse=", "),call. = FALSE)
			} else {
				warning("Could not achieve to get the covariance matrix. The information matrix was singular.\nTry to manually suppress some collinear variables. They may be collinear with the dummies.",call. = FALSE)
			}
		}

		var = sd = NA
	} else sd <- sqrt(diag(var)) #add custom warnings when negative value of the diag
	zvalue <- coef/sd
	pvalue <- 2*pnorm(-abs(zvalue))

	coeftable <- cbind(coef,sd,zvalue,pvalue)
	colnames(coeftable) <- c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
	rownames(coeftable) <- params
	class(coeftable) <- "coeftest"

	mu = get_mu(coef,env)

	#calcul pseudo r2
	loglik <- -opt$objective #moins car la fonction minimise
	model0 <- ll0_nlglm(lhs,env)
	ll_null <- model0$loglik
	# degres de liberte
	df_k = length(coef)
	if(isDummy) df_k = df_k + sum(sapply(dum_all,max))
	pseudo_r2 <- 1 - (loglik-df_k)/ll_null

	#deprecated  resids
	#null.resids <- lhs-model0$constant
	#new null resids
	# null.resids <- lhs-mean(lhs)

	#Calcul residus
	expected.predictor = famFuns$expected.predictor(mu,env)
	resids = lhs - expected.predictor

	#calcul squared corr
	sq.cor = stats::cor(lhs,expected.predictor)**2
	#calcul r2 naif
	naive.r2 = 1-sum(resids**2)/sum((lhs-mean(lhs))**2)

	#The scores
	scores = ll_glm_scores(coef,env)

	if(missing(linear.fml)) linear.fml <- NULL

	res <- list(coef=coef,coeftable=coeftable,loglik=loglik,iterations=opt$iterations,n=length(lhs),k=df_k,call=call,nonlinear.fml=fml,linear.formula=linear.fml,ll_null=ll_null,pseudo_r2=pseudo_r2,naive.r2=naive.r2,message=opt$message,convStatus=convStatus,sq.cor=sq.cor,expected.predictor=expected.predictor,hessian=hessian,cov.unscaled=var,sd=sd,scores=scores,family=family)

	#Dummies
	if(isDummy){
		#NOT YET IMPLEMENTED FOR VARIOUS CLUSTERS
		mu_noDum = attr(mu,"mu_noDum")
		dummies = getDummies(mu_noDum,env,coef)
		if(FALSE & all(dum_names == "(Intercept)")){
			#get the variance of the intercept
			#and include it in the coefficients
			ptm = proc.time()
			res$Intercept = dummies
			jacob.mat = get_Jacobian(coef,env)
			sd = famFuns$dumVar(jacob.mat,mu,env,coef)
			zvalue = dummies/sd
			pvalue = 2*pnorm(-abs(zvalue))
			line_intercept = matrix(c(dummies,sd,zvalue,pvalue),1,4)
			rownames(line_intercept) = "(Intercept)"
			coeftable = rbind(line_intercept,res$coeftable)
			print(coeftable)
			cat("SD of the intercept in:",(proc.time()-ptm)[3],"\n")
		} else {
			#TODO:
			#	- get the dummies for each cluster with a reference
			res$dummies = dummies
			res$clusterNames = dummy

			id_dummies = list()
			for(i in 1:length(dummy)){
				id_dummies[[dummy[i]]] = factor(dum_all[[i]],labels=dum_names[[i]])
			}
			res$id_dummies = id_dummies
			clustSize = sapply(dum_all,max)
			names(clustSize) = dummy
			res$clusterSize = clustSize
			if(length(obs2remove)>0){
				res$obsRemoved = obs2remove
				res$clusterRemoved = dummyOmises
			}
		}
	}

	if(family == "negbin"){
		theta = coef[".theta"]
		res$theta = theta
	}

	class(res) <- "feNmlm"

	return(res)
}


femlm = function(linear.fml,data,dummy,linear.start=0,useHessian=TRUE,opt_method=c("nlminb","optim"),debug=FALSE,family=c("poisson","negbin","logit"),opt.control=list(),optim.method="BFGS",...){

	# some controls
	stopifnot(class(linear.fml)=="formula")

	# we create the call to nlglm
	fml = update(linear.fml,.~0)
	fml[[3]] = 0
	myLinearFml = linear.fml
	myLinearFml[[2]] = NULL

	# The call
	mf = match.call(expand.dots = TRUE) #; return(mf)

	# We "reshape" the call to feNmlm
	mf[[1]] = as.name("feNmlm")
	mf[["linear.fml"]] = myLinearFml
	mf[["fml"]] = fml

	# The proper call // to avoid problems of 'emboitement'
	eval(mf,sys.frame(-1))
}

ll_glm_hessian <- function(coef,env){
	# Computes the hessian
	# cat("in Hessian:",as.vector(coef),"\n")
	params <- get("params",env)
	names(coef) <- params
	nonlinear.params <- get("nonlinear.params",env)
	k <- length(nonlinear.params)
	isNL <- get("isNL",env)
	hessianArgs = get("hessianArgs",env)
	famFuns = get(".famFuns",env)
	family = get(".family",env)
	y = get(".lhs",env)
	isDummy = get("isDummy",env)
	mu = get_savedMu(coef,env)

	jacob.mat = get_Jacobian(coef,env)

	#hessVar = getHessianLinear(jacob.mat,y,mu,env,coef)

	ll_d2 = famFuns$ll_d2(y,mu,coef)
	if(isDummy){
		dxi_dbeta = deriv_xi(jacob.mat,ll_d2,env,coef)
		jacob.mat = jacob.mat + dxi_dbeta
	} else dxi_dbeta = 0

	hessVar = crossprod(jacob.mat, jacob.mat * ll_d2)

	if(isNL){
		#we get the 2nd derivatives
		z = numDeriv::genD(evalNLpart,coef[nonlinear.params],env=env,method.args = hessianArgs)$D[,-(1:k),drop=FALSE]
		ll_dl = famFuns$ll_dl(y=y,mu=mu,coef=coef,env=env)
		id_r = rep(1:k,1:k)
		id_c = c(sapply(1:k,function(x) 1:x),recursive=TRUE)
		H = matrix(0,nrow=k,ncol=k)
		H[cbind(id_r,id_c)] = H[cbind(id_r,id_c)] = colSums(z*ll_dl)
	} else H = 0

	#on ajoute la partie manquante
	if(isNL) hessVar[1:k,1:k] = hessVar[1:k,1:k] + H

	if(family=="negbin"){
		theta = coef[".theta"]
		ll_dx_dother = famFuns$ll_dx_dother(y,mu,coef,env)

		if(isDummy) dxi_dother = deriv_xi_other(ll_d2,ll_dx_dother,env,coef)
		else dxi_dother = 0

		#calcul des derivees secondes vav de theta
		h.theta.L = famFuns$hess.thetaL(theta,jacob.mat,y,dxi_dbeta,dxi_dother,ll_d2,ll_dx_dother)
		hessVar = cbind(hessVar,h.theta.L)
		h.theta = famFuns$hess_theta_part(theta,y,mu,dxi_dother,ll_dx_dother,ll_d2)
		hessVar = rbind(hessVar,c(h.theta.L,h.theta))

# 		theta = attr(mu,".theta")
# 		print(theta)
		#Jacob.mat est la derivee totale de mu vav de beta (ie avec les dummies)
# 		hessVar = hessVar + famFuns$hess_theta_part(theta,jacob.mat,y,mu)
	} else if(family=="tobit"){
		sigma = coef[".sigma"]
		h.sigma.L = famFuns$hess.sigmaL(sigma,jacob.mat,y,mu,dxi_dbeta,ll_d2)
		hessVar = cbind(hessVar,h.sigma.L)
		h.sigma = famFuns$hess.sigma(sigma,y,mu)
		hessVar = rbind(hessVar,c(h.sigma.L,h.sigma))
		print(hessVar)
	}
# print(hessVar) ; print(class(hessVar))
	- hessVar
}

ll_glm_gradient <- function(coef,env){
	# cat("gradient:\n") ; print(as.vector(coef))

	params = get("params",env)
	names(coef) = params
	nonlinear.params = get("nonlinear.params",env)
	linear.params = get("linear.params",env)
	famFuns = get(".famFuns",env)
	family = get(".family",env)
	y = get(".lhs",env)
	mu = get_savedMu(coef,env)

	#calcul de la jacobienne
	res <- list() #stocks the results

	# cat("\tgetting jacobian")
	# ptm = proc.time()
	jacob.mat = get_Jacobian(coef,env)
	# cat("in",(proc.time()-ptm)[3],"s.\n")
# 	isDummy = get("isDummy",env)
# 	if(isDummy){
# 		ll_d2 = famFuns$ll_d2(y,mu,coef)
# 		dxi_dbeta = deriv_xi(jacob.mat,ll_d2,env,coef)
# 		print(dim(dxi_dbeta))
# 		jacob.mat = jacob.mat + dxi_dbeta
# 	}

	# cat("\tComputing gradient ")
	# ptm = proc.time()
	# res = famFuns$grad(jacob.mat,y,mu,env,coef)
	res = getGradient(jacob.mat,y,mu,env,coef)
	# cat("in",(proc.time()-ptm)[3],"s.\n")
	names(res) = c(nonlinear.params,linear.params)

	if(family=="negbin"){
		theta = coef[".theta"]
		res[".theta"] = famFuns$grad.theta(theta,y,mu)
	}
# print(res)
	return(-unlist(res[params]))
}

ll_glm_scores <- function(coef,env){
	#Computes the scores (Jacobian)
	params = get("params",env)
	names(coef) <- params
	famFuns = get(".famFuns",env)
	family = get(".family",env)
	y = get(".lhs",env)
	mu = get_savedMu(coef,env)

	jacob.mat = get_Jacobian(coef,env)
	scores = getScores(jacob.mat,y,mu,env,coef)

	if(family=="negbin"){
		theta = coef[".theta"]
		score.theta = famFuns$scores.theta(theta,y,mu)
		scores = cbind(scores,score.theta)
		#DEPREC (theta-conditionned)
# 		isDummy = get("isDummy",env)
# 		if(isDummy){
# 			ll_d2 = famFuns$ll_d2(y,mu,coef)
# 			dxi_dbeta = deriv_xi(jacob.mat,ll_d2,env,coef)
# 			jacob.mat = jacob.mat + dxi_dbeta
# 		}
# 		theta = attr(mu,".theta")
# 		scores = scores + famFuns$scores_theta_part(theta,jacob.mat,y,mu)
	}

	return(scores)
}

ll_glm <- function(coef,env){
	# Log likelihood
	# cat("LL:\n") ; print(coef)
	# misc funs
	iter = get("iter",env) + 1
	assign("iter",iter,env)
	debug = get("debug",env)
	if(debug) cat("Iter",iter,"- Evaluation LL:",as.vector(coef),"\n")

	# computing the LL
	famFuns = get(".famFuns",env)
	family = get(".family",env)
	y <- get(".lhs",env)

	if(any(is.na(coef))) stop("Divergence... (some coefs are NA)\nTry option debug=TRUE to see the problem.")

	mu = get_mu(coef,env)

	#for the NEGBIN, we add the coef
	ll = famFuns$ll(y,mu,env,coef)

	if(debug) cat("LL =",ll,"\n")
	if(ll==(-Inf)) return(1e308)
	return(-ll) #je retourne -ll car la fonction d'optimisation minimise
}

ll_glm_TEST_score <- function(coef,env){
	#Used to compute the scores numerically
	#Not user oriented

	debug <- get("debug",env)
	if(debug) print(coef)

	# computing the LL
	famFuns = get(".famFuns",env)
	y <- get(".lhs",env)

	if(any(is.na(coef))) stop("Divergence... (some coefs are NA)\nTry option debug=TRUE to see the problem.")

	mu = get_mu(coef,env)

	ll = famFuns$ll_TEST_score(y,mu,env,coef)

	return(ll)
}

evalNLpart = function(coef,env){
	# cat("Enter evalNLpart : ",as.vector(coef),"\n")
	#fonction qui evalue la partie NL
	isNL = get("isNL",env)
	if(!isNL) return(0)

	nonlinear.params <- get("nonlinear.params",env)
	nl.call <- get(".nl.call",env)
	nbSave = get(".nbSave",env)
	nbMaxSave = get(".nbMaxSave",env)
	savedCoef = get(".savedCoef",env)
	savedValue = get(".savedValue",env)
	if(!is.null(names(coef))) coef = coef[nonlinear.params]
	else if (length(coef)!=length(nonlinear.params)) stop("Problem with the length of the NL coefficients.")

	if(nbMaxSave == 0){
		for(var in nonlinear.params) assign(var,coef[var],env)
		y_nl <- eval(nl.call,envir= env)
		return(y_nl)
	}

	for(i in nbSave:1){
		#les valeurs les + recentes sont en derniere position
		if(all(coef == savedCoef[[i]])){
			return(savedValue[[i]])
		}
	}

	#Si la valeur n'existe pas, on la sauvegarde
	#on met les valeurs les plus recentes en derniere position
	for(var in nonlinear.params) assign(var,coef[var],env)
	y_nl = eval(nl.call,envir = env)

	if(nbSave<nbMaxSave){
		savedCoef[[nbSave+1]] = coef
		savedValue[[nbSave+1]] = y_nl
		assign(".nbSave",nbSave+1,env)
	} else if(nbMaxSave>1){
		tmp = list()
		tmp[[nbSave]] = coef
		tmp[1:(nbSave-1)] = savedCoef[2:nbSave]
		savedCoef = tmp

		tmp = list()
		tmp[[nbSave]] = y_nl
		tmp[1:(nbSave-1)] = savedValue[2:nbSave]
		savedValue = tmp
	} else{
		savedCoef = list(coef)
		savedValue = list(y_nl)
	}

	# cat("computed NL part:",as.vector(coef),"\n")

	assign(".savedCoef",savedCoef,env)
	assign(".savedValue",savedValue,env)
	return(y_nl)
}

get_mu = function(coef,env){
	#This function computes the RHS of the equation
	#mu_L => to save one matrix multiplication
	isNL = get("isNL",env)
	isLinear = get("isLinear",env)
	isDummy = get("isDummy",env)
	nobs = get("nobs",env)
	params = get("params",env)
	family = get(".family",env)
	names(coef) = params

	#for managing mu:
	coefMu = get(".coefMu",env)
	valueMu = get(".valueMu",env)
	wasUsed = get(".wasUsed",env)
	if(wasUsed){
		coefMu = valueMu = list()
		assign(".wasUsed",FALSE,env)
	}
	if(length(coefMu)>0) for(i in 1:length(coefMu)) if(all(coef==coefMu[[i]])) return(valueMu[[i]])

	if(isNL){
		muNL = evalNLpart(coef,env)
	} else muNL = 0

	if(isLinear){
		linear.params = get("linear.params",env)
		linear.mat = get("linear.mat",env)
		mu_L = c(linear.mat%*%coef[linear.params])
	} else mu_L = 0

	mu_noDum = muNL + mu_L

	if(isDummy){
		#we get back the last dummy
		mu_dummies = getDummies(mu_noDum,env,coef)
	} else mu_dummies = 0

	mu = mu_noDum + mu_dummies
	if(isDummy) attr(mu,"mu_noDum") = mu_noDum

	if(length(mu)==0) mu = rep(mu,nobs)

	#DEPREC (theta-conditionned)
# 	if(family == "negbin"){
# 		#If isDummy, theta has already been computed with the dummies and been
# 		#sent to the environment
# 		theta = get(".theta",env)
# 		if(!isDummy){
# 			famFuns =  get(".famFuns",env)
# 			y = get(".lhs",env)
# 			# theta = qNR(theta,famFuns$ratio_fx_dfx_theta,y=y,mu=mu_noDum)
# 			#cat("before, theta=",theta,"\n")
# 			theta = famFuns$get_theta(theta,y,mu)
# 			#cat("theta=",theta,"\n")
# 			assign(".theta",theta,env)
# 		}
# 		attr(mu,".theta") = theta
# 	}

	#we save the value of mu:
	coefMu = append(coefMu,list(coef))
	valueMu = append(valueMu,list(mu))
	assign(".coefMu",coefMu,env)
	assign(".valueMu",valueMu,env)

	return(mu)
}

get_savedMu = function(coef,env){
	#This function gets he mu without computation
	#It follows a LL evaluation
	coefMu = get(".coefMu",env)
	valueMu = get(".valueMu",env)
	assign(".wasUsed",TRUE,env)

	if(length(coefMu)>0) for(i in 1:length(coefMu)) if(all(coef==coefMu[[i]])){
		# cat("coef nb:",i,"\n")
		return(valueMu[[i]])
	}

	stop("Problem in \"get_savedMu\":\n gradient did not follow LL evaluation.")
}

get_Jacobian = function(coef,env){
		# retrieves the Jacobian of the "rhs"
		params <- get("params",env)
		names(coef) <- params
		isNL <- get("isNL",env)
		isLinear <- get("isLinear",env)
		isGradient = get("isGradient",env)

		if(isNL){
			nonlinear.params = get("nonlinear.params",env)
			jacob.mat = get_NL_Jacobian(coef[nonlinear.params],env)
		} else jacob.mat = c()

		if(isLinear){
			linear.mat = get("linear.mat",env)
			jacob.mat = cbind(jacob.mat,linear.mat)
		}

		return(jacob.mat)
}

get_NL_Jacobian = function(coef,env){
	# retrieves the Jacobian of the non linear part
	#cat("In NL JAC:\n")
	#print(coef)
	nbSave = get(".JC_nbSave",env)
	nbMaxSave = get(".JC_nbMaxSave",env)
	savedCoef = get(".JC_savedCoef",env)
	savedValue = get(".JC_savedValue",env)

	nonlinear.params <- get("nonlinear.params",env)
	coef = coef[nonlinear.params]

	if(nbSave>0) for(i in nbSave:1){
		#les valeurs les + recentes sont en derniere position
		if(all(coef == savedCoef[[i]])){
			# cat("Saved value:",as.vector(coef),"\n")
			return(savedValue[[i]])
		}
	}

	#Si la valeur n'existe pas, on la sauvegarde
	#on met les valeurs les plus recentes en derniere position
	isGradient <- get("isGradient",env)
	if(isGradient){
		call_gradient <- get(".call_gradient",env)
		#we send the coef in the environment
		for(var in nonlinear.params) assign(var,coef[var],env)
		jacob.mat <- eval(call_gradient,envir=env)
		jacob.mat <- as.matrix(as.data.frame(jacob.mat[nonlinear.params]))
	} else {
		jacobian.method <- get("jacobian.method",env)
		jacob.mat <- numDeriv::jacobian(evalNLpart,coef,env=env,method=jacobian.method)
	}

	#Controls:
	if(anyNA(jacob.mat)){
		qui <- which(apply(jacob.mat,2,function(x) anyNA(x)))
		variables <- nonlinear.params[qui]
		stop("ERROR: The Jacobian of the nonlinear part has NA!\nThis concerns the following variables:\n",paste(variables,sep=" ; "))
	}

	#Sauvegarde
	if(nbSave<nbMaxSave){
		savedCoef[[nbSave+1]] = coef
		savedValue[[nbSave+1]] = jacob.mat
		assign(".JC_nbSave",nbSave+1,env)
	} else if(nbMaxSave>1){
		tmp = list()
		tmp[[nbSave]] = coef
		tmp[1:(nbSave-1)] = savedCoef[2:nbSave]
		savedCoef = tmp

		tmp = list()
		tmp[[nbSave]] = jacob.mat
		tmp[1:(nbSave-1)] = savedValue[2:nbSave]
		savedValue = tmp
	} else{
		savedCoef = list(coef)
		savedValue = list(jacob.mat)
	}

	# print(colSums(jacob.mat))

	# cat("computed NL Jacobian:",as.vector(coef),"\n")
	# print(savedCoef)

	assign(".JC_savedCoef",savedCoef,env)
	assign(".JC_savedValue",savedValue,env)
	return(jacob.mat)
}

ll0_nlglm <- function(lhs,env){
	#I have the closed form of the ll0
	famFuns = get(".famFuns",env)
	family = get(".family",env)
	y = get(".lhs",env)

	if(family=="negbin"){
		start = c(0,1)
		lower = c(-Inf,1e-4)
	} else {
		start = 0
		lower = NULL
	}

	opt <- nlminb(start=start,objective=famFuns$ll0,y=y,gradient=famFuns$grad0,lower=lower)
	return(list(loglik=-opt$objective,constant=opt$par[1]))
}

getGradient = function(jacob.mat,y,mu,env,coef,...){
	famFuns = get(".famFuns",env)
	ll_dl = famFuns$ll_dl(y=y,mu=mu,coef=coef,env=env)
	c(crossprod(jacob.mat, ll_dl))
}

getScores = function(jacob.mat,y,mu,env,coef,...){
	famFuns = get(".famFuns",env)
	isDummy = get("isDummy",env)

	ll_dl = famFuns$ll_dl(y=y,mu=mu,coef=coef,env=env)
	scores = jacob.mat* ll_dl

	if(isDummy){
		ll_d2 = famFuns$ll_d2(y=y,mu=mu,coef=coef,env=env)
		dxi_dbeta = deriv_xi(jacob.mat,ll_d2,env,coef)
		scores = scores + dxi_dbeta * ll_dl
	}

	return(as.matrix(scores))
}

getHessianLinear = function(jacob.mat,y,mu,env,coef,...){
	isDummy = get("isDummy",env)
	famFuns = get(".famFuns",env)
	ll_d2 = famFuns$ll_d2(y,mu,coef,env)
	H = crossprod(jacob.mat, jacob.mat * ll_d2)
	if(isDummy){
		dxi_dbeta = deriv_xi(jacob.mat,ll_d2,env,coef)
		H = H + crossprod(dxi_dbeta,jacob.mat * ll_d2)
# 		S_all = get(".S",env)
# 		g = length(S_all)
# 		if(g==1){
# 			S = S_all[[1]]
# 			S_Jmu = S%*%(jacob.mat*ll_d2)
# 			S_mu = as.vector(S%*%ll_d2)
# 			H = crossprod(jacob.mat, jacob.mat * ll_d2) - crossprod(S_Jmu,S_Jmu/S_mu)
# 		} else if(g==2) {
# 			#Dans le cas de deux clusters, on a une forme fermee
# 			#de la derivee seconde
# 			H = crossprod(jacob.mat, jacob.mat * ll_d2)
# 			d_xi_d_beta = deriv_xi(jacob.mat,ll_d2,env)
# 			H = H + crossprod(d_xi_d_beta,jacob.mat * ll_d2)
# 		} else {
# 			H = crossprod(jacob.mat, jacob.mat * ll_d2)
# 			#a mettre en plus propre:
# 			f = function(x,env){
# 				mu_noDum = get_mu_noDum(x,env)
# 				dum = getDummies(mu_noDum,env,x)
# 				dum
# 			}
# 			d_xi_d_beta = numDeriv::jacobian(f,coef,env=env)
# 			H = H + crossprod(d_xi_d_beta,jacob.mat * ll_d2)
# 		}
# 	} else {
# 		H = crossprod(jacob.mat, jacob.mat * ll_d2 )
	}

	return(as.matrix(H))
}

getDummies = function(mu,env,coef){
	#function built to get all the dummy variables
	#We retrieve past dummies (that are likely to be good
	# starting values)
	mu_dummies = get(".savedDummy",env)
	family = get(".family",env)
	#DEPREC (theta-conditionned)
# 	if(family=="negbin"){
# 		famFuns =  get(".famFuns",env)
# 		theta = get(".theta",env)
# 		okTheta = FALSE
# 	} else okTheta = TRUE
# 	if(family == "tobit"){
# 		assign(".sigma",sigma,env)
# 	}

	dum_all = get(".dummy",env)
	S_all = get(".S",env)
	g = length(dum_all)
	mu_in = mu + mu_dummies
	for(iter in 1:100){
		# cat("iter",iter,"\n")
		mu0 = mu_in
		#DEPREC (theta-conditionned)
# 		if(family=="negbin"){
# 			#print(head(mu_in))
# 			#attr(mu_in,".theta") = theta
# 			#ll = famFuns$ll(y,mu_in,env,coef)
# 			#cat("in LL:",sprintf("%.10f",ll),"\n")
# 			#we need to find the dummies AND the theta that max the LL
# 			#at the same time!
# 			# theta = qNR(theta,famFuns$ratio_fx_dfx_theta,y=y,mu=mu_noDum)
# 			theta_old = theta
# 			theta = famFuns$get_theta(theta,y,mu_in)
# 			okTheta = abs(theta_old - (theta_old<-theta))<1e-5
# 			#cat("iter",iter,"theta=",theta,"\n")
# 			#we send theta in the environment in order to compute the dummies properly
# 			assign(".theta",theta,env)
# 		}
		for(i in 1:g){
			dum = dum_all[[i]]
			S = S_all[[i]]
			#get the dummies
			mu_dum = computeDummies(dum,S,mu_in,env,coef)
			#add them to the stack
			# cat("group",i,":",mu_dum,"\n")
			mu_in = mu_in + mu_dum[dum]
			mu_dummies = mu_dummies + mu_dum[dum]
		}
		diff <- max(abs(mu0-mu_in)) #max error
# 		print(diff)
		if(anyNA(mu_in)) stop("Dummies could not be computed.")
# 		if(diff<1e-6 & okTheta) break
		if(diff<1e-6) break
		#if there's only one group, no need for 2 iterations
# 		if(g==1 & family!="negbin") break
		if(g==1) break
	}
	# print(iter)
	if(iter==100) warning("[getting dummies] iteration limit reached.",call. = FALSE,immediate. = TRUE)

	#we save the dummy:
	assign(".savedDummy",mu_dummies,env)
	#DEPREC (theta-conditionned)
# 	if(family=="negbin"){
# 		assign(".theta",theta,env)
# 		attr(mu_dummies,".theta") = theta
# 	}

	mu_dummies
}

computeDummies = function(dum,S,mu,env,coef){
	family = get(".family",env)
	famFuns = get(".famFuns",env)
	y = get(".lhs",env)

	#if family is poisson or gaussian, there is a closed form
	if(family%in%c("poisson","gaussian")){
		return(famFuns$closedFormDummies(dum,S,y,mu,env))
	}

	#For non closed-form dummies:
	#Basic Newton Raphson
	init = famFuns$initDummy(S,y,mu,env,coef)
	x1 = qNR(init,famFuns$ratio_fx_dfx,dum=dum,S=S,y=y,mu=mu,env=env,coef=coef)
	#x1 = rep(NA,nrow(S))
	if(anyNA(x1)){
		#It means that the basic N-R diverged for some dummies (can happen and it's perfectly normal)
		#We combine it with dichotomy to solve the problem
		quiNA = which(is.na(x1))
		idNA = which(dum%in%quiNA)

		S = S[quiNA,idNA,drop=FALSE]
		y = y[idNA]
		dum = unclass(factor(dum[idNA])) # necessary
		mu = mu[idNA]

		x1_corrected = dichoNR(famFuns,S,y,dum,mu,env,coef)
		x1[quiNA] = x1_corrected
	}
	x1
}

get_mu_noDum = function(coef,env){
	#TO REMOVE
	#This function computes the RHS of the equation
	#mu_L => to save one matrix multiplication
	isNL <- get("isNL",env)
	isLinear <- get("isLinear",env)
	params <- get("params",env)
	names(coef) <- params

	if(isNL){
		muNL = evalNLpart(coef,env)
	} else muNL = 0

	if(isLinear){
		linear.params <- get("linear.params",env)
		linear.mat <- get("linear.mat",env)
		mu_L <- c(linear.mat%*%coef[linear.params])
	} else mu_L = 0

	mu_noDum = muNL + mu_L

	return(mu_noDum)
}

qNR = function(init,fx_dfx,...){
	#cat("Enter QNR:",init,"\n")
	#quick Newton-Raphson
	precision = 1e-5
	itermax = 100
	iter = 0
	ok = TRUE
	x1 = init
	while(ok){
		x0 = x1
		x1 = x0 - fx_dfx(x0,...)
		#cat("x1:",x1," diff:",max(abs(x0-x1),na.rm = TRUE),"fx_dfx=",fx_dfx(x0,...),"\n")
		if((iter <- iter+1) ==itermax) stop("The algorithm getting the dummies diverged.\nThe error is unknown.")
		# print(max(abs(x0-x1),na.rm = TRUE))
		if(max(abs(x0-x1),na.rm = TRUE)<precision) ok = FALSE
	}
	x1
}

dichoNR = function(famFuns,S,y,dum,mu,env,coef){
	# function that combines dichotomy and N-R
	#if there are NAs: the NR did not work (starting values likely to be at the tails of the sigmoid)
	#=> we combine NR and dichotomy. We are sure that there is a zero, so dichotomy will work 100%
	# 1st we get the boundaries

	# Get the init, the inf and sup boundaries
	sum_y = as.vector(S%*%y)
	n_group = rowSums(S)
	#cat("nb need dicho:",length(n_group),"\n")
	muMean = tapply(mu,dum,mean)
	muMin = tapply(mu,dum,min)
	muMax = tapply(mu,dum,max)
	init = famFuns$guessDummy(sum_y,n_group,muMean) # evalue a la moyenne
	#plus les mu sont eleves, plus la dummy doit etre petite pour compenser
	borne_inf = famFuns$guessDummy(sum_y,n_group,muMax) # Value > 0, (on fait moins muMax)
	borne_sup = famFuns$guessDummy(sum_y,n_group,muMin) # Value < 0
	precision = 1e-5 ; itermax = 100
	iter = 0 ; ok = TRUE
	x1 = init

	while(ok){

		# 1st step: initialisation des bornes
		Value = famFuns$dum_fx(x1,sum_y,S,mu,dum,y,coef,env)
		borne_inf[Value>0] = x1[Value>0]
		borne_sup[Value<0] = x1[Value<0]

		# cat("\n\nStart iter:")
		# cat("\nx1:",x1,"\nValues:",Value,"\n")
		# cat("borne inf:",borne_inf,"\nborne sup:",borne_sup,"\n")

		# 2nd step: NR iteration
		x0 = x1
		Derivee = famFuns$dum_dfx(x0,S,mu,dum,y,coef,env) # AJOUTER COEF (pour neg bin et y)
		x1 = x0 - Value / Derivee

		# cat("X new (NR):",x1,"\n")

		#3rd step: dichotomy (if necessary)
		# Update of the value if it goes out of the boundaries
		qui = which(x1>=borne_sup | x1<=borne_inf)
		x1[qui] = (borne_inf[qui] + borne_sup[qui])/2

		# cat("X new (dicho):",x1,"\n")

		if((iter<-iter+1)==itermax) stop("The algorithm getting the dummies diverged.\nThe error is unknown.")
		if(anyNA(x1)) stop("The algorithm getting the dummies diverged.\nThe error is unknown.")
		if(max(abs(x0-x1))<precision) ok = FALSE
	}
	x1
}

deriv_xi = function(jacob.mat,ll_d2,env,coef){
	#Derivee des dummies
	dum_all = get(".dummy",env)
	S_all = get(".S",env)
	g = length(S_all)

	if(g==1){
		S = S_all[[1]]
		dum = dum_all[[1]]
		S_Jmu = S%*%(jacob.mat*ll_d2)
		S_mu = as.vector(S%*%ll_d2)
		dxi_dbeta = - S_Jmu[dum,] / S_mu[dum]
	} else if(g==2){
		#xi: la somme des dummies
		dxi_dbeta = deriv_xi_2cluster(jacob.mat,ll_d2,env)
	} else{
		f = function(x,env){
			mu_noDum = get_mu_noDum(x,env)
			dum = getDummies(mu_noDum,env,x)
			dum
		}
		dxi_dbeta = numDeriv::jacobian(f,coef,env=env)
	}
	as.matrix(dxi_dbeta)
}

deriv_xi_2cluster = function(jacob.mat,ll_d2,env){
	#This works only for sets of two clusters, no more
	S_all = get(".S",env)
	dum_all = get(".dummy",env)
	famFuns = get(".famFuns",env)
	y = get(".lhs",env)

	#comme on va inverser une matrice, on cherche a ce quelle soit
	#la plus petite possible
	if(nrow(S_all[[1]])>nrow(S_all[[2]])) {
		i=1 ; j=2
	} else {
		i=2 ; j=1
	}

	S_i = S_all[[i]]
	dum_i = dum_all[[i]]
	S_i_ll = as.vector(S_i%*%ll_d2)
	n_i = nrow(S_i)

	S_t = S_all[[j]]
	dum_t = dum_all[[j]]
	S_t_ll = as.vector(S_t%*%ll_d2)
	n_t = nrow(S_t)

	PHI = ll_d2 / S_t_ll[dum_t]
	THETA = ll_d2 / S_i_ll[dum_i]

	PHI_mat = matrix(0,nrow = n_t,ncol = n_i)
	PHI_mat[cbind(dum_t,dum_i)] = PHI

	THETA_mat = matrix(0,nrow = n_i,ncol = n_t)
	THETA_mat[cbind(dum_i,dum_t)] = THETA

	A_t = - S_t%*%(jacob.mat*PHI)
	B_t = PHI_mat%*%(S_i%*%(jacob.mat*THETA))
	invmat = MASS::ginv(diag(n_t)-PHI_mat%*%THETA_mat)
	d_gamma_t = invmat%*%(A_t+B_t)

	#On reconstruit les xi a partir des d_gamma_t
	dxi_hat = -(S_i%*%( (jacob.mat+d_gamma_t[dum_t,])*THETA ))[dum_i,] + d_gamma_t[dum_t,]

	return(dxi_hat)
}

deriv_xi_other = function(ll_d2,ll_dx_dother,env,coef){
	# derivative of the dummies wrt an other parameter
	family = get(".family",env)
	dum_all = get(".dummy",env)
	S_all = get(".S",env)
	g = length(S_all)

	if(g==1){
		S = S_all[[1]]
		dum = dum_all[[1]]
		S_Jmu = S%*%ll_dx_dother
		S_mu = as.vector(S%*%ll_d2)
		dxi_dother = - S_Jmu[dum,] / S_mu[dum]
	} else if(g==2){
		#xi: la somme des dummies
		dxi_dother = deriv_xi_2cluster_other(ll_d2,ll_dx_dother,env)
	} else{
		if(family=="negbin") other = ".theta"
		else if(family=="tobit") other = ".sigma"
		f = function(x,env,coef){
			coef[other] = x
			mu_noDum = get_mu_noDum(coef,env)
			dum = getDummies(mu_noDum,env,coef)
			dum
		}
		dxi_dother = numDeriv::jacobian(f,coef[other],env=env,coef=coef)
	}
	as.matrix(dxi_dother)
}

deriv_xi_2cluster_other = function(ll_d2,ll_dx_dother,env){
	#This works only for sets of two clusters, no more
	# For the
	S_all = get(".S",env)
	dum_all = get(".dummy",env)
	famFuns = get(".famFuns",env)
	y = get(".lhs",env)

	#comme on va inverser une matrice, on cherche a ce quelle soit
	#la plus petite possible
	if(nrow(S_all[[1]])>nrow(S_all[[2]])) {
		i=1 ; j=2
	} else {
		i=2 ; j=1
	}

	#Les matrices d'indices
	S_i = S_all[[i]]
	dum_i = dum_all[[i]]
	S_i_ll = as.vector(S_i%*%ll_d2)
	n_i = nrow(S_i)

	S_t = S_all[[j]]
	dum_t = dum_all[[j]]
	S_t_ll = as.vector(S_t%*%ll_d2)
	n_t = nrow(S_t)

	#Les matrices principales
	PHI = ll_d2 / S_t_ll[dum_t]
	MU = ll_d2 / S_i_ll[dum_i]

	PHI_mat = matrix(0,nrow = n_t,ncol = n_i)
	PHI_mat[cbind(dum_t,dum_i)] = PHI

	MU_mat = matrix(0,nrow = n_i,ncol = n_t)
	MU_mat[cbind(dum_i,dum_t)] = MU

	# Les memes avec la derivee croisee
	PHI_other = ll_dx_dother / S_t_ll[dum_t]
	MU_other = ll_dx_dother / S_i_ll[dum_i]

	PHI_other_mat = matrix(0,nrow = n_t,ncol = n_i)
	PHI_other_mat[cbind(dum_t,dum_i)] = PHI_other

	MU_other_mat = matrix(0,nrow = n_i,ncol = n_t)
	MU_other_mat[cbind(dum_i,dum_t)] = MU_other

	A_t = - S_t%*%(PHI_other)
	B_t = PHI_mat%*%(S_i%*%(MU_other))
	invmat = MASS::ginv(diag(n_t) - PHI_mat%*%MU_mat)
	d_gamma_t = invmat%*%(A_t+B_t)

	#On reconstruit les xi a partir des d_gamma_t
	dxi_hat = -(S_i%*%(MU_other + d_gamma_t[dum_t,]*MU))[dum_i,] + d_gamma_t[dum_t,]

	return(dxi_hat)
}

gt=function(ptm) (proc.time() - ptm)[[3]]
