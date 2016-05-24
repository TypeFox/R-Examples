
# Faire aide => textable et summary

# ADD other diagnostics (like odd-ratio, APE, PAE for probit/logit + other R2)
# ADD model0 with dummies

# REegler le probleme des cov.unscaled

print.feNmlm <- function(x,sd=c("standard","white","cluster","twoway"),cluster,...){
	#Ajouter le sandwich estimator des variances dans le print
	sd.val = match.arg(sd)

	x = summary(x,sd.val,cluster,...)

	# if(!is.null(x$clusterNames)) cat("Standard errors clustered by \"",x$clusterNames[1],"\"\n",sep="")

	coeftable = x$coeftable
	cat("Non linear glm, family =",x$family,"\n")
	cat("Observations:",x$n,"\n")
	if(!is.null(x$clusterSize)) cat("Clusters size: ",paste0(x$clusterNames,": ",x$clusterSize,collapse=", "),"\n",sep="")
	cat("Standard errors type:",sd.val,"\n")
	# The matrix of coefficients
	if(x$family=="negbin") {
		stats::printCoefmat(coeftable[-nrow(coeftable),])
		cat("\nDispersion parameter: theta =",coeftable[".theta",1],"\n")
	} else if(x$family=="tobit"){
		stats::printCoefmat(coeftable[-nrow(coeftable),])
		cat("\nSigma =",coeftable[".sigma",1],"\n")
	} else stats::printCoefmat(coeftable)

	cat("\n# Evaluations:",x$iterations,"\n")
	bic <- -2*x$loglik+x$df*log(x$n)
	aic <- -2*x$loglik+2*x$df
	cat("Log-likelihood:",x$loglik,"\nBIC:",bic,"\nAIC:",aic,"\n")
	cat("Pseudo-R2:",x$pseudo_r2,"\n")
	cat("Squared Cor.:",x$sq.cor,"\n")
	cat("Convergence state:",x$message,"\n")
}

##

summary.feNmlm <- function(object,sd=c("standard","white","cluster","twoway"),cluster,dof_correction=TRUE,...){
	#computes the clustered SD and return the modified vcov and coeftable
	#computes the BIC/AIC
	x = object

	sd.val = match.arg(sd)

	x$bic <- -2*x$loglik+x$df*log(x$n)
	x$aic <- -2*x$loglik+2*x$df

	if(anyNA(x$cov.unscaled)){
		warning("Standard errors set to NA because of collinearity.",call. = FALSE)
		return(x)
	}

	if(sd.val == "white"){
		vcov = crossprod(x$scores%*%x$cov.unscaled)
	} else if(sd.val == "cluster"){

		if(missing(cluster) && !is.null(x$id_dummies)){
			cluster = x$id_dummies[1] #we cluster with the first dummy
			if(length(x$clusterNames)>1) warning("Standard-errors clustered w.r.t. ",x$clusterNames[1],call. = FALSE,immediate. = TRUE)
		}
		if(missing(cluster)) stop("To display clustered standard errors, you must provide the cluster.")
		if(! (is.list(cluster) && length(cluster)==1 && length(cluster)!=x$n) ) stop("The 'cluster' must be a list containing one element being the vector of IDs of each observation.")

		vcov = vcovClust(cluster[[1]],x$cov.unscaled,x$scores,dof_correction)

	} else if(sd.val=="twoway"){

		if(missing(cluster) & !is.null(x$id_dummies)) cluster = x$id_dummies[1:2]
		if(missing(cluster)) stop("To display twoway standard errors, you must provide the clusters.")
		if(!is.list(cluster) && length(cluster)==2) stop("To display twoway standard errors, you must provide the clusters.\nThe 'cluster' argument must be a list with two elements.")
		if(!all(sapply(cluster,length)==x$n)) stop("The length of the elements of the argument 'clusters' must be the same as in the data.")

		# The clustered vcovs
		vcv_clu1 = vcovClust(cluster[[1]],x$cov.unscaled,x$scores,dof_correction)
		vcv_clu2 = vcovClust(cluster[[2]],x$cov.unscaled,x$scores,dof_correction)
		vcv_clu12 = vcovClust(paste0(cluster[[1]],"_",cluster[[2]]),x$cov.unscaled,x$scores,dof_correction)

		# The two way vcov:
		vcov = vcv_clu1 + vcv_clu2 - vcv_clu12

	} else vcov=x$cov.unscaled

	# Eigenfix if the vcov is not semi-positive definite
	# TODO better
	if(any(diag(vcov)<0)){
		warning("Variance was Eigenfixed.")
		e = eigen(vcov)
		val = e$values
		if(is.complex(val)) val = Re(val)
		vect = e$vectors
		val[val<0] = 0
		vcov = Re(vect%*%diag(val)%*%t(vect))
	}

	sd = sqrt(diag(vcov))

	#The coeftable is modified accordingly
	coeftable = x$coeftable
	coeftable[,2:4] = cbind(sd,x$coef/sd,2*pnorm(-abs(x$coef/sd)))

	attr(coeftable,"type") = attr(vcov,"type") = attr(sd,"type") = sd.val
	x$cov.unscaled = vcov
	x$coeftable = coeftable
	x$sd = sd

	return(x)
}

##

res2tex <- function(...,sd=c("standard","white","cluster","twoway"),cluster,digits=4,pseudo=TRUE,title,sdBelow=TRUE,drop,order,dict,file,append=TRUE,convergence=TRUE){
	#drop: a vector of regular expressions
	#order: a vector of regular expressions
	#dict: a 'named' vector
	#file: a character string
	# mettre les etoiles 'a la stata', avec le sym

	# To reimplement:
	keepFactors = FALSE

	#Options:
	dots <- list(...)
	#problem handling
	for(i in 1:length(dots)) if(!"feNmlm"%in%class(dots[[i]])) dots[[i]] <- NULL
	if(length(dots)==0) stop("Not any proper model (feNmlm) as argument!")

	sdType = match.arg(sd)
	n_models <- length(dots)

	var_list <- coef_list <- coef_below <- sd_below <- list()
	depvar_list <- obs_list <- list()
	r2_list <- aic_list <- bic_list <- convergence_list <- list()

	if(!keepFactors){
		factorNames = c()
		isFactor = vector(mode = "list",n_models)
	}

	for(m in 1:n_models){
		x <- summary(dots[[m]],sd=sdType,cluster)
		#variable dependante:
		depvar <- gsub(" ","",as.character(x$nonlinear.fml)[[2]])
		#we change the name of the dep. vari. if requested
		#if(!is.null(dict) & depvar%in%names(dict)) depvar[qui] = dict[depvar[qui]]

		a <- x$coeftable
		class(a) <- NULL

		#on enleve les facteurs des variables a garder
		if(!keepFactors){
			fact = rownames(a)
			qui_drop = grepl("factor(",fact,fixed = TRUE)
			a = a[!qui_drop,]
			b = fact[qui_drop]
			c = sapply(b,function(x) strsplit(x,"factor(",fixed=TRUE)[[1]][2])
			d = sapply(c,function(x) strsplit(x,")",fixed=TRUE)[[1]][1])
			factor_var = unique(d)
			if(!is.null(x$clusterNames)) factor_var = c(factor_var,x$clusterNames,recursive=TRUE)
		} else factor_var = c()

		#on enleve les espaces dans les noms de variables
		var <- c(gsub(" ","",row.names(a)))

		coef = as.character(round(a[,1],digits))
		sd = as.character(round(a[,2],digits))
		pval = cut(a[,4],breaks = c(-1,0.01,0.05,0.1,1),labels = c("\\sym{***}","\\sym{**}","\\sym{*}",""))
		structured_coef = c(paste0(coef,pval," (",sd,")"))

		if(!keepFactors){
			lFactor = rep("YES",length(factor_var))
			names(lFactor) = factor_var
			isFactor[[m]] = lFactor
			factorNames = unique(c(factorNames,factor_var,recursive=TRUE))
		}

		#saving the infos
		var_list[[m]] <- var
		names(structured_coef) <- var
		coef_list[[m]] <- structured_coef
		if(sdBelow){
			cb = c(paste0(coef,pval))
			sb = c(paste0("(",sd,")"))
			names(cb) = names(sb) = var
			coef_below[[m]] = cb
			sd_below[[m]] = sb
		}

		#la depvar
		depvar_list[[m]] <- depvar

		#statistics
		#Pseudo-R2 // AIC // BIC // N
		n <- x$n
		obs_list[[m]] <- n
		convergence_list[[m]] = x$convStatus

		K <- x$df #nb params
		ll <- x$loglik
		bic_list[[m]] <- round(-2*ll+K*log(n),3)
		aic_list[[m]] <- round(-2*ll+2*K,3)
		if(pseudo){
			r2_list[[m]] <- round(x$pseudo_r2,5)
		}
	}

	#prompting the infos gathered

	#Starting the table
	start_table = paste0("\\begin{table}[htbp]\\centering\n\\caption{", ifelse(!missing(title),title,"no title"),"}\n")
	end_table = "\\end{table}"

	#intro and outro Latex tabular
	intro_latex <- paste0("\\begin{tabular}{l|",paste0(rep("c",n_models),collapse="|"),"}\n")
	outro_latex <- "\\end{tabular}\n"

	#1st line
	first_line <- paste0("Variables&",paste0(depvar_list,collapse="&"),"\\\\\n\\hline\n\\hline\n")

	#Coefficients, the tricky part
	coef_lines <- list()
	all_vars <- unique(c(var_list,recursive=TRUE))

	# dropping some coefs
	if(!missing(drop)){
		if(!is.character(drop)) stop("the arg. 'drop' must be a character vector of regular expression (see help regexp).")
		for(var2drop in drop) all_vars = all_vars[!grepl(var2drop,all_vars)]
	}

	#ordering the coefs
	if(!missing(order)){
		if(!is.character(order)) stop("the arg. 'order' must be a character vector of regular expression (see help regexp).")
		for(var2order in rev(order)){
			who = grepl(var2order,all_vars)
			all_vars = c(all_vars[who],all_vars[!who])
		}
	}

	#changing the names of the coefs
	aliasVars = all_vars
	names(aliasVars) = all_vars
	if(!missing(dict)){
		if(!is.character(dict)|| is.null(names(dict))) stop("the arg. 'dict' must be a named character vector.")
		qui = which(all_vars%in%names(dict))
		who = aliasVars[qui]
		aliasVars[qui] = dict[who]
	}

	coef_mat <- all_vars
	for(m in 1:n_models) coef_mat <- cbind(coef_mat,coef_list[[m]][all_vars])
	coef_mat[is.na(coef_mat)] <- "  "
	if(sdBelow){
		coef_lines = c()
		for(v in all_vars){
			myCoef = mySd= myLine = c()
			for(m in 1:n_models){
				myCoef = c(myCoef,coef_below[[m]][v])
				mySd = c(mySd,sd_below[[m]][v])
			}
			myCoef[is.na(myCoef)] = "  "
			mySd[is.na(mySd)] = "  "
			myCoef = paste0(aliasVars[v],"&",paste0(myCoef,collapse="&"))
			mySd = paste0("  &",paste0(mySd,collapse="&"))
			myLines = paste0(myCoef,"\\\\\n",mySd,"\\\\\n")
			coef_lines = c(coef_lines,myLines)
		}
		coef_lines = paste0(coef_lines,collapse="")
	} else {
		coef_lines = paste0(paste0(apply(coef_mat,1,paste0,collapse="&"),collapse="\\\\\n"),"\\\\\n")
	}

	#Factors (if needed)
	if(!keepFactors && length(factorNames)>0){
		dumIntro = paste0("\\hline\n\\emph{Dummies}& ",paste(rep(" ",n_models),collapse="&"),"\\\\\n")
		for(m in 1:n_models) {
			quoi = isFactor[[m]][factorNames]
			quoi[is.na(quoi)] = "NO"
			isFactor[[m]] = quoi
		}
		allFactors = matrix(c(isFactor,recursive=TRUE),nrow = length(factorNames))
		allFactors = cbind(factorNames,allFactors)
		factor_lines <- paste0(paste0(apply(allFactors,1,paste0,collapse="&"),collapse="\\\\\n"),"\\\\\n")
	} else {
		factor_lines = NULL
		dumIntro = NULL
	}

	#Fit statistics
	fit_part <- paste0("\\hline\n\\emph{Fit statistics}& ",paste(rep(" ",n_models),collapse="&"),"\\\\\n")
	#Misc
	info_aic <- paste0("AIC & ",paste(aic_list,collapse="&"),"\\\\\n")
	info_bic <- paste0("BIC & ",paste(bic_list,collapse="&"),"\\\\\n")
	info_obs <- paste0("N & ",paste(obs_list,collapse="&"),"\\\\\n")
	info_r2 <- paste0("Adj-pseudo $R^2$ &",paste(r2_list,collapse="&"),"\\\\\n")
	info_convergence = paste0("Convergence &",paste(convergence_list,collapse="&"),"\\\\\n")

	if(!pseudo) info_r2 <- ""
	if(!convergence) info_convergence = ""

	if(!missing(file)) sink(file = file,append = append)

	cat(paste0(start_table,
				  intro_latex,
				  first_line,
				  coef_lines,
				  dumIntro,
				  factor_lines,
				  fit_part,
				  info_obs,
				  info_convergence,
				  info_r2,
				  info_aic,
				  info_bic,
				  outro_latex,
				  end_table))

	if(!missing(file)) sink()

}

##

res2table <- function(...,sd=c("standard","white","cluster","twoway"),cluster,digits=4,pseudo=TRUE,sdBelow=TRUE,drop,order,convergence=TRUE){

	# To be reimplemented:
	keepFactors = FALSE

	#Options:
	dots <- list(...)
	#problem handling
	for(i in 1:length(dots)) if(!"feNmlm"%in%class(dots[[i]])) dots[[i]] <- NULL
	if(length(dots)==0) stop("Not any proper model (feNmlm) as argument!")

	sdType = match.arg(sd)
	n_models <- length(dots)

	var_list <- coef_list <- coef_below <- sd_below <- list()
	depvar_list <- obs_list <- list()
	r2_list <- aic_list <- bic_list <- convergence_list <- list()

	if(!keepFactors){
		factorNames = c()
		isFactor = vector(mode = "list",n_models)
	}

	for(m in 1:n_models){
		x = summary(dots[[m]],sd=sdType,cluster)
		#variable dependante:
		depvar <- gsub(" ","",as.character(x$nonlinear.fml)[[2]])

		a <- x$coeftable
		class(a) <- NULL

		#on enleve les facteurs des variables a garder
		if(!keepFactors){
			fact = rownames(a)
			qui_drop = grepl("factor(",fact,fixed = TRUE)
			a = a[!qui_drop,]
			b = fact[qui_drop]
			c = sapply(b,function(x) strsplit(x,"factor(",fixed=TRUE)[[1]][2])
			d = sapply(c,function(x) strsplit(x,")",fixed=TRUE)[[1]][1])
			factor_var = unique(d)
			if(!is.null(x$clusterNames)) factor_var = c(factor_var,x$clusterNames,recursive=TRUE)
		} else factor_var = c()

		#on enleve les espaces dans les noms de variables
		var <- c(gsub(" ","",row.names(a)))

		coef = as.character(round(a[,1],digits))
		sd = as.character(round(a[,2],digits))
		pval = cut(a[,4],breaks = c(-1,0.01,0.05,0.1,1),labels = c("***","**","*",""))
		structured_coef = c(paste0(coef,pval," (",sd,")"))

		if(!keepFactors){
			lFactor = rep("YES",length(factor_var))
			names(lFactor) = factor_var
			isFactor[[m]] = lFactor
			factorNames = unique(c(factorNames,factor_var,recursive=TRUE))
		}

		#saving the infos
		var_list[[m]] <- var
		names(structured_coef) <- var
		coef_list[[m]] <- structured_coef
		if(sdBelow){
			cb = c(paste0(coef,pval))
			sb = c(paste0("(",sd,")"))
			names(cb) = names(sb) = var
			coef_below[[m]] = cb
			sd_below[[m]] = sb
		}

		#la depvar
		depvar_list[[m]] <- depvar

		#statistics
		#Pseudo-R2 // AIC // BIC // N
		n <- x$n
		obs_list[[m]] <- n
		convergence_list[[m]] <- x$convStatus

		K <- x$df #nb params
		ll <- x$loglik
		bic_list[[m]] <- round(-2*ll+K*log(n),3)
		aic_list[[m]] <- round(-2*ll+2*K,3)
		if(pseudo){
			r2_list[[m]] <- round(x$pseudo_r2,5)
		}
	}

	# The coefficients
	all_vars <- unique(c(var_list,recursive=TRUE))

	# dropping some coefs
	if(!missing(drop)){
		if(!is.character(drop)) stop("the arg. 'drop' must be a character vector of regular expression (see help regexp).")
		for(var2drop in drop) all_vars = all_vars[!grepl(var2drop,all_vars)]
	}

	#ordering the coefs
	if(!missing(order)){
		if(!is.character(order)) stop("the arg. 'order' must be a character vector of regular expression (see help regexp).")
		for(var2order in rev(order)){
			who = grepl(var2order,all_vars)
			all_vars = c(all_vars[who],all_vars[!who])
		}
	}

	coef_mat <- all_vars
	for(m in 1:n_models) coef_mat <- cbind(coef_mat,coef_list[[m]][all_vars])
	coef_mat[is.na(coef_mat)] <- "  "
	res = coef_mat

	# Used to draw a line
	myLine = "______________________________________"
	longueur = apply(res,2,function(x) max(nchar(as.character(x))))
	theLine = sapply(longueur,function(x) sprintf("%.*s",x,myLine))

	if(!keepFactors && length(factorNames)>0){
		#TO REWRITE: from old code
		for(m in 1:n_models) {
			quoi = isFactor[[m]][factorNames]
			quoi[is.na(quoi)] = "NO"
			isFactor[[m]] = quoi
		}
		allFactors = matrix(c(isFactor,recursive=TRUE),nrow = length(factorNames))
		allFactors = cbind(factorNames,allFactors)
		factor_lines <- paste0(paste0(apply(allFactors,1,paste0,collapse="&"),collapse="\\\\\n"),"\\\\\n")

		myLine = "-------------------------------"
		# res = rbind(res,theLine)
		res = rbind(res,c("Dummies:",sprintf("%.*s",longueur[-1],myLine)))
		factmat = matrix(c(strsplit(strsplit(factor_lines,"\n")[[1]],"&"),recursive = TRUE),ncol=n_models+1,byrow=TRUE)
		factmat[,ncol(factmat)]=gsub("\\","",factmat[,ncol(factmat)],fixed = TRUE)
		res = rbind(res,factmat)
	}

	res <- rbind(res,theLine)
	res <- rbind(res,c("N",c(obs_list,recursive = TRUE)))
	res <- rbind(res,c("Convergence",c(convergence_list,recursive = TRUE)))
	res <- rbind(res,c("Adj-pseudo R2",c(r2_list,recursive = TRUE)))
	res <- rbind(res,c("AIC",c(aic_list,recursive = TRUE)))
	res <- rbind(res,c("BIC",c(bic_list,recursive = TRUE)))

	res <- as.data.frame(res)
	names(res) <- c("variables",paste0("model ",1:n_models))
	row.names(res) = res$variables
	res$variables = NULL

	return(res)
}


getFE = function(x){
	#x is a feglm object
	#This function retrieves the dummies

	#TODO:
	#	- add the possibility to set the reference

	xi = x$dummies
	if(is.null(xi)) stop("There is no dummy to be retrieved.")
	clustNames = x$clusterNames
	g = length(x$clusterNames)
	#if there are more than two dummies, we set some references
	n = length(xi)
	#The casqe of the last group:
	dum = unclass(factor(x$id_dummies[[g]]))
	k = max(dum)
	qui = sapply(1:k,function(x) which.max(dum==x))
	xi_all = rep(0,nlevels(x$id_dummies[[g]]))

	alpha = list()
	for(i in 1:g){
		if(i==g){
			alpha[[clustNames[i]]] = xi[qui] - xi_all
		}else{
			dum = unclass(factor(x$id_dummies[[i]]))
			k = max(dum)
			S = Matrix::Matrix(0,k,n,sparse = TRUE)
			S[cbind(dum,1:n)] = 1
			sum_xi = as.vector(S%*%xi)
			n_group = rowSums(S)
			alpha[[clustNames[i]]] = a = (sum_xi - sum_xi[1])/n_group
			xi_all = xi_all + a[dum[qui]]
		}
	}
	alpha
}

#
# To compute clustered standard errors
#

vcovClust <- function (cluster,myBread,scores,dof_correction=FALSE){
	# Internal function: no need for controls, they come beforehand
	#	- cluster: the vector of dummies
	#	- myBread: the naive variance covariance matrix
	# - scores
	#Note: if length(unique(cluster)) == n (i.e. White correction), then the dof are such that vcovClust is equivalent to vcovHC(res,type="HC1")

	n <- NROW(scores)
	k <- NCOL(scores)

	#Control for cluster type
	cluster <- unclass(as.factor(cluster))
	nid <- max(cluster)
	#ind: the matrix used to create the sum of the scores per cluster (avoid loops)
	ind = Matrix::Matrix(0,nid,n,sparse = TRUE)
	myIndex = cbind(cluster,1:n)
	ind[myIndex] = 1

	#Compute the chocolate_bar and then the sandwich:
	RightScores <- as.matrix(ind %*% scores)
	chocobar <- crossprod(RightScores)

	#Finite sample correction:
	if(dof_correction) dof  <- nid / (nid-1) * (n-1) / (n-k)
	else dof = 1

	return(myBread %*% chocobar %*% myBread * dof)
}

