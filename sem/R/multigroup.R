### multigroup SEMs  
# last modified J. Fox 2013-06-14

## model definition

multigroupModel <- function(..., groups=names(models), allEqual=FALSE){
	models <- list(...)
	if (length (models) == 1){
		if (missing(groups)) stop("group names are missing, yet only 1 model is defined")
		ngroups <- length(groups)
		model <- models[[1]]
		models <- vector(ngroups, mode="list")
		for (g in 1:ngroups) models[[g]] <- model
		names(models) <- groups
		if (!allEqual){
			for (g in 1:ngroups){
				model <- models[[g]]
				models[[g]][!is.na(model[, 2]), 2] <- paste(model[!is.na(model[, 2]), 2], ".", groups[g], sep="")
			}
		}
	}
	else {
		if (is.null(names(models)) && missing(groups)) groups <- paste("Group.", 1:length(models), sep="")
		names(models) <- groups
	}
	class(models) <- "semmodList"
	models
}

print.semmodList <- function(x, ...){
	cat("\nMultigroup Structural Equation Model\n")
	groups <- names(x)
	for (group in groups){
		cat("\n Group: ", group, "\n\n")
		print(x[[group]])
	}
	invisible(x)
}


## sem() method for semmodList objects

sem.semmodList <- function(model, S, N, data, raw=FALSE, fixed.x=NULL, robust=!missing(data), formula, group="Group", debug=FALSE, ...){
    data.out <- NULL
    if (missing(S)){
        if (missing(data)) stop("S and data cannot both be missing")
        data.df <- inherits(data, "data.frame")
        if (data.df && missing(group)) stop("S and group cannot both be missing")
        if (data.df){
            if (!is.factor(data[, group])) stop("Groups variable, ", group, ", is not a factor")
            levels <- levels(data[, group])
            if (missing(formula)) formula <- as.formula(paste("~ . -", group))
        }
        else {
            if (!all(sapply(data, function(d) inherits(d, "data.frame")))) stop("data must be a data frame or list of data frames")
            levels <- names(data)
            if (is.null(levels)) levels <- paste("Group", seq(along=data), sep=".")
            if (missing(formula)) formula <- ~ .
        }
        G <- length(levels)
        if (is.list(formula) && length(formula) != G) stop("number of formulas, ", length(formula), ", not equal to number of groups, ", G, sep="")
        if (is.list(formula)){    
            if (!all(names(model) == names(formula))) warning("names of groups (", paste(names(model), collapse=", "), 
                ") is not the same as names of formulas in formula argument (", 
                paste(names(formula), collapse=", "), ")")
        }
        S <- vector(G, mode="list")
        names(S) <- levels
        N <- numeric(G)
        data.out <- vector(G, mode="list")
        for (g in 1:G){
            data.group <- if (data.df) subset(data, subset = data[, group] == levels[g]) else data[[g]]
            N.all <- nrow(data.group)
            form <- if (is.list(formula)) formula[[g]] else formula
            data.group <- model.matrix(form, data=data.group)
            colnames(data.group)[colnames(data.group) == "(Intercept)"] <- "Intercept"
            N[g] <- nrow(data.group)
            if (N[g] < N.all) warning(N.all - N[g]," observations removed due to missingness in group ", levels[g])
            S[[g]] <- if (raw) rawMoments(data.group) else{
                data.group <- data.group[, colnames(data.group) != "Intercept"]
                cov(data.group)
            }
            data.out[[g]] <- data.group
        }
    }
    else G <- length(S)
    if (length(model) != G) stop("number of group models, ", length(model), ", not equal to number of moment/data matrices, ", G, sep="")
    pars <-  unique(na.omit(unlist(lapply(model, function(mod) mod[, 2]))))
    vars <- rams <- vector(length(model), mode="list")
    all.par.names <- character(0)
    all.pars <- numeric(0)
    for (i in 1:G){
        obs.variables <- colnames(S[[i]])
        mod <- model[[i]]
        if ((!is.matrix(mod)) | ncol(mod) != 3) stop("model argument must be a 3-column matrix")
        startvalues <- as.numeric(mod[, 3])
        par.names <- mod[, 2]
        n.paths <- length(par.names)
        heads <- from <- to <- rep(0, n.paths)
        for (p in 1:n.paths){
            path <- parse.path(mod[p, 1])
            heads[p] <- abs(path$direction)
            to[p] <- path$second
            from[p] <- path$first
            if (path$direction == -1) {
                to[p] <- path$first
                from[p] <- path$second
            }
        }
        ram <- matrix(0, n.paths, 5)
        all.vars <- unique(c(to, from))
        latent.vars <- setdiff(all.vars, obs.variables)
        not.used <- setdiff(obs.variables, all.vars)
        if (length(not.used) > 0){
            rownames(S[[i]]) <- colnames(S[[i]]) <- obs.variables
            obs.variables <- setdiff(obs.variables, not.used)
            S[[i]] <- S[[i]][obs.variables, obs.variables]
            data.out[[i]] <- data.out[[i]][, obs.variables]
            warning("The following observed variables are in the input covariance or raw-moment matrix for group ", i,
                " but do not appear in the model:\n",
                paste(not.used, collapse=", "), "\n")
        }
        vars[[i]] <- c(obs.variables, latent.vars)
        ram[,1] <- heads
        ram[,2] <- apply(outer(vars[[i]], to, "=="), 2, which)
        ram[,3] <- apply(outer(vars[[i]], from, "=="), 2, which)   
        par.nos <- apply(outer(pars, par.names, "=="), 2, which)
        if (length(par.nos) > 0)
            ram[,4] <- unlist(lapply(par.nos, function(x) if (length(x) == 0) 0 else x))
        ram[,5]<- startvalues
        colnames(ram) <- c("heads", "to", "from", "parameter", "start value")
        rams[[i]] <- ram
        all.pars <- c(all.pars, par.nos)
        all.par.names <- c(all.par.names, par.names)
    }
    all.pars <- unique(unlist(all.pars))
    all.par.names <- unique(na.omit(all.par.names))	
    class(rams) <- "msemmod"
    result <- sem(rams, S, N, group=group, groups=names(model), raw=raw, fixed.x=fixed.x, param.names=all.par.names[all.pars], var.names=vars, debug=debug, ...)
    result$semmodList <- model
    result$data <- if(missing(data)) NULL else data.out
    if (robust && !missing(data) && inherits(result, "msemObjectiveML")){
        res <- robustVcovMsem(result)
        result$robust.vcov <- res$vcov
        result$chisq.scaled <- res$chisq.scaled
        result$adj.objects <- res$adj.objects
    }
    result
}

parse.path <- function(path) {                                           
	path.1 <- gsub("-", "", gsub(" ","", path))
	direction <- if (regexpr("<>", path.1) > 0) 2 
			else if (regexpr("<", path.1) > 0) -1
			else if (regexpr(">", path.1) > 0) 1
			else stop(paste("ill-formed path:", path))
	path.1 <- strsplit(path.1, "[<>]")[[1]]
	list(first=path.1[1], second=path.1[length(path.1)], direction=direction)
}


## sem() method for msemmod objects

sem.msemmod <- function(model, S, N, start.fn=startvalues, group="Group", groups=names(model), raw=FALSE, fixed.x, param.names, var.names, debug=FALSE, analytic.gradient=TRUE, warn=FALSE,
    maxiter=5000, par.size = c("ones", "startvalues"), start.tol = 1e-06, start=c("initial.fit", "startvalues"), initial.maxiter=1000,
    optimizer = optimizerMsem, objective = msemObjectiveML, ...){
    par.size <- match.arg(par.size)
    start <- match.arg(start)
    G <- length(groups)
    if (length(model) != G || length(N) != G) 
        stop("inconsistent number of groups in model (", length(model), "), S (", G, "), and N (", length(N), ") arguments")
    if (is.null(names(S))) names(S) <- groups
    if (is.null(names(N))) names(N) <- groups
    if (is.null(names(model))) names(model) <- groups
    if (!all(groups == names(model))) warning("names of groups (", paste(groups, collapse=", "), 
        ") is not the same as names of models in model argument (", 
        paste(names(model), collapse=", "), ")")
    if (!all(groups == names(S))) warning("names of groups (", paste(groups, collapse=", "),
        ") is not the same as names of moment matrices in S argument (", 
        paste(names(S), collapse=", "), ")")
    if (!all(groups == names(N))) warning("names of groups (", paste(groups, collapse=", "),
        ") is not the same as names of sample sizes in N argument (", 
        paste(names(N), collapse=", "), ")")
    if (length(fixed.x) == 1) fixed.x <- lapply(1:G, function(g) fixed.x)
    n.fix <- 0 
    if (!is.null(fixed.x)){
        n.fix <- numeric(G)
        for (g in 1:G){
            fx <- fixed.x[[g]]
            n.fix[g] <- length(fx)
            if (n.fix[g] == 0) next
            fx <- which(rownames(S[[g]]) %in% fx)
            mod <- model[[g]]
            for (i in 1:n.fix[g]){
                for (j in 1:i){
                    mod <- rbind(mod, c(2, fx[i], fx[j],
                        0, S[[g]][fx[i], fx[j]]))
                }
            }
            model[[g]] <- mod
        }
    }
    t <- max(sapply(model, function(r) max(r[, 4])))
    if(missing(param.names)) param.names <- paste("Parameter", 1:t, sep=".")
    if (missing(var.names)) var.names <- lapply(model, function(mod) paste("Variable", 1:max(mod[, c(2,3)]), sep="."))
    n <- sapply(S, nrow)
    m <- sapply(model, function(r)  max(r[, c(2, 3)]))
    logdetS <- sapply(S, function(s) log(det(unclass(s))))
    sel.free.2 <- sel.free.1 <- arrows.2.free <- arrows.1.free <- arrows.2t <- arrows.2 <- arrows.1 <- 
        two.free <- one.free <- one.head <- sel.free <- fixed <- par.posn <- correct <- J <- vector(mode="list", length=G)  
    initial.iterations <- if (start == "initial.fit") numeric(G) else NULL
    for (g in 1:G){
        mod <- model[[g]]
        #         tt <- sum(mod[, 4] != 0)
        #         mod[mod[, 4] != 0, 4] <- 1:tt
        initial.pars <- mod[, 4]
        unique.pars <- unique(initial.pars)
        unique.pars <- unique.pars[unique.pars != 0]
        tt <- length(unique.pars)
        if (tt > 0) for (i in 1:tt) mod[initial.pars == unique.pars[i], 4] <- i
        startvals <- if (start == "initial.fit"){
            prelim.fit <- sem(mod, S[[g]], N=N[[g]], raw=raw, param.names=if(tt > 0) as.character(1:tt) else character(0), var.names=as.character(1:m[[g]]), 
                maxiter=initial.maxiter)
            initial.iterations[g] <- prelim.fit$iterations
            coef(prelim.fit)
        }
        else start.fn(S[[g]], mod)
        model[[g]][mod[, 4] != 0, 5] <- ifelse(is.na(mod[mod[, 4] != 0, 5]), startvals, mod[mod[, 4] != 0, 5])
        J[[g]] <- matrix(0, n[g], m[g])
        correct[[g]] <- matrix(2, m[g], m[g])
        diag(correct[[g]]) <- 1
        observed <- 1:n[g]
        J[[g]][cbind(observed, observed)] <- 1
        par.posn[[g]] <-  sapply(1:t, function(i) which(model[[g]][,4] == i)[1])
        colnames(model[[g]]) <- c("heads", "to", "from", "parameter", "start value")
        rownames(model[[g]]) <- rep("", nrow(model[[g]]))
        fixed[[g]] <- model[[g]][, 4] == 0
        sel.free[[g]] <- model[[g]][, 4]
        sel.free[[g]][fixed[[g]]] <- 1
        one.head[[g]] <- model[[g]][, 1] == 1
        one.free[[g]] <- which( (!fixed[[g]]) & one.head[[g]] )
        two.free[[g]] <- which( (!fixed[[g]]) & (!one.head[[g]]) )
        arrows.1[[g]] <- model[[g]][one.head[[g]], c(2, 3), drop=FALSE]
        arrows.2[[g]] <- model[[g]][!one.head[[g]], c(2, 3), drop=FALSE]
        arrows.2t[[g]] <- model[[g]][!one.head[[g]], c(3 ,2), drop=FALSE]
        arrows.1.free[[g]] <- model[[g]][one.free[[g]], c(2, 3), drop=FALSE]
        arrows.2.free[[g]] <- model[[g]][two.free[[g]], c(2, 3), drop=FALSE]
        sel.free.1[[g]] <- sel.free[[g]][one.free[[g]]]
        sel.free.2[[g]] <- sel.free[[g]][two.free[[g]]]
    }
    unique.free.1 <- lapply(sel.free.1, unique)
    unique.free.2 <- lapply(sel.free.2, unique)
    startvals <- numeric(t)
    for (j in 1:t) startvals[j] <- mean(unlist(sapply(model, function(r) r[r[, 4] == j, 5])), na.rm=TRUE)
    model.description <- list(G=G, m=m, n=n, t=t, fixed=fixed, ram=model, sel.free=sel.free, arrows.1=arrows.1, 
        one.head=one.head, arrows.2=arrows.2, arrows.2t=arrows.2t, J=J, S=S, logdetS=logdetS, 
        N=N, raw=raw, correct=correct, unique.free.1=unique.free.1, unique.free.2=unique.free.2, 
        arrows.1.free=arrows.1.free, arrows.2.free=arrows.2.free, param.names=param.names, 
        var.names=var.names)
    result <- optimizer(start=startvals, objective=objective, gradient=analytic.gradient,
        maxiter=maxiter, debug=debug, par.size=par.size, model.description=model.description, warn=warn, ...)
    if (!is.na(result$iterations)) if(result$iterations >= maxiter) warning("maximum iterations exceeded")
    result <- c(result, list(ram=model, param.names=param.names, var.names=var.names, group=group, groups=groups,
        S=S, N=N, J=J, n=n, m=m, t=t, raw=raw, optimizer=optimizer, objective=objective, fixed.x=fixed.x, n.fix=n.fix,
        initial.iterations=initial.iterations))
    cls <- gsub("\\.", "", deparse(substitute(objective)))
    cls <- gsub("2", "", cls)
    class(result) <- c(cls, "msem")
    result
}




## ML objective function for multigroup SEMs

msemObjectiveML2 <- function(gradient=TRUE){
	result <- list(
			objective = function(par, model.description){
				with(model.description, {
							f <- numeric(G)
							AA <- PP <- CC <- vector(G, mode="list")
							grad.all <- if (gradient) rep(0, t) else NULL
							for (g in 1:G){
								A <- P <- matrix(0, m[g], m[g])
								val <- ifelse (fixed[[g]], ram[[g]][, 5], par[sel.free[[g]]])
								A[arrows.1[[g]]] <- val[one.head[[g]]]
								P[arrows.2t[[g]]] <- P[arrows.2[[g]]] <- val[!one.head[[g]]]
								I.Ainv <- solve(diag(m[g]) - A)
								C <- J[[g]] %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J[[g]])
								Cinv <- solve(C)
								f[g] <- sum(diag(S[[g]] %*% Cinv)) + log(det(C)) - n[[g]] - logdetS[g]
								CC[[g]] <- C
								AA[[g]] <- A
								PP[[g]] <- P
								if (gradient){
									grad.P <- correct[[g]] * t(I.Ainv) %*% t(J[[g]]) %*% Cinv %*% (C - S[[g]]) %*% Cinv %*% J[[g]] %*% I.Ainv
									grad.A <- grad.P %*% P %*% t(I.Ainv)        
									grad <- rep(0, t)
									grad[sort(unique.free.1[[g]])] <- tapply(grad.A[arrows.1.free[[g]]], ram[[g]][ram[[g]][,1]==1 & ram[[g]][,4]!=0, 4], sum)
									grad[sort(unique.free.2[[g]])] <- tapply(grad.P[arrows.2.free[[g]]], ram[[g]][ram[[g]][,1]==2 & ram[[g]][,4]!=0, 4], sum)
									grad.all <- grad.all + ((N[g] - (!raw))/(sum(N) - (!raw)*G))*grad
								}
							}
							ff <- f
							f <- sum((N - (!raw))*f)/(sum(N) - (!raw)*G)
							attributes(f) <- list(gradient=grad.all, A=AA, P=PP, C=CC, f=ff)
							f
						})
			}
	)
	if (gradient)
		result$gradient <- function(par, model.description){
			with(model.description, {
						grad.total <- rep(0, t)
						for (g in 1:G){
							A <- P <- matrix(0, m[g], m[g])
							val <- ifelse (fixed[[g]], ram[[g]][,5], par[sel.free[[g]]])
							A[arrows.1[[g]]] <- val[one.head[[g]]]
							P[arrows.2t[[g]]] <- P[arrows.2[[g]]] <- val[!one.head[[g]]]
							I.Ainv <- solve(diag(m[g]) - A)
							C <- J[[g]] %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J[[g]])
							Cinv <- solve(C)
							grad.P <- correct[[g]] * t(I.Ainv) %*% t(J[[g]]) %*% Cinv %*% (C - S[[g]]) %*% Cinv %*% J[[g]] %*% I.Ainv
							grad.A <- grad.P %*% P %*% t(I.Ainv)        
							grad <- rep(0, t)
							grad[sort(unique.free.1[[g]])] <- tapply(grad.A[arrows.1.free[[g]]], ram[[g]][ram[[g]][,1]==1 & ram[[g]][,4]!=0, 4], sum)
							grad[sort(unique.free.2[[g]])] <- tapply(grad.P[arrows.2.free[[g]]], ram[[g]][ram[[g]][,1]==2 & ram[[g]][,4]!=0, 4], sum)
							grad.total <- grad.total + ((N[g] - (!raw))/(sum(N) - (!raw)*G))*grad
						}
						grad.total
					})
		}
	class(result) <- "msemObjective"
	result
}

msemObjectiveML <- function(gradient=TRUE){
	result <- list(
			objective = function(par, model.description){
				with(model.description, {
							res <- msemCompiledObjective(par=par, model.description=model.description, objective="objectiveML")
							AA <- PP <- CC <- vector(G,  mode="list")
							for(g in 1:model.description$G)
							{
								AA[[g]] <- res$A[[g]]
								PP[[g]] <- res$P[[g]]
								CC[[g]] <- res$C[[g]]
							}
							
							f <- res$f
							attributes(f) <- list(gradient=res$gradient, A=AA, P=PP, C=CC, f=res$ff)
							f
						})
			}
	)
	if (gradient)
		result$gradient <- function(par, model.description){
			with(model.description, {
						
						res <- msemCompiledObjective(par=par, model.description=model.description, objective="objectiveML")
						res$gradient
					})
		}
	class(result) <- "msemObjective"
	result
}

msemObjectiveGLS <- function(gradient=FALSE){
	result <- list(
			objective = function(par, model.description){
				with(model.description, {
							
							res <- msemCompiledObjective(par=par, model.description=model.description, objective="objectiveGLS")
							AA <- PP <- CC <- vector(G,  mode="list")
							for(g in 1:model.description$G)
							{
								AA[[g]] <- res$A[[g]]
								PP[[g]] <- res$P[[g]]
								CC[[g]] <- res$C[[g]]
							}
							
							f <- res$f
							attributes(f) <- list(A=AA, P=PP, C=CC, f=res$ff)
							f
						})
			}
	)
	
	class(result) <- "msemObjective"
	result
}

msemObjectiveFIML <- function(gradient=FALSE){
	result <- list(
			objective = function(par, model.description){
				with(model.description, {
							
							res <- msemCompiledObjective(par=par, model.description=model.description, objective="objectiveFIML")
							AA <- PP <- CC <- vector(G,  mode="list")
							for(g in 1:model.description$G)
							{
								AA[[g]] <- res$A[[g]]
								PP[[g]] <- res$P[[g]]
								CC[[g]] <- res$C[[g]]
							}
							
							f <- res$f
							attributes(f) <- list(A=AA, P=PP, C=CC, f=res$ff)
							f
						})
			}
	)
	
	class(result) <- "msemObjective"
	result
}

##  nlm()-based optimizer for multigroup SEMs

optimizerMsem <- function(start, objective=msemObjectiveML, gradient=TRUE,
		maxiter, debug, par.size, model.description, warn=FALSE, ...){
	with(model.description, {
				obj <- objective(gradient=gradient)$objective
				typsize <- if (par.size == 'startvalues') abs(start) else rep(1, t)
				
				if(identical(objective, msemObjectiveML)) objectiveCompiled <- "objectiveML"
				else if (identical(objective, msemObjectiveGLS)) {
					objectiveCompiled <- "objectiveGLS"
					gradient <- FALSE
				}
				else if (identical(objective, msemObjectiveFIML)) {
					objectiveCompiled <- "objectiveFIML"
					gradient <- FALSE
				}
				else stop("optimizerMsem requires the msemObjectiveML, msemObjectiveGLS or msemObjectiveFIML objective function")
				
				if (!warn) save.warn <- options(warn=-1)
				
				res <- msemCompiledSolve(model.description=model.description, start=start, objective=objectiveCompiled, 
						typsize=typsize, debug=debug, maxiter=maxiter)
				
				if (!warn) options(save.warn)
				result <- list(covergence=NULL, iterations=NULL, coeff=NULL, vcov=NULL, criterion=NULL, C=NULL, A=NULL, P=NULL)
				result$convergence <- res$code <= 2
				result$iterations <- res$iterations
				par <- res$estimate
				names(par) <- param.names
				result$coeff <- par
				if (!result$convergence)
					warning(paste('Optimization may not have converged; nlm return code = ',
									res$code, '. Consult ?nlm.\n', sep=""))
				vcov <- matrix(NA, t, t)
				qr.hess <- try(qr(res$hessian), silent=TRUE)
				if (class(qr.hess) == "try-error"){
					warning("Could not compute QR decomposition of Hessian.\nOptimization probably did not converge.\n")
				}
				else if (qr.hess$rank < t){
					warning(' singular Hessian: model is probably underidentified.\n')
					which.aliased <- qr.hess$pivot[-(1:qr.hess$rank)]
					result$aliased <- param.names[which.aliased]
				}
				else {
					df <- sum(N) - (!raw)*G
					vcov <-  (2/df) * solve(res$hessian)
					if (any(diag(vcov) < 0)) {
						result$aliased <- param.names[diag(vcov) < 0]
						warning("Negative parameter variances.\nModel may be underidentified.\n")
					}
				}
				colnames(vcov) <- rownames(vcov) <- param.names
				result$vcov <- vcov
				result$criterion <- res$minimum
				C <- res$C
				A <- res$A
				P <- res$P
				for (g in 1:G){
					rownames(C[[g]]) <- colnames(C[[g]]) <-rownames(S[[g]])
					rownames(A[[g]]) <- colnames(A[[g]]) <- var.names[[g]]
					rownames(P[[g]]) <- colnames(P[[g]]) <- var.names[[g]]
				}
				result$C <- C
				result$A <- A
				result$P <- P
				result$group.criteria <- res$ff
				class(result) <- "msemResult"
				result
			})
}

msemOptimizerNlm <- function(start, objective=msemObjectiveML, gradient=TRUE,
		maxiter, debug, par.size, model.description, warn=FALSE, ...){
	with(model.description, {
				obj <- objective(gradient=gradient)$objective
				typsize <- if (par.size == 'startvalues') abs(start) else rep(1, t)
				if (!warn) save.warn <- options(warn=-1)
				res <- nlm(obj, start, iterlim=maxiter, print.level=if(debug) 2 else 0,
						typsize=typsize, hessian=TRUE, model.description=model.description, ...)
				if (!warn) options(save.warn)
				result <- list(covergence=NULL, iterations=NULL, coeff=NULL, vcov=NULL, criterion=NULL, C=NULL, A=NULL, P=NULL)
				result$convergence <- res$code <= 2
				result$iterations <- res$iterations
				par <- res$estimate
				names(par) <- param.names
				result$coeff <- par
				if (!result$convergence)
					warning(paste('Optimization may not have converged; nlm return code = ',
									res$code, '. Consult ?nlm.\n', sep=""))
				vcov <- matrix(NA, t, t)
				qr.hess <- try(qr(res$hessian), silent=TRUE)
				if (class(qr.hess) == "try-error"){
					warning("Could not compute QR decomposition of Hessian.\nOptimization probably did not converge.\n")
				}
				else if (qr.hess$rank < t){
					warning(' singular Hessian: model is probably underidentified.\n')
					which.aliased <- qr.hess$pivot[-(1:qr.hess$rank)]
					result$aliased <- param.names[which.aliased]
				}
				else {
					df <- sum(N) - (!raw)*G
					vcov <-  (2/df) * solve(res$hessian)
					if (any(diag(vcov) < 0)) {
						result$aliased <- param.names[diag(vcov) < 0]
						warning("Negative parameter variances.\nModel may be underidentified.\n")
					}
				}
				colnames(vcov) <- rownames(vcov) <- param.names
				result$vcov <- vcov
				result$criterion <- res$minimum
				obj <- obj(par, model.description)
				C <- attr(obj, "C")
				A <- attr(obj, "A")
				P <- attr(obj, "P")
				for (g in 1:G){
					rownames(C[[g]]) <- colnames(C[[g]]) <-rownames(S[[g]])
					rownames(A[[g]]) <- colnames(A[[g]]) <- var.names[[g]]
					rownames(P[[g]]) <- colnames(P[[g]]) <- var.names[[g]]
				}
				result$C <- C
				result$A <- A
				result$P <- P
				result$group.criteria <- attr(obj, "f")
				class(result) <- "msemResult"
				result
			})
}

# methods for msem and msemObjectiveML objects


print.msemObjectiveML <- function(x, ...){
	n <- x$n
	n.fix <- x$n.fix
	df <- sum(n*(n + 1)/2) - x$t - sum(n.fix*(n.fix + 1)/2)
	chisq <- sum(x$N - !x$raw)*x$criterion
	cat("\n Model Chisquare =", chisq, " Df =", df, "\n\n")
	print(x$coeff)
	invisible(x)
}

print.msemObjectiveGLS <- function(x, ...) print.msemObjectiveML(x, ...)
print.msemObjectiveFIML <- function(x, ...) print.msemObjectiveML(x, ...)

summary.msemObjectiveML <- function(object, digits=getOption("digits"), conf.level=.90, robust=FALSE, 
		analytic.se=object$t <= 500,
        fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc", "BIC"),
        ...){
    fit.indices <- if (is.null(fit.indices)) ""
    else {
        if (missing(fit.indices)){
            if (is.null(opt <- getOption("fit.indices"))) c("AIC", "BIC") else opt
        }
        else match.arg(fit.indices, several.ok=TRUE)
    }
	if(inherits(object, "msemObjectiveGLS")) analytic.se <- FALSE
	else if(inherits(object, "msemObjectiveFIML")) analytic.se <- FALSE
	groups <- object$groups
	G <- length(groups)
	par <- object$coeff
	vcov <- vcov(object, robust=robust, analytic.se=analytic.se)
	n <- object$n
	m <- object$m
	S <- object$S
	C <- object$C
	A <- object$A
	P <- object$P
	N <- object$N
	J <- object$J
	n.fix <- object$n.fix
	if (length(n.fix) == 1) n.fix <- rep(n.fix, G)
	group.criteria <- object$group.criteria
	var.names <- object$var.names
	param.names <- object$param.names
	semmod <- object$semmodList
	ram <- object$ram
	group.summaries <- list(G, model="list")
	for (g in 1:G){
		par.names <- param.names[ram[[g]][, 4]]
		par.gr <- par[par.names]
		group <- list(coeff=par.gr, vcov=vcov[par.names, par.names], n.fix=n.fix[g], n=n[[g]], m=m[[g]], S=S[[g]], C=C[[g]], N=N[[g]],
				t=length(par.gr), raw=object$raw, var.names=var.names[[g]], ram=ram[[g]],
				J=J[[g]], A=A[[g]], P=P[[g]], criterion=group.criteria[g], par.posn=ram[[g]][, 4] != 0, 
				iterations=object$iterations, semmod=semmod[[g]], adj.obj=object$adj.objects[[g]], 
				robust.vcov=object$robust.vcov[par.names, par.names])
		class(group) <- if(inherits(object, "msemObjectiveGLS")) c("objectiveGLS", "sem") 
				else if(inherits(object, "msemObjectiveFIML")) c("objectiveFIML", "sem")
				else c("objectiveML", "sem")
		group.summaries[[g]] <- if(inherits(object, "msemObjectiveGLS"))
					summary(group, digits=digits, conf.level=conf.level, robust=FALSE, fit.indices=fit.indices, ...)
			else if(inherits(object, "msemObjectiveFIML")) 
					summary(group, digits=digits, conf.level=conf.level, robust=FALSE, fit.indices=fit.indices, ...)
			else summary(group, digits=digits, conf.level=conf.level, robust=robust, analytic.se=FALSE, fit.indices=fit.indices, ...)
		group.summaries[[g]]$iterations <- NA
	}
	df <- sum(n*(n + 1)/2) - object$t - sum(n.fix*(n.fix + 1)/2)
	chisq <- sum(N - !object$raw)*object$criterion
	wt <- (N - object$raw)/(sum(N - object$raw))
	SRMR <-if (!object$raw && "SRMR" %in% fit.indices) sum(wt*sapply(group.summaries, function(s) s$SRMR)) else NA
	GFI <- if (!object$raw && "GFI" %in% fit.indices) sum(wt*sapply(group.summaries, function(s) s$GFI)) else NA
	chisqNull <- if(!object$raw) sum(sapply(group.summaries, function(s) s$chisqNull)) else NA
	dfNull <- sum(n*(n - 1)/2)
	if (!object$raw && df > 0){
		AGFI <- if ("AGFI" %in% fit.indices) 1 - (sum(n*(n * 1))/(2*df))*(1 - GFI) else NA
		NFI <- if ("NFI" %in% fit.indices) (chisqNull - chisq)/chisqNull else NA
		NNFI <- if ("NNFI" %in% fit.indices) (chisqNull/dfNull - chisq/df)/(chisqNull/dfNull -1) else NA
        L1 <- max(chisq - df, 0)
    	L0 <- max(L1, chisqNull - dfNull)
        CFI <- if ("CFI" %in% fit.indices) 1 - L1/L0 else NA
        RNI <- if ("RNI" %in% fit.indices) 1 - (chisq - df)/(chisqNull - dfNull) else NA
        IFI <- if ("IFI" %in% fit.indices) (chisqNull - chisq)/(chisqNull - df) else NA
        if ("RMSEA" %in% fit.indices){
    		RMSEA <- sqrt(G*max(object$criterion/df - 1/(sum(N - 1)), 0))
    		tail <- (1 - conf.level)/2
    		max <- sum(N)
    		while (max > 1) {
    			res <- optimize(function(lam) (tail - pchisq(chisq, 
    											df, ncp = lam))^2, interval = c(0, max))
    			if (is.na(res$objective) || res$objective < 0) {
    				max <- 0
    				warning("cannot find upper bound of RMSEA")
    				break
    			}
    			if (sqrt(res$objective) < tail/100) 
    				break
    			max <- max/2
    		}
    		lam.U <- if (max <= 1) 
    					NA
    				else res$minimum
    		max <- max(max, 1)
    		while (max > 1) {
    			res <- optimize(function(lam) (1 - tail - pchisq(chisq, 
    											df, ncp = lam))^2, interval = c(0, max))
    			if (sqrt(res$objective) < tail/100) 
    				break
    			max <- max/2
    			if (is.na(res$objective) || res$objective < 0) {
    				max <- 0
    				warning("cannot find lower bound of RMSEA")
    				break
    			}
    		}
    		lam.L <- if (max <= 1) 
    					NA
    				else res$minimum
    		RMSEA.U <- sqrt(G*lam.U/(sum(N - 1) * df))
    		RMSEA.L <- sqrt(G*lam.L/(sum(N - 1) * df))
            }
        else RMSEA.U <- RMSEA.L <- RMSEA <- NA
	}
	else RMSEA.U <- RMSEA.L <- RMSEA <- NFI <- NNFI <-IFI <- RNI <- CFI <- AGFI <- NA
	if (robust){
		chisq.adjusted <- sum(sapply(group.summaries, function(x) x$chisq.adjusted))
		if (!object$raw && df > 0){
			chisqNull.adjusted <- sum(sapply(group.summaries, function(x) x$chisqNull.adjusted))
			NFI.adjusted <- if ("NFI" %in% fit.indices) (chisqNull.adjusted - chisq.adjusted)/chisqNull.adjusted else NA
			NNFI.adjusted <- if ("NNFI" %in% fit.indices) (chisqNull.adjusted/dfNull - chisq.adjusted/df)/(chisqNull.adjusted/dfNull - 1) else NA
			L1 <- max(chisq.adjusted - df, 0)
			L0 <- max(L1, chisqNull.adjusted - dfNull)
			CFI.adjusted <- if ("CFI" %in% fit.indices) 1 - L1/L0 else NA
            RNI.adjusted <- if ("RNI" %in% fit.indices) 1 - (chisq.adjusted - df)/(chisqNull.adjusted - dfNull) else NA
            IFI.adjusted <- if ("IFI" %in% fit.indices) (chisqNull.adjusted - chisq.adjusted)/(chisqNull.adjusted - df) else NA
		}
	}
	if (object$raw) cat("\nModel fit to raw moment matrix.\n")	
	if (robust && !is.null(object$robust.vcov)){
		cat("\nSatorra-Bentler Corrected Fit Statistics:\n")
		cat("\n Corrected Model Chisquare = ", chisq.adjusted, "  Df = ", df, 
				"Pr(>Chisq) =", if (df > 0) pchisq(chisq.adjusted, df, lower.tail=FALSE)
						else NA)
		if (!object$raw) {		
			cat("\n Corrected Chisquare (null model) = ", chisqNull.adjusted,  "  Df = ", dfNull)
		}
		if (df > 0 && !object$raw){
			if (!is.na(NFI.adjusted)) cat("\n Corrected Bentler-Bonett NFI = ", NFI.adjusted)
			if (!is.na(NNFI.adjusted)) cat("\n Corrected Tucker-Lewis NNFI = ", NNFI.adjusted)
			if (!is.na(CFI.adjusted)) cat("\n Corrected Bentler CFI = ", CFI.adjusted)
            if (!is.na(RNI.adjusted)) cat("\n Corrected Bentler RNI = ", RNI.adjusted)
            if (!is.na(IFI.adjusted)) cat("\n Corrected Bollen IFI = ", IFI.adjusted)
		}
		cat("\n\nUncorrected Fit Statistics:\n")
	}
	if (inherits(object, "msemObjectiveGLS") || inherits(object, "msemObjectiveFIML")) {
		AIC <- AICc <- BIC <- NA
	}
	else {
		AIC <- if ("AIC" %in% fit.indices) AIC(object) else NA
		AICc <- if ("AICc" %in% fit.indices) AICc(object) else NA
		BIC <- if ("BIC" %in% fit.indices) BIC(object) else NA
	}
	cat("\n Model Chisquare =", chisq, " Df =", df, " Pr(>Chisq) =", pchisq(chisq, df, lower.tail=FALSE))
	if (!is.na(chisqNull)) cat("\n Chisquare (null model) =", chisqNull, " Df =", dfNull)
	if (!is.na(GFI)) cat("\n Goodness-of-fit index =", GFI)
	if (!is.na(AGFI)) cat("\n Adjusted goodness-of-fit index =", AGFI)
	if (!is.na(RMSEA)) cat("\n RMSEA index = ", RMSEA, " ", 100*conf.level, "% CI: (", RMSEA.L, ", ", RMSEA.U, ")", sep="")
	if (!is.na(NFI)) cat("\n Bentler-Bonett NFI =", NFI)
	if (!is.na(NNFI)) cat("\n Tucker-Lewis NNFI =", NNFI)
	if (!is.na(CFI)) cat("\n Bentler CFI =", CFI)
    if (!is.na(RNI)) cat("\n Bentler RNI = ", RNI)
    if (!is.na(IFI)) cat("\n Bollen IFI = ", IFI)
	if (!is.na(SRMR)) cat("\n SRMR =", SRMR)
	if (!is.na(AIC)) cat("\n AIC =", AIC)
	if (!is.na(AICc)) cat("\n AICc =", AICc)
	if (!is.na(BIC)) cat("\n BIC =", BIC)
	cat("\n\n")
	if (is.null(object$initial.iterations)) cat("Iterations:", object$iterations, "\n\n")
	else cat("Iterations: initial fits,", object$initial.iterations, "  final fit,", object$iterations, "\n\n")
	for (g in 1:G){
		cat("\n  ", object$group, ": ", groups[g], "\n", sep="")
		print(group.summaries[[g]])
	}
	invisible(object)
}

summary.msemObjectiveGLS <- function(object, digits=getOption("digits"), conf.level=.90,
    fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR"), ...){
    fit.indices <- if (missing(fit.indices)){
       getOption("fit.indices")
    }
    else match.arg(fit.indices, several.ok=TRUE)
    summary.msemObjectiveML(object, digits=digits, conf.level=conf.level, robust=FALSE, fit.indices=fit.indices, ...)
	invisible(object)
}

summary.msemObjectiveFIML <- function(object, digits=getOption("digits"), conf.level=.90, ...){
	summary.msemObjectiveML(object, digits=digits, conf.level=conf.level, 
			robust=FALSE, ...)
	invisible(object)
}

deviance.msemObjectiveML <- function(object, ...) 
	object$criterion * sum(object$N - (!object$raw))

AIC.msemObjectiveML <- function(object, ..., k) {
	deviance(object) + 2*object$t
}

AICc.msemObjectiveML <- function(object, ...) {
	deviance(object) + 2*object$t*(object$t + 1)/(sum(object$N) - object$t - 1)
}

BIC.msemObjectiveML <- function(object, ...) {
	# deviance(object) + object$t*log(sum(object$N))
    deviance(object) - df.residual(object)*log(sum(object$N))
}

residuals.msem <- function(object, ...){
	result <- mapply(function(S, C) S - C, S=object$S, C=object$C, SIMPLIFY=FALSE)
	result <- lapply(result, function(res) {attr(res, "N") <- NULL; res})
	result
}

coef.msem <- function(object, ...) object$coeff

df.residual.msem <- function (object, ...){
	n <- object$n
	n.fix <- object$n.fix
	sum(n*(n + 1)/2) - object$t - sum(n.fix*(n.fix + 1)/2)
}

anova.msemObjectiveML <- function(object, model.2, ...) {
	dev.1 <- deviance(object)
	df.1 <- df.residual(object)
	dev.2 <- deviance(model.2)
	df.2 <- df.residual(model.2)
	name.1 <- deparse(substitute(object))
	name.2 <- deparse(substitute(model.2))
	df <- abs(df.1 - df.2)
	if (df == 0) 
		stop("the models have the same Df")
	if (all(object$N != model.2$N))
		stop("the models are fit to different numbers of observations")
	if (any(sapply(object$S, nrow) != sapply(model.2$S, nrow)) || !all.equal(object$S, model.2$S))
		stop("the models are fit to different moment matrices")
	chisq <- abs(dev.1 - dev.2)
	table <- data.frame(c(df.1, df.2), c(dev.1, dev.2), 
			c(NA, df), c(NA, chisq), c(NA, pchisq(chisq, df, lower.tail = FALSE)))
	names(table) <- c("Model Df", "Model Chisq", "Df", "LR Chisq", "Pr(>Chisq)")
	rownames(table) <- c(name.1, name.2)
	structure(table, heading = c("LR Test for Difference Between Models", ""), class = c("anova", "data.frame"))
}

logLik.msemObjectiveML <- function(object, ...) -0.5*deviance(object)

effects.msem <- function(object, ...){
	eff <- function(A, m, semmod){
		I <- diag(m)
		endog <- classifyVariables(semmod)$endogenous
		AA <- -A
		diag(AA) <- 1
		Total <- solve(AA) - I
		Indirect <- Total - A
		result <- list(Total = Total[endog, ], Direct = A[endog, ], Indirect = Indirect[endog, ])
		class(result) <- "semeffects"
		result
	}
	G <- length(object$groups)
	A <- object$A
	m <- object$m
	semmod <- object$semmodList
	result <- vector(G, mode="list")
	for (g in 1:G){
		result[[g]] <- eff(A[[g]], m[g], semmod[[g]])
	}
	names(result) <- object$groups
	class(result) <- "semeffectsList"
	result
}

print.semeffectsList <- function (x, digits = getOption("digits"), ...) {
	groups <- names(x)
	for (group in groups){
		cat("\n\n Group: ", group, "\n")
		print(x[[group]], digits=digits)
	}
	invisible(x)
}

standardizedCoefficients.msem <- function (object, ...){
	groups <- object$groups
	G <- length(groups)
	param.names <- object$param.names
	ram <- object$ram
	A <- object$A
	P <- object$P
	par <- coef(object)
	for (g in 1:G){
		par.names <- param.names[ram[[g]][, 4]]
		par.gr <- par[par.names]
		t <- length(par.gr)
		par.posn <- ram[[g]][, 4] != 0
		ram[[g]][par.posn, 4] <- 1:t
		group <- list(coeff=par.gr, t=t, ram=ram[[g]], A=A[[g]], P=P[[g]], par.posn=par.posn, param.names=par.names)
		class(group) <- "sem"
		cat("\n\n Group: ", groups[g], "\n")
		print(standardizedCoefficients(group, ...))
	}
}

standardizedResiduals.msem <- function (object, ...) {
	res <- residuals(object)
	S <- object$S
	for (g in 1:length(S)){
		s <- diag(S[[g]])
		res[[g]] <- res[[g]]/sqrt(outer(s, s))
	}
	res
}

normalizedResiduals.msemObjectiveML <- function (object, ...) {
	res <- residuals(object)
	N <- object$N - (!object$raw)
	C <- object$C
	for (g in 1:length(res)){
		c <- diag(C[[g]])
		res[[g]] <- res[[g]]/sqrt((outer(c, c) + C[[g]]^2)/N[g])
	}
	res
}

fscores.msem <- function (model, data = model$data, center = TRUE, scale = FALSE, ...) {
	m <- model$m
	P <- model$P
	A <- model$A
	var.names <- model$var.names
	C <- model$C
	group <- model$group
	groups <- model$groups
	G <- length(groups)
	scores <- B <- vector(G, mode = "list")
	names(scores) <- names(B) <- groups
	for (g in 1:G) {
		observed <- var.names[[g]] %in% rownames(C[[g]])
		if (all(observed)) {
			warning("there are no latent variables in group ", 
					groups[g])
		}
		IAinv <- solve(diag(m[g]) - A[[g]])
		Sigma <- IAinv %*% P[[g]] %*% t(IAinv)
		B[[g]] <- solve(Sigma[observed, observed]) %*% Sigma[observed, 
				!observed]
		rownames(B[[g]]) <- var.names[[g]][observed]
		colnames(B[[g]]) <- var.names[[g]][!observed]
		if (!is.null(data)) {
			X <- data[[g]][, var.names[[g]][observed]]
			if (center || scale) 
				X <- scale(X, center = center, scale = scale)
			scores[[g]] <- X %*% B[[g]]
		}
	}
	if (is.null(data)) 
		return(B)
	else return(scores)
}

vcov.msem <- function (object, robust=FALSE, analytic = inherits(object, "msemObjectiveML") && object$t <= 500, ...) { 
	if(robust){
		if (is.null(object$robust.vcov)) stop("robust coefficient covariance matrix not available")
		return(object$robust.vcov)
	}
	if (!analytic) 
		return(object$vcov)
	if (!inherits(object, "msemObjectiveML")) 
		stop("analytic coefficient covariance matrix unavailable")
	hessian <- function(model) {
		A <- model$A
		P <- model$P
		S <- model$S
		C <- model$C
		J <- model$J
		m <- model$m
		N <- model$N
		rams <- model$ram
		groups <- model$groups
		G <- length(groups)
		nms <- names(coef(model))
		raw <- model$raw
		Hessian <- matrix(0, nrow=length(nms), ncol=length(nms))
		rownames(Hessian) <- colnames(Hessian) <- nms
		wts <- (N - !raw)/(sum(N) - G*!raw)
		for (g in 1:G){
			I.Ainv <- solve(diag(m[g]) - A[[g]])
			Cinv <- solve(C[[g]])
			AA <- t(I.Ainv) %*% t(J[[g]])
			BB <- J[[g]] %*% I.Ainv %*% P[[g]] %*% t(I.Ainv)
			CC <- t(I.Ainv) %*% t(J[[g]])
			DD <- J[[g]] %*% I.Ainv
			dF.dBdB <- accumulate(AA %*% Cinv %*% t(AA), t(BB) %*% 
							Cinv %*% BB, AA %*% Cinv %*% BB, t(BB) %*% Cinv %*% 
							t(AA), m[g])
			dF.dPdP <- accumulate(CC %*% Cinv %*% t(CC), t(DD) %*% 
							Cinv %*% DD, CC %*% Cinv %*% DD, t(DD) %*% Cinv %*% 
							t(CC), m[g])
			dF.dBdP <- accumulate(AA %*% Cinv %*% t(CC), t(BB) %*% 
							Cinv %*% DD, AA %*% Cinv %*% DD, t(BB) %*% Cinv %*% 
							t(CC), m[g])
			ram <- rams[[g]]
			fixed <- ram[, 4] == 0
			sel.free <- ram[, 4]
			sel.free[fixed] <- 0
			one.head <- ram[, 1] == 1
			one.free <- which((!fixed) & one.head)
			two.free <- which((!fixed) & (!one.head))
			two.free.cov <- which((!fixed) & (!one.head) & (ram[, 2] != ram[, 3]))
			arrows.1 <- ram[one.head, c(2, 3), drop = FALSE]
			arrows.2 <- ram[!one.head, c(2, 3), drop = FALSE]
			arrows.2t <- ram[!one.head, c(3, 2), drop = FALSE]
			arrows.1.free <- ram[one.free, c(2, 3), drop = FALSE]
			arrows.2.free <- ram[two.free, c(2, 3), drop = FALSE]
			sel.free.1 <- sel.free[one.free]
			sel.free.2 <- sel.free[two.free]
			unique.free.1 <- unique(sel.free.1)
			unique.free.2 <- unique(sel.free.2)
			posn.matrix <- matrix(1:(m[g]^2), m[g], m[g])
			posn.free <- c(posn.matrix[arrows.1.free], (m[g]^2) + posn.matrix[arrows.2.free])
			DBB <- dF.dBdB[posn.matrix[arrows.1.free], posn.matrix[arrows.1.free], 
					drop = FALSE]
			DPP <- dF.dPdP[posn.matrix[arrows.2.free], posn.matrix[arrows.2.free], 
					drop = FALSE]
			DBP <- dF.dBdP[posn.matrix[arrows.1.free], posn.matrix[arrows.2.free], 
					drop = FALSE]
			hessian <- rbind(cbind(DBB, DBP), cbind(t(DBP), DPP))
			n1 <- length(one.free)
			n2 <- length(two.free)
			nn <- rep(c(sqrt(2), sqrt(2)/2), c(n1, n2))
			nn[c(one.free, two.free) %in% two.free.cov] <- sqrt(2)
			hessian <- hessian * outer(nn, nn)
			pars <- ram[, 4][!fixed]
			all.pars <- ram[, 4]
			t <- length(pars)
			Z <- outer(sort(unique(pars)), pars, function(x, y) as.numeric(x == y))
			hessian <- Z %*% hessian %*% t(Z)
			par.names <- c(nms[all.pars[one.free]], nms[all.pars[two.free]])
			rownames(hessian) <- colnames(hessian) <- par.names
			Hessian[par.names, par.names] <- Hessian[par.names, par.names] + wts[g]*hessian
		}
		Hessian
	}
	h <- hessian(object)
	t <- object$t
	N <- sum(object$N)
	raw <- object$raw
	G <- length(object$groups)
	param.names <- rownames(h)
	vcov <- matrix(NA, t, t)
	qr.hess <- try(qr(h), silent = TRUE)
	if (class(qr.hess) == "try-error") {
		warning("Could not compute QR decomposition of Hessian.\nOptimization probably did not converge.\n")
	}
	else if (qr.hess$rank < t) {
		warning(" singular Hessian: model is probably underidentified.\n")
		which.aliased <- qr.hess$pivot[-(1:qr.hess$rank)]
		attr(vcov, "aliased") <- param.names[which.aliased]
	}
	else {
		vcov <- (2/(N - G*(!raw))) * solve(h)
		if (any(diag(vcov) < 0)) {
			attr(vcov, "aliased") <- param.names[diag(vcov) <  0]
			warning("Negative parameter variances.\nModel may be underidentified.\n")
		}
	}
	vcov
}

robustVcovMsem <- function(model){
	G <- length(model$groups)
	t <- model$t  
	N <- model$N
	raw <- model$raw
	wts <- (N - !raw)/(sum(N) - !raw*G)
	hessian <- matrix(0, t, t)
	rownames(hessian) <- colnames(hessian) <- model$param.names
	chisq <- 0
	adj.objects <- vector(G, mode="list")
	for (g in 1:G){
		ram <- model$ram[[g]]
		parameters <- ram[, 4]
		unique.pars <- unique(parameters[parameters != 0])
		par.posn <- sapply(unique.pars, function(x) which(x == parameters)[1])
		unique.posn <- which(parameters %in% unique.pars)
		rownames(ram)[unique.posn] <- unique(model$param.names[ram[, 4]])
		ram[unique.posn, 4] <- unlist(apply(outer(unique.pars, parameters, "=="), 2, which))
		mod.g <- list(var.names=model$var.names[[g]], ram=ram,J=model$J[[g]], n.fix=model$n.fix, n=model$n[[g]], 
				N=model$N[g], m=model$m[[g]], t=length(unique.pars),
				coeff=model$coeff[parameters],  criterion=model$group.criteria[g],  S=model$S[[g]], raw=model$raw,
				C=model$C[[g]], A=model$A[[g]], P=model$P[[g]])
		adj.objects[[g]] <- sbchisq(sem.obj=mod.g, sem.data=model$data[[g]])
		hess <- solve((N[g] - 1)*robustVcov(mod.g, adj.obj=adj.objects[[g]]))
		pars <- rownames(hess)
		hessian[pars, pars] <- hessian[pars, pars] + wts[g]*hess
		chisq <- chisq + adj.objects$chisq.scaled
	}
	return(list(vcov=solve(hessian)/(sum(N) - !raw*G), chisq.scaled=chisq, adj.objects=adj.objects))
}
