modSearch <- function(formula, x, y, endpoint, method, optimal, force.in, nvmax){
	source <- data.frame(x,y)
	if (NCOL(x)==1){ xnames <- as.character(formula[[3]]) 
	} else {
		xnames <- colnames(x)
		if (is.null(force.in)){
			xnamesSort <- xnames
		} else {
			xnamesSort <- xnames[c(force.in,(1:NCOL(x))[-force.in])]
		}
	}
	colnames(source) <- c(xnames, as.character(formula[[2]]))
	family <- ifelse(endpoint=="quantitative", "gaussian", "binomial")
	if (NCOL(x)>1){
		if (endpoint=="quantitative"){
			search <- summary(regsubsets(x=x, y=y, method=method, force.in=force.in, 
			nvmax=nvmax, really.big=NCOL(x)>50))
		} else {
			w <- glm(formula, family=binomial, data=source)$weights
			search <- summary(regsubsets(x=x, y=y, weights=w, method=method, force.in=force.in,
			nvmax=nvmax, really.big=NCOL(x)>50))
		}
		if (optimal!="rsq"){ opt <- which(search[[optimal]]==min(search[[optimal]]))
		} else { opt <- which(search[[optimal]]==max(search[[optimal]])) }
		idx <- search$which[opt,-1]
		names <- as.character(na.omit(ifelse(idx, xnamesSort, NA)))
		rsq <- search$rsq[opt]
		mod <- glm(as.formula(paste(formula[[2]],"~",paste(names,collapse="+"))), family=family, data=source)
	} else {
		mod <- glm(formula, family=family, data=source)
		rsq <- NULL
		names <- xnames
	}
	list(mod=mod, rsq=rsq, names=names)
}
