beta.list <- function(v) {
	p <- max(as.numeric(gsub("(beta\\[|,[[:digit:]]+\\])","",grep("beta", names(v), value=TRUE))))
	lapply(1:p, function(j) c(v[grep(paste("beta\\[", j, ",[[:digit:]]+\\]", sep=""), names(v))]))
}

gpcm.bug <- function(v,cats,mdl,gibbs=c("bugs","jags")){
  gibbs <- match.arg(gibbs)
	# All models have theta
	theta.pars <- v[grep("theta", names(v))]
	# Models with alpha
      if(mdl=="gpcm.bug" || mdl=="tpl.bug"){
	      alpha.pars <- v[grep("alpha", names(v))]		
	} else {
		alpha.pars <- rep(1,length(cats))
	}
	#Dichotomous or polytomous
	if(mdl=="tpl.bug" || mdl=="rasch.bug"){
		beta.pars <- as.numeric(v[grep("delta", names(v))])
		beta.pars <- cbind(0,beta.pars)	
	} else {
		pars <- beta.list(v)
		p <- length(pars)
		#Dump beta parameters into matrix
		beta.pars <- matrix(0,ncol=max(cats),nrow=p)
		if(gibbs=="bugs"){
			for(i in 1:p){
				ed <- cats[i]-1
				beta.pars[i,2:cats[i]] <- pars[[i]][1:ed]
			}
		} else {
			for(i in 1:p){
				ed <- cats[i]
				beta.pars[i,2:cats[i]] <- pars[[i]][2:ed]
			}
		}
	}
	return(list(beta.pars,theta.pars,alpha.pars))
}

