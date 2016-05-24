#===================================================================================
# ----------------------------------------------------------------------------------
# ---------- Estimating adjusted prevalence ratio and interval confidence ----------
# ----------------------------------------------------------------------------------
#
# The prcorrlog.R script are a set functions created by 
#
# Raydonal Ospina M. 
# Department of Statistics 
# Federal University of Pernambuco, Recife - Brazil
# contact: rayospina@gmail.com, raydonal@de.ufpe.br
# 
# Leila D. Amorim 
# Department of Statistics 
# Federal University of Bahia, Salvador - Brazil
# contact: leiladen@ufba.br
#
# Originaly created 29  march  2010. 
# Last modification: 27 may 2011 at 12:55 Version 1.0
# Last modification: 25 Oct 2011 at 09:10 Version 1.1
# Last modification: 19 Set 2013 at 11:10 Version 1.2
#
# We used some functions to lmer package. In this case, copyright (C)
# Douglas Bates <bates@stat.wisc.edu> and 
# Martin Maechler <maechler@R-project.org> All rights reserved.
#
#====================================================================================


#=====================================================================================
#                 print.cin ( is a method to print IC using classes)
#=====================================================================================
"print.cin" <-
function (x, digits = 5, ...) #max(4, getOption("digits")-3)
{
    cat("\n")
    cat(x$head, "\n")
    rq <- structure(x$ci, names = x$ends)
    print(rq, digits = digits)
    cat("\n")
    invisible(x)
}

#=====================================================================================



#=====================================================================================
#   Return the pairs of expressions that separated by vertical bars. ( lmer package )
#=====================================================================================
"findbars" <- function(term)
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}
#=====================================================================================


#=====================================================================================
#  Return the list of '/'-separated terms in an expression that contains slashes
#  ( lmer package )
#=====================================================================================
"slashTerms" <- function(x)
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}
#=====================================================================================


#=====================================================================================
# from a list of length 2 return recursive interaction terms
#   ( lmer package )
#=====================================================================================
"makeInteraction" <- function(x)
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}
#=====================================================================================


#=====================================================================================
# expand any slashes in the grouping factors returned by findbars
#   ( lmer package )
#=====================================================================================
"expandSlash" <- function(bb)
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}
#=====================================================================================


#=====================================================================================
#                 is.glm (is a glm object?)
#=====================================================================================
"is.glm" <- function (x) 
inherits(x, "glm")
#=====================================================================================


#=====================================================================================
#                 is.glmer ( is a glmer object?)
#=====================================================================================
"is.glmer" <- function (x) 
inherits(x, "glmerMod")
#=====================================================================================


#=====================================================================================
# pr.conditional() Estimate prevalence ratios and their confidence intervals using
# logistic models and conditional procedure		
#=====================================================================================
"pr.conditional" <- function(object, conf = conf,...)
{
	if(is.glm(object))
	{
		fit.coef = object$coefficients 
		fit.cov = as.matrix(vcov(object))
		m.aux = cbind(1, diag(rep(1,length(fit.coef)-1)))
		
		pr = (1+exp(-fit.coef[1])) /  (1+exp(-(m.aux%*%fit.coef)))
		
		p1.hat = exp(fit.coef[1])/(1+exp(fit.coef[1]))		
		p0.hat = exp((m.aux%*%fit.coef)) / (1+exp((m.aux%*%fit.coef)))

		m0.aux= cbind(1, diag(rep(0,length(fit.coef)-1)))

		xstart = as.matrix((1-p1.hat)*m0.aux-(diag(as.vector(1-p0.hat)) %*%  m.aux))
		
		
		lnrp = log(pr)
		varlogrp = diag((xstart) %*%  fit.cov %*% t(xstart))
		sdlnrp = sqrt(varlogrp)
		quant.cond = abs(qnorm((1-conf)/2)) 

		liminfIC.cond = exp(lnrp -  quant.cond*sdlnrp )
		limsupIC.cond = exp(lnrp +  quant.cond*sdlnrp)
		
		 conditional.rp = cbind(pr, liminfIC.cond , limsupIC.cond )



		nam = names(fit.coef[-1])

		dimnames(conditional.rp) = list(NULL, c("Estimate", paste(as.character(c((1 - conf)/2, 
        1 - ((1 - conf)/2)) * 100), "%", sep = "")))
		rownames(conditional.rp) = nam

		conditional.rp		

	}
	else
	{
		if(is.glmer(object))
		{
		fit.coef= object@beta
		fit.cov = as.matrix(vcov(object))

		m.aux = cbind(1, diag(rep(1,length(fit.coef)-1)))
		pr = (1+exp(-fit.coef[1])) /  (1+exp(-(m.aux%*%fit.coef)))

		p1.hat = exp(fit.coef[1])/(1+exp(fit.coef[1]))		
		p0.hat = exp((m.aux%*%fit.coef)) / (1+exp((m.aux%*%fit.coef)))

		m0.aux= cbind(1, diag(rep(0,length(fit.coef)-1)))

		xstart = as.matrix((1-p1.hat)*m0.aux-(diag(as.vector(1-p0.hat)) %*%  m.aux))
		lnrp = log(pr)
		varlogrp = diag((xstart) %*%  fit.cov %*% t(xstart))
		sdlnrp = sqrt(varlogrp)
		quant.cond = abs(qnorm((1-conf)/2)) 

		liminfIC.cond = exp(lnrp -  quant.cond*sdlnrp )
		liminfIC.cond = ifelse(liminfIC.cond< 0, 0, liminfIC.cond) ## Truncate the IC for negative values

		limsupIC.cond = exp(lnrp +  quant.cond*sdlnrp)
		limsupIC.cond = ifelse(limsupIC.cond< 0, 0, limsupIC.cond) ## Truncate the IC for negative values
		
		 conditional.rp = cbind(pr, liminfIC.cond , limsupIC.cond )

		

		nam = names(fit.coef[-1])
		dimnames(conditional.rp) = list(NULL, c("Estimate", paste(as.character(c((1 - conf)/2, 
        1 - ((1 - conf)/2)) * 100), "%", sep = "")))
		rownames(conditional.rp) = nam
		conditional.rp

		}
		else
		{
		stop("There are not glm and glmer models. Revise arguments and try again")
		}
	}
conditional.rp	
}
#=====================================================================================


#=====================================================================================
# pr.marginal() Estimate prevalence ratios and their confidence intervals using
# logistic models and marginal procedure		

#=====================================================================================
"pr.marginal" <- function(object, conf = conf,...) 
{
	if(is.glm(object))
	{
	#glm class model
	coef.glm =  object$coefficients
	fit.coef = as.vector(t(coef.glm))
	fit.cov = as.matrix(vcov(object))
	m.aux = as.matrix(object$model[-1])

	
	m.covariate = as.matrix(cbind(1,m.aux))
	n.coef = dim(m.covariate)[2]
	p1.marg.aux = colMeans(m.aux)
	p0.marg.aux = 1- p1.marg.aux #mean(m.aux)

	marginal.rp = NULL


	m.ind0 = cbind(1, diag(1,(n.coef-1)))
	m.ind1 = cbind(1, diag(0,(n.coef-1)))

	for(i in 1:(n.coef-1))
	{
		m.aux.p0 =  m.covariate
		m.aux.p0[,i+1] = 0
		p0.marg = c(1,p0.marg.aux)
		p0.marg[i+1]=0

		p0.aux = (exp(m.aux.p0%*%fit.coef))/(1+exp(m.aux.p0%*%fit.coef))
		p0.med = mean(p0.aux)


		m.aux.p1 =  m.covariate
		m.aux.p1[,i+1] = 1
		p1.marg = c(1,p1.marg.aux)
		p1.marg[i+1]=0
	
		p1.aux = (exp(m.aux.p1%*%fit.coef))/(1+exp(m.aux.p1%*%fit.coef))
		p1.med = mean(p1.aux)

		rp.marg =  p1.med / p0.med #prevalence ratio

		

		# IC by delta method 
		lnrp.marg = log(rp.marg) 
		# xstar.marg = (1-p1.med)*p1.marg - (1-p0.med)*p0.marg


		xstar.marg = ((1-p1.med) * m.ind0[i,]) - ((1-p0.med)*m.ind1[i,])

		varlogrp.marg = as.vector(xstar.marg) %*% as.matrix(fit.cov) %*% as.vector(t(xstar.marg))

		stdlogrp.marg = sqrt(varlogrp.marg)
		quant.marg = abs(qnorm((1-conf)/2)) 

		liminfIC.marg = exp(lnrp.marg -  quant.marg*stdlogrp.marg )
		liminfIC.marg = ifelse(liminfIC.marg< 0, 0, liminfIC.marg) ## Truncate the IC for negative values

		limsupIC.marg = exp(lnrp.marg +  quant.marg*stdlogrp.marg )

		limsupIC.marg = ifelse(limsupIC.marg< 0, 0, limsupIC.marg) ## Truncate the IC for negative values

		marginal.rp =(cbind(marginal.rp, c(rp.marg, liminfIC.marg, limsupIC.marg) ) ) # pr whit IC

	}

		
		marginal.rp = t(marginal.rp)

		nam = names(coef.glm[-1])

	dimnames(marginal.rp) = list(NULL, c("Estimate", paste(as.character(c((1 - conf)/2, 
        1 - ((1 - conf)/2)) * 100), "%", sep = "")))
		rownames(marginal.rp) = nam
		marginal.rp # output relative prevalence whit IC by delta method
	}
	else
	{
	# glmer class model
	coef.lmer =  object@beta
	fit.coef = as.vector(t(coef.lmer))
	fit.cov = as.matrix(vcov(object))
	m.model = object@frame
	m.aux = as.matrix(m.model[c(-1, -dim(m.model)[2])])
	m.covariate = as.matrix(cbind(1,m.aux))
	n.coef = dim(m.covariate)[2]
	p1.marg.aux = colMeans(m.aux)
	p0.marg.aux = 1- p1.marg.aux #mean(m.aux)

	marginal.rp = NULL

	m.ind0 = cbind(1, diag(1,(n.coef-1)))
	m.ind1 = cbind(1, diag(0,(n.coef-1)))

	for(i in 1:(n.coef-1))
	{
		m.aux.p0 =  m.covariate
		m.aux.p0[,i+1] = 0
		p0.marg = c(1,p0.marg.aux)
		p0.marg[i+1]=0



		p0.aux = (exp(m.aux.p0%*%fit.coef))/(1+exp(m.aux.p0%*%fit.coef))
		
		p0.med = mean(p0.aux)

		m.aux.p1 =  m.covariate
		m.aux.p1[,i+1] = 1
		p1.marg = c(1,p1.marg.aux)
		p1.marg[i+1]=0

		p1.aux = (exp(m.aux.p1%*%fit.coef))/(1+exp(m.aux.p1%*%fit.coef))
		p1.med = mean(p1.aux)

		rp.marg =  p1.med / p0.med #prevalence ratio



		# IC by delta method 
		lnrp.marg = log(rp.marg)

		xstar.marg = ((1-p1.med) * m.ind0[i,]) - ((1-p0.med)*m.ind1[i,])
		#xstar.marg = (1-p1.med)*p1.marg - (1-p0.med)*p0.marg

		varlogrp.marg = as.vector(xstar.marg) %*% as.matrix(fit.cov) %*% as.vector(t(xstar.marg))
		stdlogrp.marg = sqrt(varlogrp.marg)
		quant.marg = abs(qnorm((1-conf)/2))

		liminfIC.marg = exp(lnrp.marg -  quant.marg*stdlogrp.marg )
		limsupIC.marg = exp(lnrp.marg +  quant.marg*stdlogrp.marg )

		marginal.rp =(cbind(marginal.rp, c(rp.marg,liminfIC.marg, limsupIC.marg) ) ) # pr whit IC


	}


		marginal.rp = t(marginal.rp)
		nam = names(coef.lmer[-1])
		dimnames(marginal.rp) = list(NULL, c("Estimate", paste(as.character(c((1 - conf)/2, 
        1 - ((1 - conf)/2)) * 100), "%", sep = "")))
		rownames(marginal.rp) = nam

		marginal.rp # output relative prevalence whit IC by delta method


	}
}
#=====================================================================================


#=====================================================================================
# boot.conditional() Auxiliar bootstrap function to calculate confidence intervals to 
# prevalence ratio using the conditional model
#=====================================================================================
"boot.conditional" <- function(data, indices, object, conf,...){
  data<-data[indices,]
  fit = update(object, data=data)

	obj.aux = pr.conditional(fit,conf)

	as.vector(t(obj.aux[,1]))
  }
#=====================================================================================


#=====================================================================================
# boot.marg() Auxiliar bootstrap function to calculate confidence intervals to 
# prevalence ratio using the marginal model
#=====================================================================================
"boot.marg" <- function(data, indices, object, conf,...){
  data<-data[indices,]
  fit = update(object, data=data)

	obj.aux = pr.marginal(fit,conf)

	as.vector(t(obj.aux[,1]))
  }
#=====================================================================================


#=====================================================================================
# prLogisticBootMarg() Confidence intervals to prevalence ratio using the marginal model
# by default the Bootstrap IC are obtained by using the normal and percentile method 
# whit confidence nivel at 95 %
#=====================================================================================
prLogisticBootMarg <- function(
		object,		# A fitted model
		data,		# data set used by estimate the model
		conf = 0.95,	# A scalar or vector containing the confidence level(s) of the required interval(s). 
		R = 99,		# The number of bootstrap replicates
		... 		# ... Any extra arguments
		)
{
#require(boot)

	if(is.glm(object))
	{
	nam = names(object$coefficients[-1])
	}
	else
	{
	if(is.glmer(object))
	{
	nam = names(object@beta[-1])
	}
	else{
	stop("There are not glm and glmer models. Review the arguments function and try again")
	}
	}

output= boot(data = data, statistic = boot.marg, R = R, object=object, conf=conf,...)

n = length(output$t0)

ic.boot = NULL
for(i in 1:n)
{
ic.boot.aux = boot.ci(output, type=c("norm", "perc"),index=i)
ic.boot.normal = ic.boot.aux$normal[c(2,3)]
ic.boot.percent = ic.boot.aux$percent[c(4,5)]
ic.boot = cbind(ic.boot, c(ic.boot.normal, ic.boot.percent))
}

	ic.boot = ifelse(ic.boot< 0, 0, ic.boot) ## Truncate the IC for negatve values	
	ic.boot=cbind(output$t0,t(ic.boot))
	


	dimnames(ic.boot) = list(NULL, c("Estimate", paste(as.character(c(
	(1 - conf)/2, 1 - ((1 - conf)/2), (1 - conf)/2, 1 - ((1 - conf)/2)) * 100), 
	"%", sep = "")))

	rownames(ic.boot) = nam

	ic.boot # We can modified this function to obtain other statistics

	res <- list() # Output

	head1 <- paste(paste(as.character(conf * 100), "%", sep = ""),  
	c("Confidence Interval using Bootstrap Method\n"))
	head2 <- c("\t\t      Normal\t  Percentile")
	res$head <- rbind(head1, head2)
	res$ci = ic.boot

	class(res) <- "cin"
	res

}
#=====================================================================================


#=====================================================================================
# prLogisticBootCond() Confidence intervals to prevalence ratio using the conditional model
# by default the Bootstrap IC are obtained by using the normal and percentile method 
# whit confidence nivel at 95 %
#=====================================================================================
prLogisticBootCond =function(
		object,		# A fitted model
		data,		# data set used by estimate the model
		conf = 0.95,	# A scalar or vector containing the confidence level(s) of the required interval(s). 
		R = 99,		# The number of bootstrap replicates
		... 		# ... Any extra arguments
		)
{
#require(boot)

	if(is.glm(object))
	{
	nam = names(object$coefficients[-1])
	}
	else
	{
	if(is.glmer(object))
	{
	nam = names(object@beta[-1])
	}
	else{
	stop("There are not glm and glmer models. Review the arguments functions and try again")
	}
	}

output= boot(data = data, statistic = boot.conditional, R = R, object=object, conf=conf,...)

n = length(output$t0)

ic.boot = NULL
for(i in 1:n)
{
ic.boot.aux = boot.ci(output, type=c("norm", "perc"),index=i)
ic.boot.normal = ic.boot.aux$normal[c(2,3)]
ic.boot.percent = ic.boot.aux$percent[c(4,5)]
ic.boot = cbind(ic.boot, c(ic.boot.normal, ic.boot.percent))
}


	ic.boot=cbind(output$t0,t(ic.boot))
	ic.boot = ifelse(ic.boot< 0, 0, ic.boot) ## Truncate the IC for negatve values

	dimnames(ic.boot) = list(NULL, c("Estimate", paste(as.character(c(
	(1 - conf)/2, 1 - ((1 - conf)/2), (1 - conf)/2, 1 - ((1 - conf)/2)) * 100), 
	"%", sep = "")))

	rownames(ic.boot) = nam

	ic.boot # We can modified this function to obtain other statistics

	res <- list() # Output

	head1 <- paste(paste(as.character(conf * 100), "%", sep = ""),  
	c("Confidence Interval using Bootstrap Method\n"))
	head2 <- c("\t\t      Normal\t  Percentile")
	res$head <- rbind(head1, head2)
	res$ci = ic.boot

	class(res) <- "cin"
	res

}
#=====================================================================================


#=====================================================================================
# prLogisticDelta()	Estimate prevalence ratios and their confidence intervals using
#		logistic models		
#=====================================================================================

prLogisticDelta <- function (
		formula,   # A symbolic description of the model to be fitted. 
		cluster = FALSE, # A logical value indicating if the model have a random effects
		pattern = c("conditional", "marginal"), # Method for estimate prevalence ratio 
		conf = 0.95, # A scalar or vector containing the confidence level(s) of the required interval(s). 
		dataset, # data set used to estimate the model 
		... # ... Any extra arguments
		)
{

	#require(lme4)    
        
	if (!cluster) 
	{
	bars <- expandSlash(findbars(formula[[3]]))
        if ((!length(bars)==FALSE))  
		{
		stop("Random effects terms specified in formula. Use cluster = TRUE ")
		}
		else
		{
		fit = glm(formula, family=binomial, data=dataset)	
		}
	}
	else
	{
	 if (cluster)
	 bars <- expandSlash(findbars(formula[[3]]))
	 if (!length(bars)) 
	 {
	 stop("No random effects terms specified in formula. Use cluster = FALSE")	  
         }
	 else  
	 {
	 fit = glmer(formula, family = binomial, data = dataset)		
         }
	}

	fit

pattern <- match.arg(pattern)

	if((pattern=="conditional") )
	{
	out = pr.conditional(fit, conf=conf,...)
	}
	else
	{
		if((pattern=="marginal") )
		{
		out = pr.marginal(fit, conf=conf,...)
		}
		else
		{
		stop("There are some problem whit the arguments. Review and try again")
		}
	}
	



 res <- list()
    res$head <- paste(paste(as.character(conf * 100), "%", sep = ""),  
	c("Confidence Interval using Delta method"))
    res$ci <- out
    class(res) <- "cin"
    res


}
#=====================================================================================





