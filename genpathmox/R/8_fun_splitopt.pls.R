#' @title Defining optimum partition for a specific variable.
#' @details
#' Internal function. \code{splitopt.pls} is called by \code{all.part.pls}.
#' @param x matrix or data frame containing the data.
#' @param inner A square (lower triangular) boolean matrix representing the inner 
#' model (i.e. the path relationships between latent variables).
#' @param outer list of vectors with column indices or column names from Data indicating 
#' the sets of manifest variables forming each block (i.e. which manifest variables correspond to each block).
#' @param scaling optional list of string vectors indicating the type of 
#' measurement scale for each manifest variable specified in \code{blocks}.
#' \code{scaling} must be specified when working with non-metric variables.
#' Possible values: \code{"num"} (numeric), \code{"raw"}, \code{"nom"} (nominal), 
#' and \code{"ord"} (ordinal).
#' @param mode character vector indicating the type of measurement for each
#' block. Possible values are: \code{"A", "B", "newA", "PLScore", "PLScow"}. 
#' The length of \code{mode} must be equal to the length of \code{outer}.
#' @param scheme string indicating the type of inner weighting
#' scheme. Possible values are \code{"centroid"}, \code{"factorial"}, or
#' \code{"path"}.
#' @param scaled whether manifest variables should be standardized. 
#' Only used when \code{scaling = NULL}. When (\code{TRUE}, data is 
#' scaled to standardized values (mean=0 and variance=1). 
#' @param splits vector indicating the binary partition
#' @param fact vector indicating the variable
#' @param method string indicating the method: LM or LAD
#' @param n.node number indicating a stop condition
#' @param \dots Further arguments passed on to \code{\link{splitopt.pls}}. 
#' @return list containing information of the optimum partition for a specific variable
#' @keywords internal
#' @export

splitopt.pls <- function(x,inner, outer,mode,scheme,scaling,scaled,splits,fact,method,n.node,...) 
{	
	Fi				=	NULL                                             	
	pval.split		=	NULL                                     
	df.num			=	NULL
	df.den			=	NULL
	g1.ind			=	NULL
	g2.ind			=	NULL
	new.mod			=	list()
	modtwo.list		=	list()
	partition.list	=	list()
	
	#Ho

	pls	=	plspm(x,inner,outer,mode,scaling,scheme,scaled=scaled)
	
	LV = pls$scores

	global	=	build.block(inner,latent=LV)
	
	n.block	=	ncol(inner)+1
	
	g.resp	=	global$resp
	g.pred	=	global$x.block
	
	if (method=="lm")
	{
		reg0	=	lm(g.resp~g.pred-1)                   	
		SSR0	=	sum(reg0$residuals^2)
		df0	=	(nrow(g.pred) - ncol(g.pred))
	}
	if (method=="lad")
	{
		reg0	=	rq(g.resp~g.pred-1,method="fn")                   	
		SSR0	=	sum(abs(reg0$residuals))    
		df0	=	(nrow(g.pred) - ncol(g.pred))
	}
	
	for (i in 1:nrow(splits)) 
	{                 
		modnum									=	as.numeric(fact)                        
		split									=	splits[i,]
		partition.list[[length(partition.list)+1]]	=	split
		modtwo 									= 	split[modnum]                             
		
		bin.lev									=	bin.levels(fact,split)
		new.mod[[length(new.mod)+1]]				=	bin.lev
		modtwo.list[[length(modtwo.list)+1]]			=	modtwo
		
		g1.latent 								= 	subset(LV,modtwo==1)                    
		g2.latent 								= 	subset(LV,modtwo==2)
		
		g1.ind[i]									=	nrow(g1.latent)
		g2.ind[i]									=	nrow(g2.latent)
		
		if (g1.ind[i] < n.block || g2.ind[i] < n.block) next
		
		if (g1.ind[i] <= n.node || g2.ind[i] <= n.node) next
	
		g1											=	build.block(inner,latent=g1.latent)
		g2											=	build.block(inner,latent=g2.latent)
		
		g1.resp										=	g1$resp
		g1.pred										=	g1$x.block

		g2.resp										=	g2$resp
		g2.pred										=	g2$x.block
		
		if (nrow(g1.pred) <= ncol(g1.pred) || nrow(g2.pred) <= ncol(g2.pred)) next
		
		df1	= (nrow(g1.pred)+nrow(g2.pred))-(ncol(g1.pred)+ncol(g2.pred))

		if(method == "lm")
		{
			reg11		=	lm(g1.resp~g1.pred-1)        	
			SSR11		=	sum(reg11$residuals^2)   
			reg22		=	lm(g2.resp~g2.pred-1)      	
			SSR22		=	sum(reg22$residuals^2)    
		}
		if(method == "lad")
		{
			reg11		=	rq(g1.resp~g1.pred-1,method="fn")        	
			SSR11		=	sum(abs(reg11$residuals))  
			reg22		=	rq(g2.resp~g2.pred-1,method="fn")      	
			SSR22		=	sum(abs(reg22$residuals))    
		}
		
		SSR1		=	SSR11+SSR22
                                                                                       
		Fi[i]		=	((SSR0-SSR1)/(df0-df1))/(SSR1/df1)         
		pval.split	=	pf(Fi,(df0-df1),df1,lower.tail=FALSE)  
		df.num[i]		=	df0-df1
		df.den[i]		=	df1
		}
	
	fun.min 		= f.min(pval.split)
	
	pval.opt			=	fun.min$v.min                           	
	pos.opt			=	fun.min$p.min               
	F.opt			=	Fi[pos.opt]                                           	                        		
	mod.opt			=	unlist(new.mod[pos.opt])
	df.num.opt		=	df.num[pos.opt]
	df.den.opt		= 	df.den[pos.opt] 
	ind1.opt			=	g1.ind[pos.opt] 
	ind2.opt			=	g2.ind[pos.opt]
	modtwo.opt		=	unlist(modtwo.list[pos.opt])
	partition.opt		=	unlist(partition.list[pos.opt])
	
	list(Fi=F.opt,
		pval=pval.opt,
		mod=mod.opt,
		df.num=df.num.opt,
		df.den=df.den.opt,
		g1.ind=ind1.opt,
		g2.ind=ind2.opt,
		modtwo=modtwo.opt,
		partition=partition.opt)
		
}
