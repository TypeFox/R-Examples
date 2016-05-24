#' @title Defining optimum partition for a specific variable (regression model)
#' @details
#' Internal function. \code{splitopt.reg} is called by \code{all.part.reg}.
#' @param x matrix or data.frame with the data.
#' @param splits vector indicating the binary partition
#' @param fact vector indicating the variable
#' @param method string indicating the method: LM or LAD
#' @param \dots Further arguments passed on to \code{\link{splitopt.reg}}.
#' @return list containing information of the optimum partition for a specific variable
#' @keywords internal
#' @export

splitopt.reg <- function(x,splits,fact,method,...)
 {
	Fi				=	NULL                                             	
	pval.split		=	NULL                                     
	df.num			=	NULL
	df.den			=	NULL
	ind1			=	NULL
	ind2			=	NULL
	new.mod			=	list()
	modtwo.list		=	list()
	particion.list	=	list()
	q				= 	ncol(x)                                          	
	n				= 	nrow(x)                                           	

	g.resp	=	as.matrix(x[,1])
	g.pred	=	cbind(rep(1,nrow(as.matrix(x[, -1]))),as.matrix(x[,-1]))
	
	if(method=="lm")
	{
		reg0	=	lm(g.resp~g.pred-1)                   	
		SSR0	=	sum(reg0$residuals^2)
		df0		=	(nrow(g.pred) - ncol(g.pred))
	}
	if(method=="lad")
	{
		reg0	=	rq(g.resp~g.pred-1,method="fn")                   	
		SSR0	=	sum(abs(reg0$residuals))    
		df0		=	(nrow(g.pred) - ncol(g.pred))
	}
	for (i in 1:nrow(splits)) 
	{                    
		modnum	=	as.numeric(fact)                        
		split	=	splits[i,]
		particion.list[[length(particion.list)+1]]	=	split
		modtwo = split[modnum]                             
		
		bin.lev	=	bin.levels	(fact,split)
		new.mod[[length(new.mod)+1]]	=	bin.lev
		modtwo.list[[length(modtwo.list)+1]]	=	modtwo
		
		x1 = subset(x,modtwo==1)                        
		x2 = subset(x,modtwo==2)
	
		if(nrow(x1)<=ncol(x1)||nrow(x2)<=ncol(x2)) next
		
		if(nrow(x1)<=15||nrow(x2)<=15) next
		
		ind1[i]	=	nrow(x1)
		ind2[i]	=	nrow(x2)
				
		l1.resp	=	as.matrix(x1[,1])
		l1.pred	=	cbind(rep(1,nrow(as.matrix(x1[, -1]))),as.matrix(x1[,-1]))

		l2.resp	=	as.matrix(x2[,1])
		l2.pred	=	cbind(rep(1,nrow(as.matrix(x2[,-1]))),as.matrix(x2[,-1]))
		
		df1		=	(nrow(l1.pred)+nrow(l2.pred))-(ncol(l1.pred)+ncol(l2.pred))

		if(method=="lm")
		{
			reg11		=	lm(l1.resp~l1.pred-1)        	
			SSR11		=	sum(reg11$residuals^2)   
			reg22		=	lm(l2.resp~l2.pred-1)      	
			SSR22		=	sum(reg22$residuals^2)    
		}
		if(method=="lad")
		{
			reg11		=	rq(l1.resp~l1.pred-1,method="fn")        	
			SSR11		=	sum(abs(reg11$residuals))  
			reg22		=	rq(l2.resp~l2.pred-1,method="fn")      	
			SSR22		=	sum(abs(reg22$residuals))    
		}
		
		SSR1			=	SSR11+SSR22
                                                                                       
		Fi[i]			=	((SSR0-SSR1)/(df0-df1))/(SSR1/df1)         
		pval.split[i]	=	pf(Fi[i],(df0-df1),df1,lower.tail=FALSE)  
		df.num[i]		=	df0-df1
		df.den[i]		=	df1
		} 
	
	fun.min 		= 	f.min(pval.split)
	pval.opt		=	fun.min$v.min                           	
	pos.opt			=	fun.min$p.min               
	F.opt			=	Fi[pos.opt]                                           	                        		
	mod.opt			=	unlist(new.mod[pos.opt])
	df.num.opt		=	df.num[pos.opt]
	df.den.opt		= 	df.den[pos.opt] 
	ind1.opt		=	ind1[pos.opt] 
	ind2.opt		=	ind2[pos.opt]
	modtwo.opt		=	unlist(modtwo.list[pos.opt])
	particion.opt	=	unlist(particion.list[pos.opt])
	list(Fi=F.opt,pval=pval.opt,mod=mod.opt,df.num=df.num.opt,df.den=df.den.opt,ind1=ind1.opt,ind2=ind2.opt,modtwo=modtwo.opt,particion=particion.opt)
}
