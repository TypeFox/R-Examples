# 
# Author: JonathanRosenblatt
###############################################################################

#----------- Validating input of compare procedure------------#



mu.test.class<- function(classes){
	if (any(classes!='Mutoss')) stop('Input is not a "Motoss" class object')
	##TODO: will this cause prolems for inherited class objects?
}

mu.test.type<- function(types){
	if( any(types != types[1]) ){
		message(' Notice:You are comparing methods for different error types. \n These should not be compared!	\n Output will be generated nevertheless. \n')
	}
}


mu.test.rates<-function(rates){
	if( any(rates != rates[1]) ){
		message(' Notice:You are comparing methods with different error rates. \n These should not be compared!	\n Output will be generated nevertheless. \n')
	}
}

mu.test.same.data<- function(pvals){
	pvals.different<- any( apply(pvals,1, function(x) any(x!=x[1])))
	if(pvals.different) stop('Different data was used for suppied procedures.')		
}


mu.test.name<- function(hyp.names){
	names.different<- any( apply(hyp.names,1, function(x) any(x!=x[1])))
	if(names.different) message('Notice: Hypotheses have different names. Can they be compared?')		
}


#-------- Create comparison list for Mutoss objects--------------#


compareMutoss<-function(...){
	objects<-list(...)
	
	classes<- sapply(objects, function(x) class(x) 	)#getting object classes
	mu.test.class(classes) #testing for compatible object classes
	
	types<-	sapply(objects, function(x) x@errorControl@type)#getting error control types
	mu.test.type(types) #testing for compatible error control types
	
	rates<-sapply(objects, function(x) x@errorControl@alpha)#extracting error rates
	
	pi.nulls<- as.numeric(lapply(objects, function(x) x@pi0))#extracting pi0 estimates
	
	pvalues<- sapply(objects, function(x) x@pValues )# getting adjusted pvals
	mu.test.same.data(pvalues)
	m<- nrow(pvalues)
	
	raw.hyp.names<- lapply(objects, function(x) x@hypNames )# getting hypothesis names
	if(all(sapply(raw.hyp.names, function(x) identical(x, character(0))))){
		hyp.names<- paste('hyp', 1:m, sep='') 
	}	
	else if(all(sapply(raw.hyp.names, function(x) length(x)==m))) {
		mu.test.name(raw.hyp.names)
		hyp.names<- raw.hyp.names[,1]
	}
	
	##TODO: [JR] Deal with mising hyp names only in a subset of objects
	
	# Preparing Raw Pvalues
	raw.pvals<- pvalues[,1]
	pval.order<- order(raw.pvals)
	pval.ranks<- rank(raw.pvals)
	raw.pvals.frame<-data.frame(
			pValue=raw.pvals,
			order=pval.order,
			ranks=pval.ranks)
	row.names(raw.pvals.frame)<-hyp.names 
	
	#Preparing adjusted pvalues
	adj.pvals<- lapply(objects, function(x) x@adjPValues )
	method.names<- unlist(lapply(adj.pvals, function(x) attributes(x)[1]))
	adj.pvals.frame<- data.frame(adj.pvals)
	colnames(adj.pvals.frame)<- method.names
	row.names(adj.pvals.frame)<- hyp.names
	
	#Preparing critical values
	critical<- lapply(objects, function(x) x@criticalValues )
	method.names<- unlist(lapply(critical, function(x) attributes(x)[1]))
	critical.frame<- data.frame(critical)
	colnames(critical.frame)<- method.names
	row.names(critical.frame)<- hyp.names	
	
	#Preparing decisions
	rejections<- lapply(objects, function(x) x@rejected )
	method.names<- unlist(lapply(rejections, function(x) attributes(x)[1]))
	rejections.frame<- data.frame(rejections)
	colnames(rejections.frame)<- method.names
	row.names(rejections.frame)<- hyp.names
	
	##TODO: [JR] Add groud truth to comparison method
	
	comparing<- list(
			types=types,
			rates=rates,
			pi.nulls=pi.nulls,
			raw.pValues=raw.pvals.frame,
			adjusted.pvals=adj.pvals.frame,
			criticalValue=critical.frame,
			rejections=rejections.frame
	)
	return(comparing)		
	
}
#For testing purposes
#source('~/workspace/mutoss/src/BasicFunctions/DummyBigObjects.R')
#test<- list(mu.test.obj.1, mu.test.obj.2)

#-------------- Comparison of adjusted p values -----------------#

mu.compare.adjusted<- function(comparison.list, identify.check=F){
	adjPValues<- comparison.list[['adjusted.pvals']]
	hyp.num<- nrow(adjPValues)
	method.num<- ncol(adjPValues)
	method.names<- factor(colnames(adjPValues))
	method.index<- as.numeric(method.names)
	raw.pValues<-comparison.list[['raw.pValues']]
	
	pvalue.ranks<-raw.pValues$ranks
	
	method.type<-comparison.list[['types']]
	stacked.adjPValues<- unlist(adjPValues, use.names=F)
	x<- rep(pvalue.ranks, method.num)
	method.labels<- rep(method.names, each=hyp.num)
	hyp.labels<- rep(row.names(adjPValues), method.num)
	point.charachters<- rep(method.index, each=hyp.num) #for plotting purposes only
	
	point.size<- hyp.num^(-0.1)
	plot(stacked.adjPValues~x,
			pch=point.charachters,
			ylim=c(0,1),
			cex=point.size,
			xlab='')
	
	the.title<- paste('Adjusted p-values for ',unique(method.type),' controlling procedures')
	title(the.title)			
	
	par(xpd=T)
	
	legend(x=0, y=-0.15,
			horiz=T,
			legend=method.names,
			pch=method.index,
			cex=method.num ^ (-1/4)	)
	par(xpd=F) #reset par to default value	
	
	if(identify.check) {
		identify(stacked.adjPValues~x, labels=hyp.labels )
	}
	##TODO: [JR] Plotting method using colors?	
}

#For testing purposes
#source('~/workspace/mutoss/src/BasicFunctions/DummyBigObjects.R')

#----------- Comparison of critical vales---------- #

mu.compare.critical<- function(comparison.list, identify.check=F){
	method.type<-comparison.list[['types']] #extracting method type
	
	criticalValues<- comparison.list[['criticalValue']]#extracting critical values
	hyp.num<- nrow(criticalValues)
	method.num<- ncol(criticalValues)
	method.names<- factor(colnames(criticalValues))
	method.index<- as.numeric(method.names)
	
	raw.pValues<-comparison.list[['raw.pValues']] #extracting raw palues	
	pvalue.ranks<-raw.pValues$ranks	
	
	stacked.criticalValues<- unlist(criticalValues, use.names=F)
	x<- rep(pvalue.ranks, method.num)
	method.labels<- rep(method.names, each=hyp.num)
	hyp.labels<- rep(row.names(criticalValues), method.num)
	point.charachters<- rep(method.index, each=hyp.num) #for plotting purposes only
	
	point.size<- hyp.num^(-0.1)
	
	plot(stacked.criticalValues~x, 
			pch=point.charachters, 
			ylim=c(0,1),
			xlab='',
			cex=point.size)
	
	the.title<- paste('Critical Values for ',unique(method.type),' controlling procedures')
	title(the.title)			
	
	par(xpd=T)
	legend(
			x=0,
			y=-0.15,
			horiz=T,
			legend=method.names,
			pch=method.index,
			cex=method.num ^ (-1/4)	)
	par(xpd=F) #reset par to default value	
	
	if(identify.check) {
		identify(stacked.criticalValues~x, labels=hyp.labels )
	}
}

#For testing purposes:
#source('~/workspace/mutoss/src/BasicFunctions/DummyBigObjects.R')
#mu.compare.critical(1)
#mu.compare.critical(compare.3, T)

#----- Sumary of comparison-----------#
mu.compare.summary<- function(comparison.list){
	method.type<-comparison.list[['types']] #extracting method type
	error.rates<-comparison.list[['rates']] #extracting error rates
	pi.nulls<-comparison.list[['pi.nulls']] #extracting pi0
	rejections<-comparison.list[['rejections']]
	hyp.num<- nrow(rejections)
	method.num<- ncol(rejections)
	method.names<- factor(colnames(rejections))
	method.index<- as.numeric(method.names)
	
	count.rejections<-apply(rejections, 2, sum)
	
	seperate<- rep('|', method.num)	
	summary<-data.frame(
			method.type, seperate, 
			error.rates, seperate, 
			count.rejections, seperate,
			pi.nulls)
	colnames(summary) <-c(
			'Error Type'," ",
			'Error Rate'," ",
			'Rejections Count', " ",
			'pi_0') 
	
	cat('\n Comparing multipe hypothesis procedures.\n',
			hyp.num, 'hypotheses tested.\n\n')
	
	print(summary)	
}

#For testing purposes :
#source('~/workspace/mutoss/src/BasicFunctions/DummyBigObjects.R')
#mu.compare.summary(1)
#mu.compare.summary(compare.3)