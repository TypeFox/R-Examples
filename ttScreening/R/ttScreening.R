ttScreening<- function(y=y,formula,imp.var,data,iterations=100,sva.method=c("two-step","irw"),
					cv.cutoff=50,n.sv=NULL,train.alpha=0.05,test.alpha=0.05, 
					FDR.alpha=0.05,Bon.alpha=0.05,percent=(2/3),
					linear= c("robust","ls"),vfilter = NULL, B = 5, numSVmethod = "be",rowname=NULL){


	 family <- gaussian()
       if (missing(data)) 
       data <- environment(formula)
 
      mf <- match.call(expand.dots = FALSE)
      m <- match(c("formula", "data", "na.action", "etastart", "mustart", "offset"), names(mf), 0L)
    	mf <- mf[c(1L, m)]
   	mf$drop.unused.levels <- TRUE
	    mf[[1L]] <- quote(stats::model.frame)
	    mf <- eval(mf, parent.frame())
         control <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
    	mt <- attr(mf, "terms")

      Y <- model.response(mf, "any")
	if (is.null(rowname)){
	rownames(y) <- seq(1:nrow(y))
	} else {rownames(y) <- rownames(y)}	

       X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts=NULL,na.action=NULL)
        else matrix(, NROW(Y), 0L)

	
	'%ni%'<-Negate('%in%')
	x.omit<-which(rownames(data) %ni% attributes(X)$dimnames[[1]])
	if(length(x.omit)< 1){y=na.omit(y)}else{y=na.omit(y[,-x.omit])}	
	if(is.null(imp.var)){imp.var<-length(attr(terms(formula,keep.order=TRUE),"term.labels"))}
	imp.var.location<-match(attr(terms(formula,keep.order=FALSE),"term.labels"), attr(terms(formula,keep.order=TRUE),"term.labels"))
	remove.imp.var<-which(imp.var.location %in% imp.var)
	final.imp.var<-which(attributes(X)$assign %in% remove.imp.var)

	if(as.numeric(max(y)<1 && min(y)>0) == 1){edata <- as.matrix(log2(y/(1-y)))}else{edata<-as.matrix(y)}

	train.length<-test.length<-TT.output<-FDR.output<-Bon.output<-final.temp<-freq.temp<-cpg.select<-NULL
	selection<-matrix(rep(0,nrow(edata)*iterations),nrow=nrow(edata),ncol=iterations)
	pvalue.matrix<-matrix(rep(NA,nrow(edata)*iterations),nrow=nrow(edata),ncol=iterations)
	
	mod<-X
	mod0<-X[,-final.imp.var]

	if(is.null(n.sv) || n.sv > 0){
	svobj = sva2(dat=as.matrix(edata),mod,mod0,n.sv=n.sv,method=sva.method)
	n.sv.temp<-svobj$n.sv
	temp=data.frame(svobj$sv)
	location<-NULL
		if(n.sv.temp > 0){
		for (w in 1:n.sv.temp){
			if(length(unique(temp[,w])) == 2){
			location<-c(location,w)}else{location<-location}
		}
		}

		if(is.null(location)){svobj$sv<-temp
		}else {svobj$sv<- as.numeric(temp[,-location])
		print(paste(length(location),"surrogate variables were not included in model fitting due to their lack of information"))
		svobj$n.sv<-ncol(svobj$sv)}

		if(n.sv.temp == 0 | is.null(svobj$n.sv) == TRUE){
		modSv = as.matrix(mod)
		svobj$n.sv <- 0}else{
		modSv = cbind(mod,as.matrix(svobj$sv))
		colnames(modSv)<-c(colnames(mod),paste("SV_",colnames(as.matrix(svobj$sv)),sep=''))
		}
	}else{
		modSv = mod
	}

	#TT#
	if(linear=='robust'){print("warning: robust regression may slow down screening process")}
	x <- 1:ncol(edata)
	length <- round(percent * ncol(edata))
	for (i in 1:iterations){
		set.seed(i)		
		iter<-0
		repeat{
			iter<-iter+1
        		train <- sample(x, size=length, replace=FALSE)
	  		test <- x[-train]

				train.sub<-capture.output({lmfit<-tryCatch(eBayes(lmFit(edata[,train], design=modSv[train,], method=linear))$p.value, warning=function(w)return(train.warn=NULL),error=function(e)return(test.warn=NULL))})
				if(is.null(lmfit)==TRUE){train.temp = NULL;lmfit.test=NULL}else{train.temp<-which(apply(as.matrix(lmfit[,final.imp.var]),1,min)<= train.alpha)}			
				train.length<-c(train.length,length(train.temp))	#number of sites selected at training level
				
				if(length(train.temp) == 1){lmfit.test<-tryCatch(matrix(summary(glm(matrix(edata[train.temp,test],ncol=1)~modSv[test,-1],family=gaussian))$coefficients[,4],nrow=1,ncol(modSv)),warning=function(w)return(test.warn=NULL))}
				if(length(train.temp)> 1){test.sub<-capture.output({lmfit.test<-tryCatch(eBayes(lmFit(edata[train.temp,test], design=modSv[test,], method=linear))$p.value,warning=function(w)return(test.warn=NULL),error=function(e)return(test.warn=NULL))})}
				if(is.null(lmfit.test)==TRUE || length(train.temp) < 1){test.temp=NULL}else{test.temp<-train.temp[(apply(as.matrix(lmfit.test[,final.imp.var]),1,min)<= test.alpha)]
																	test.pvalue<-apply(as.matrix(lmfit.test[,final.imp.var]),1,min)}
			if((is.null(lmfit) == FALSE & is.null(lmfit.test) == FALSE)){break.indicator<-'nobreak';break}	
			if(iter == 10){warning("Not enough variation in variables to ensure valid training and testing samples");break.indicator<-'break';break}															 			
		}

			if(break.indicator == 'break'){break}
			test.length<-c(test.length,length(test.temp))		#number of sites selected at testing level
			if(length(test.temp) < 1){selection[,i]<- 0}else{selection[test.temp,i]<- 1
						  					pvalue.matrix[test.temp,i]<- test.pvalue[(test.pvalue <= test.alpha)]}

	}


		cutoff <- ((cv.cutoff/100)*iterations)
		cpg.select<-rowSums(selection)
		freq.temp<-rowSums(selection)[rowSums(selection) >= cutoff]
		final.temp<-which(rowSums(selection) >= cutoff)


	if(length(final.temp) > 1){
		tt<-eBayes(lmFit(edata[final.temp,],design=modSv,method=linear))
		TT.output<-data.frame(rownames(edata)[final.temp],freq.temp,tt$coefficients,tt$p.value)
	}
	if(length(final.temp) == 1){
		if(linear == 'ls'){
		tt<- summary(glm(matrix(edata[final.temp,],ncol=1)~ modSv[,-1],family=gaussian))
		TT.output<-data.frame(rownames(edata)[final.temp],freq.temp,t(tt$coefficients[,1]),t(tt$coefficients[,4]))
		}else{tt<-summary(rlm(matrix(edata[final.temp,],ncol=1)~modSv[,-1],family=gaussian))
		TT.output<-data.frame(rownames(edata)[final.temp],t(tt$coefficients[-1,1]),t(2*pt(-abs(as.numeric(tt$coefficients[-1,3])),df=tt$df[2])))}
	}
	if(length(final.temp) < 1){
		TT.output<- t(rep("NA",(ncol(modSv)*2 + 2)))
	}
	colnames(TT.output)<-c("Row_Ind_Select","Selection Prop",paste("Coeff",colnames(modSv)),paste("Pvalue",colnames(modSv)))

		
	#FDR/Bon#
		lmfit<-eBayes(lmFit(edata, design=modSv, method=linear))
		lmfit.min<-apply(as.matrix(lmfit$p.value[,final.imp.var]),1,min)	

		lmFitFDR.rob<-which(p.adjust(lmfit.min,method="fdr")<= FDR.alpha)
		lmFitBon.rob<-which(p.adjust(lmfit.min,method="bonferroni")<= Bon.alpha)

	if(length(lmFitFDR.rob) > 1){
		FDR.output<-data.frame(rownames(edata)[lmFitFDR.rob], lmfit$coefficients[lmFitFDR.rob,], lmfit$p.value[lmFitFDR.rob,])}
	if(length(lmFitFDR.rob) == 1){
		if(linear == 'ls'){
		ff<- summary(glm(matrix(edata[lmFitFDR.rob,],ncol=1)~modSv[,-1],family=gaussian))
		FDR.output<-data.frame(rownames(edata)[lmFitFDR.rob], t(ff$coefficients[,1]), t(ff$coefficients[,4]))
		}else{
		ff<-summary(rlm(matrix(edata[lmFitFDR.rob,],ncol=1)~modSv[,-1],family=gaussian))
		FDR.output<-data.frame(rownames(edata)[lmFitFDR.rob],t(ff$coefficients[-1,1]),t(2*pt(-abs(as.numeric(ff$coefficients[-1,3])),df=ff$df[2])))}
		}
	if(length(lmFitFDR.rob) < 1){
		FDR.output<- t(rep("NA",(ncol(modSv)*2 + 1)))}

	if(length(lmFitBon.rob) > 1){
		Bon.output<-data.frame(rownames(edata)[lmFitBon.rob], lmfit$coefficients[lmFitBon.rob,], lmfit$p.value[lmFitBon.rob,])}
	if(length(lmFitBon.rob) == 1){
		if(linear == 'ls'){
		bb<- summary(glm(matrix(edata[lmFitBon.rob,],ncol=1)~modSv[,-1],family=gaussian))
		Bon.output<-data.frame(rownames(edata)[lmFitBon.rob], t(bb$coefficients[,1]), t(bb$coefficients[,4]))
		}else{
		bb<-summary(rlm(matrix(edata[lmFitBon.rob,],ncol=1)~modSv[,-1],family=gaussian))
		Bon.output<-data.frame(rownames(edata)[lmFitBon.rob],t(bb$coefficients[-1,1]),t(2*pt(-abs(as.numeric(bb$coefficients[-1,3])),df=bb$df[2])))}}
	if(length(lmFitBon.rob) < 1){
		Bon.output <- t(rep("NA",(ncol(modSv)*2 + 1)))}
	
	colnames(FDR.output)<-colnames(Bon.output)<-c("Row_Ind_Select",paste("Coeff",colnames(modSv)),paste("Pvalue",colnames(modSv)))


output=list(TT.cpg=rownames(edata)[final.temp],train.cpg=train.length,test.cpg=test.length,
			selection=selection,pvalue.matrix=pvalue.matrix,TT.output=TT.output, FDR.output=FDR.output, Bon.output=Bon.output,
			FDR.cpg=rownames(edata)[lmFitFDR.rob],Bon.cpg=rownames(edata)[lmFitBon.rob])

output
}
