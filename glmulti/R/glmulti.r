# S4 class definition
setClass("glmulti", representation(name="character", params="list", nbmods="integer", crits="numeric", 
		 K="integer", formulas="list", call="call", adi="list",objects="list"))


	
# S3 function definitions
# utilities
summary.glmulti<-function(object, ...)
{
	who=object@formulas[object@crits==object@crits[1]]
	ww = exp(-(object@crits - object@crits[1])/2)
	ww=ww/sum(ww)
	cucu=function(i) sum(ww[1:i])
	wwc=lapply(1:length(ww),cucu)
	ret=list(name=object@name, method=object@params$method, fitting=object@params$fitfunction, crit=object@params$crit, 
			level=object@params$level, marginality=object@params$marginality,confsetsize=object@params$confsetsize, 
			bestic=object@crits[1], icvalues=object@crits ,	bestmodel=deparse(who[[1]]), modelweights=ww)
	if (object@params$method == "g") {
		ret=c(ret, generations=object@params$generations, elapsed= (object@params$elapsed)/60)
	}
	hasobj = try(object@objects,TRUE)
	hasobj = (!inherits(hasobj,"try-error") && length(object@objects) > 0) 
	ret=c(ret, includeobjects=hasobj) 
	ret
}

print.glmulti<-function(x, ...)
{
	write(x@name, file="")
	write(paste("Method: ",x@params$method," / Fitting: ",x@params$fitfunction," / IC used: ", 
			x@params$crit,  sep=""), file="")
	write(paste("Level: ", x@params$level, " / Marginality: ", x@params$marginality,  sep=""), file="")
	write(paste("From ",length(x@crits)," models:",sep=""), file="")
	write(paste("Best IC: ",x@crits[1],sep=""), file="")
	write(paste("Best model:",sep=""), file="")
	who=x@formulas[x@crits==x@crits[1]]
	for (w in who)
		print(deparse(w))
	ww = exp(-(x@crits - x@crits[1])/2)
	ww=ww/sum(ww)
	write(paste("Evidence weight: ",ww[1], sep=""), file="")
	write(paste("Worst IC: ",x@crits[length(x@crits)],sep=""), file="")
	write(paste(length(which(x@crits<=x@crits[1]+2)), " models within 2 IC units.", sep=""), file="")
	cucu=function(i) sum(ww[1:i])
	wwc=lapply(1:length(ww),cucu)
	write(paste(length(which(wwc<=0.95)), " models to reach 95% of evidence weight.", sep=""), file="")
	if (x@params$method == "g") {
		write(paste("Convergence after ", x@params$generations, " generations.",  sep=""), file="")
		write(paste("Time elapsed: ",  (x@params$elapsed), " minutes.", sep=""), file="")
	}
}

plot.glmulti<-function(x, type="p", ...) 
{
	if ( type == "p") {
		plot(x@crits,xlab="Best models",
				ylab=paste("Support (",x@params$crit,")"),pch=19, main="IC profile", ...)
		
		abline(h=x@crits[1]+2,col="red")
	} else if (type == "w") {
		ww = exp(-(x@crits - x@crits[1])/2)
		ww=ww/sum(ww)
		plot(ww,xlab="Best models",
				ylab=paste("Evidence weight (",x@params$crit,")"),pch=19, main="Profile of model weights", ...)
		cucu=function(i) sum(ww[1:i])
		wwc=lapply(1:length(ww),cucu)
		abline(v=min(which(wwc>=0.95)),col="red")
	} else if (type=="r") {
		if (length(x@objects)) {
		# shows some diagnostics of the fit
		windows(21,7)
		par(mfrow=c(2,min(length(x@crits), 5)))
		for (k in 1:min(length(x@crits), 5)) 
			plot(x@objects[[k]],which=c(1), main=deparse(x@formulas[[k]]),...)
		for (k in 1:min(length(x@crits), 5)) 
			plot(x@objects[[k]],which=c(2),...)
		} else warning("Unavailable: use includeobjects=T when calling glmulti.")		
	} else if (type=="s") {
		# plots variable (i.e. terms) importances
		ww = exp(-(x@crits - x@crits[1])/2)
		ww=ww/sum(ww)
		# handle synonymies for interactions
		# this translates to unique notations (e.g. x:y and y:x are the same)
		clartou=function(x) {
			sort(strsplit(x, ":")[[1]])-> pieces
			if (length(pieces)>1) paste(pieces[1],":",pieces[2],sep="")
			else x
			}
		# list terms in models
		tet = lapply(x@formulas, function(x) sapply(attr(delete.response(terms(x)),"term.labels"), clartou))
		# all unique terms
		unique(unlist(tet))-> allt
		# importances
		sapply(allt, function(x) sum(ww[sapply(tet, function(t) x%in%t)]))-> imp
		# draw
		par(oma=c(0,3,0,0))
		barplot(sort(imp),xlab="",xlim=c(0,1), ylab="",horiz=T,las=2, names.arg=allt[order(imp)],main="Model-averaged importance of terms", ...)
		abline(v=0.8, col="red")
	} else	warning("plot: Invalid type argument for plotting glmulti objects.")
}


# model averaging: done through coef
coef.glmulti <- function(object, select="all", varweighting="Buckland", icmethod="Lukacs", alphaIC=0.05, ...) 
{
	ww = exp(-(object@crits - object@crits[1])/2)
	ww=ww/sum(ww)
	cucu=function(i) sum(ww[1:i])
	wwc=lapply(1:length(ww),cucu)
	# what models are to be used
	whom=c()
	if (length(select)>1)
		whom=select
	else if (select=="all")
			whom=1:length(object@crits)
		else if (is.integer(select))
				whom = 1:select
			else if (is.numeric(select)) {
				if (select <=1) 
					whom = which(wwc<=select)
				else 
					whom = which(object@crits <= (object@crits[1]+select))
		}
	mods = object@crits[whom]
	formo = object@formulas[whom]
	hasobj = try(object@objects,TRUE)
	if (!inherits(hasobj,"try-error") && length(object@objects) > 0) {
	# model objects included: no need to refit
	coffee = object@objects[whom]
	} else {
	# refit models
	coffee=list()
	for (i in formo) {
		ff=object@params$fitfunction
		cak=as.call(list(substitute(match.fun(ff)), formula=i, data=object@call$data))
		if (length(object@adi)>=1)
			for (j in 1:length(object@adi)) {
				cak[[length(names(cak))+1]] = object@adi[[j]]
				names(cak)[length(names(cak))] = names(object@adi)[j]
			}
		modf=eval(cak)
		coffee=c(coffee,list(modf))
	}
	}

	# construct list of coefficients
	if (length(coffee)==1) {
		# only one model ! Do conditional inference for continuity
		warning("Only one candidate: standard conditional inference was performed.")
		return(coef(coffee[[1]]))
	}
	coke=lapply(coffee,getfit)
   	namou=unique(unlist(lapply(coffee,function(x) dimnames(getfit(x))[[1]])))
	# this equates synonymous notations (e.g. x:y and y:x)
	unique(sapply(namou, function(x) {
		sort(strsplit(x, ":")[[1]])-> pieces
		if (length(pieces)>1) paste(pieces[1],":",pieces[2],"<>",pieces[2],":",pieces[1],sep="")
		else x
		}))-> codamou 
	namou = sapply(codamou, function(x)  {
		sort(strsplit(x, "<>")[[1]])-> pieces
		if (length(pieces)>1) (pieces[1])
		else x
		})
	namou2 = sapply(codamou, function(x)  {
		sort(strsplit(x, "<>")[[1]])-> pieces
		if (length(pieces)>1) pieces[2]
		else x
		})
		
	coconutM=matrix(0,length(formo),length(namou))
	coconutSE=matrix(0,length(formo),length(namou))
	coconutN = numeric(length(namou))
	# get values, deviations, presence/absence and sample sizes for all models
	matchou=function(quiqui) {
		match(quiqui,namou)-> w1
		if (is.na(w1)) {
			match(quiqui,namou2)
		} else w1
	}
	gettou=function(i) {
		ele=coke[[i]]
		nana = dimnames(ele)			
		if (length(nana)==2) nana=nana[[1]]
		mimi=numeric(4*length(namou))
		if (length(nana) > 0) {
			for (k in 1:(length(nana))) {
				mimi[matchou(nana[k])]=ele[k,1]
				mimi[matchou(nana[k])+length(namou)]=ele[k,2]
				mimi[matchou(nana[k])+2*length(namou)]=1
				mimi[matchou(nana[k])+3*length(namou)]=ele[k,3]
			}
		} 
		return(mimi)
	}
	lol=sapply(lapply(1:length(coke),gettou),rbind)
	

	coconutM = matrix(unlist(t(lol[1:length(namou),])),nrow=length(whom))
	coconutSE = matrix(unlist(t(lol[(1:length(namou))+length(namou),])),nrow=length(whom))
	# NA are set to zero
	coconutM[is.na(coconutM)]=0
	coconutN =  matrix(unlist(t(lol[(1:length(namou))+2*length(namou),])),nrow=length(whom))
	
	# handle degrees of freedom
	coconutDF = matrix(unlist(t(lol[(1:length(namou))+3*length(namou),])),nrow=length(whom))
	modelsdf= unlist(apply(coconutDF,1,max))
	
	nene = matrix(rep(1,length(whom))%*%coconutN, ncol=1, dimnames=list( namou, c("Nb models")))
	# construct weight vectors
	waou=ww[whom]/sum(ww[whom])
	waouv=waou
	for (i in 2:length(namou)) waouv=rbind(waouv,waou) 
	waouv= t(waouv)*coconutN
	totwaou = waou%*%coconutN
	# weight estimates
	averest = matrix((rep(1,length(whom))%*%(waouv*coconutM)), ncol=1, dimnames=list( namou, c("Estimate")))
	weighty =  matrix(totwaou, ncol=1, dimnames=list( namou, c("Importance")))
	# weight variances
	if (varweighting=="Johnson") {
		squaredevs = waou%*%(((coconutM-t(matrix(rep(averest,length(whom)), length(namou), length(whom))))^2))
		condivars =  waou%*%((coconutSE^2))
		avervar = matrix(condivars+squaredevs, ncol=1, dimnames=list( namou, c("Uncond. variance")))
	} else if (varweighting=="Buckland") {
		squaredevs = ((coconutM-t(matrix(rep(averest,length(whom)), length(namou), length(whom)))))^2
		condivars = coconutSE^2
		avervar = matrix(((waou)%*%(sqrt(squaredevs+condivars)))^2, ncol=1, dimnames=list( namou, c("Uncond. variance")))
	} else 
		avervar = matrix(rep(NA,length(namou)), ncol=1, dimnames=list( namou, c("Uncond. variance")))

	# now move on to confidence intervals
	if (icmethod=="Burnham") {
		# uses Burnham & Anderson (2002) suggestion
		stuvals = (as.numeric(lapply(modelsdf,  function(x) qt(1-alphaIC/2, x)))/qnorm(1-alphaIC/2))^2
		adjsem= (matrix(rep(stuvals, length(namou)), nrow=length(whom)))*coconutSE^2
		adjsem = adjsem + (coconutM-t(matrix(rep(averest, length(whom)), ncol=length(whom))))^2
		adjse = qnorm(1-alphaIC/2)*(waou%*%sqrt(adjsem))
		uncondIC = matrix(adjse, ncol=1,  dimnames=list(namou, c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
	} else if (icmethod=="Lukacs") {
		# uses Lukacs et al. (2008) student-like method
		# get degrees of freedom for each model
		averddf = sum(waou*modelsdf)
		uncondIC = matrix(sqrt(avervar)*qt(1-alphaIC/2,averddf), ncol=1,  dimnames=list(namou, c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
	} else {
		# uses standard gaussian interval 
		uncondIC = matrix(sqrt(avervar)*qnorm(1-alphaIC/2), ncol=1,  dimnames=list(namou, c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
	}

	averaging = cbind(averest, avervar, nene, weighty, uncondIC)
	ordonat = order(weighty[,1])
	# remove null component if necessary
	which(namou=="NULLOS")-> nooz
	if (length(nooz)>0) ordonat=setdiff(ordonat,nooz)
	averaging[ordonat,]
}


# model averaged prediction
predict.glmulti <- function(object, select="all", newdata=NA, se.fit=FALSE, varweighting="Buckland", icmethod="Lukacs", alphaIC=0.05,...) 
{
	ww = exp(-(object@crits - object@crits[1])/2)
	ww=ww/sum(ww)
	cucu=function(i) sum(ww[1:i])
	wwc=lapply(1:length(ww),cucu)
	# what models are to be used
	whom=c()
	if (length(select)>1)
		whom=select
	else if (select=="all")
			whom=1:length(object@crits)
		else if (is.integer(select))
				whom = 1:select
			else if (is.numeric(select)) {
				if (select <=1) 
					whom = which(wwc<=select)
				else 
					whom = which(object@crits <= (object@crits[1]+select))
		}
	mods = object@crits[whom]
	formo = object@formulas[whom]
	hasobj = try(object@objects,TRUE)
	if (!inherits(hasobj,"try-error") &&length(object@objects) > 0) {
	# model objects included: no need to refit
	coffee = object@objects[whom]
	} else {
	# refit models
	coffee=list()
	for (i in formo) {
		ff=object@params$fitfunction
		cak=as.call(list(substitute(match.fun(ff)), formula=i, data=object@call$data))
		if (length(object@adi)>=1)
			for (j in 1:length(object@adi)) {
				cak[[length(names(cak))+1]] = object@adi[[j]]
				names(cak)[length(names(cak))] = names(object@adi)[j]
			}
		modf=eval(cak)
		coffee=c(coffee,list(modf))
	}
	}
	
	if (length(coffee)==1) {
		# only one model ! Do conditional inference for continuity
		warning("Only one candidate: standard conditional prediction was performed.")
		return(predict(coffee[[1]],se.fit=se.fit,...))
	}
	
	waou=ww[whom]/sum(ww[whom])
	# make predictions
	predicts=list()
	for (i in 1:length(formo)) {
		if(!is.data.frame(newdata))	{
				arlette = predict(coffee[[i]],se.fit=se.fit, ...)
        		predicts = c(predicts, list(arlette))
        		} else {
        		arlette = predict(coffee[[i]], newdata=newdata, se.fit=se.fit, ...)
       			 predicts = c(predicts, list(arlette))
		}
	}
	# handle average predictions
	if (se.fit) { 
	preds=lapply(predicts,function(x) x[[1]])
	}  else preds = predicts
	nbpo = length(preds[[1]])
	all = t(matrix(unlist(preds), nrow=nbpo))
	dimnames(all)=list(c(), names(preds[[1]]) )
	# handle NA values
	nana = lapply(preds, is.na)
	nbok = numeric(nbpo)
	for (i in 1:length(preds)) {
		nbok = nbok + nana[[i]]
		preds[[i]][is.na(preds[[i]])] = 0
		}
	minou = matrix(waou%*%t(matrix(unlist(preds),nrow=nbpo)), dimnames=list(names(preds[[1]]),c() )) 
	
	# handle variances if appropriate
	mvar=NULL
	if (se.fit) {
		waouv=waou
		if (nbpo>1) for (i in 2:nbpo) waouv=rbind(waouv,waou) 
		waouv= matrix(t(waouv),nrow=length(whom),ncol=nbpo)
		predse=lapply(predicts,function(x) x[[2]])
		allse = t(matrix(unlist(predse), nrow=nbpo))
		dimnames(allse)=list(c(), names(predse[[1]]) )
		# handle NA values
		nana = lapply(predse, is.na)
		nbok2 = numeric(nbpo)
		for (i in 1:length(predse)) {
			nbok2 = nbok2 + nana[[i]]
			predse[[i]][is.na(predse[[i]])] = 0
		}
		# get degrees of freedom
		modelsdf = unlist(lapply(coffee,function(x) max(getfit(x)[,3])))
		
		# compute variance components
			if (varweighting=="Buckland") {
			squaredevs = ((all-t(matrix(rep(minou,length(whom)), nbpo, length(whom)))))^2
			condivars = allse^2
			avervar = matrix(((waou)%*%(sqrt(squaredevs+condivars)))^2, ncol=1, dimnames=list( names(predse[[1]]), c("Uncond. variance")))
		} else if (varweighting=="Johnson") {
			squaredevs = waou%*%(((all-t(matrix(rep(minou,length(whom)), nbpo, length(whom))))^2))
			condivars =  waou%*%((allse^2))
			avervar = matrix(condivars+squaredevs, ncol=1, dimnames=list(  names(predse[[1]]), c("Uncond. variance")))
		}	
		# move on to confidence intervals
		if (icmethod=="Burnham") {
			# uses Burnham & Anderson (2002) suggestion
			stuvals = (as.numeric(lapply(modelsdf,  function(x) qt(1-alphaIC/2, x)))/qnorm(1-alphaIC/2))^2
			adjsem= (matrix(rep(stuvals, nbpo), nrow=length(whom)))*allse^2
			adjsem = adjsem + (all-t(matrix(rep(minou, length(whom)), ncol=length(whom))))^2
			adjse = qnorm(1-alphaIC/2)*(waou%*%sqrt(adjsem))
			uncondIC = matrix(adjse, ncol=1,  dimnames=list( names(predse[[1]]), c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
		} else if (icmethod=="Lukacs") {
			# uses Lukacs et al. (2010) student-like method
			# get degrees of freedom for each model
			averddf = sum(waou*modelsdf)
			uncondIC = matrix(sqrt(avervar)*qt(1-alphaIC/2,averddf), ncol=1,  dimnames=list( names(predse[[1]]), c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
		} else {
			# uses standard gaussian interval 
			uncondIC = matrix(sqrt(avervar)*qnorm(1-alphaIC/2), ncol=1,  dimnames=list( names(predse[[1]]), c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
		}
		mvar = cbind(avervar, uncondIC)
		
	}
	
	list(averages = t(minou), variability = mvar, omittedNA = sum(nbok))
	
}




#  S4 function definitions
		
# write redefinition
setGeneric("write", function(x, file = "data",   ncolumns=if(is.character(x)) 1 else 5, append = FALSE, sep = " ") standardGeneric("write"))
setMethod("write","glmulti", function(x, file , ncolumns, append, sep)
{
		if (!missing(file)&&length(grep(file, "|object"))==1) {
			ir = gsub("\\|object","",file)
			if (ir=="")
				saveRDS(x, file=x@name)
			else 
				saveRDS(x, file=ir)
		} else {
			concato=cbind(data.frame(K=x@K),data.frame(IC=x@crits), data.frame(Models=as.character(x@formulas)))
			if (missing(file))
				write.table(concato, file = paste(x@name, ".txt"),  append = append, sep = sep )
			else
				write.table(concato, file = file,  append = append, sep = sep )
		}
})
		
# accessing fitted models for coefficients
setGeneric("getfit", function(object, ...) standardGeneric("getfit"))

setMethod("getfit","ANY", function(object, ...)
{
	summ = summary(object)
	summ1 = summ$coefficients
	didi=dimnames(summ1)
	if (is.null(didi[[1]])) {
			summ1 = matrix(rep(0,2), nrow=1, ncol=2, dimnames=list(c("NULLOS"),list("Estimate","Std. Error")))
			return(cbind(summ1, data.frame(df=c(0))))
	}
	summ1=summ1[,1:2]
	if (length(dim(summ1))==0) {
		didi = dimnames(summ$coefficients)
		summ1=matrix(summ1, nrow=1, ncol=2, dimnames=list(didi[[1]],didi[[2]][1:2]))
	}
	return(cbind(summ1, data.frame(df=rep(summ$df[2], length(summ$coefficients[,1])))))
	
})

# consensus method
setGeneric("consensus", function (xs, confsetsize=NA, ...)  standardGeneric("consensus"))

setMethod("consensus", signature(xs="list"), function (xs, confsetsize, ...)
{
	lespaul = list()
	for (i in xs) {
	if (class(i)=="glmulti")
		lespaul=c(lespaul, i)
	else if (class(i)=="character") {
			paul = .readRDS(file=i)
			if (class(paul)=="glmulti")
				lespaul=c(lespaul,paul)
		}
	}
	takeobjects = (min(sapply(lespaul, function (x) length(x@objects)))>0) 
	neo = new ("glmulti")
	neo@objects=list()
	conca = function (x) {
		cbind(data.frame(K=x@K),data.frame(IC=x@crits), data.frame(formulas=as.character(x@formulas)))
	}
	tot=lapply(lespaul, conca)
	tota = tot[[1]]
	for (h in 2: length(tot))
		tota = rbind(tota, tot[[h]])
	tobekept = !duplicated(tota$formulas)
	rearrange = order(tota$IC[tobekept])
	tot = tota
	tot = tot[tobekept,]
	tot = tot[rearrange,]
	
	if (takeobjects) {
	concab = function (x) {
		x@objects
	}
	tob = lapply(lespaul, concab)
	toba = tob[[1]]
	for (h in 2: length(tob))
		toba = c(toba, tob[[h]])
	tob=toba
	tob = tob[tobekept]
	tob = tob[rearrange]
	}
	
	if (is.na(confsetsize) || length(tot$K)<confsetsize) {
		if (!is.na(confsetsize)&&length(tot$K)<confsetsize)
			warning("Could not gather enough models.")
		neo@K = tot$K
		neo@formulas =  as.list(lapply(unlist(lapply(tot$formulas, function(j) as.character(j))), function(ff) as.formula(c(ff))))
		neo@crits = tot$IC
		if (takeobjects) neo@objects = tob
	} else {
		tot = tot[1:confsetsize, ]
		neo@K = tot$K
		neo@formulas = as.list(lapply(unlist(lapply(tot$formulas, function(j) as.character(j))), function(ff) as.formula(c(ff))))
		neo@crits = tot$IC
		if (takeobjects) neo@objects = tob[1:confsetsize]
	}
	neo@call = lespaul[[1]]@call
	neo@adi = lespaul[[1]]@adi
	neo@name = paste("consensus of ", length(lespaul), "-", lespaul[[1]]@name, sep="")
	neo@params = lespaul[[1]]@params
	return (neo)
	})


# support for coxph
logLik.coxph <- function(object,...) {
    y <-  object$loglik[length(object$loglik)]
    class(y) <- "logLik"
    attr(y,'df') <- if(is.null(object$coef)) 0 else sum(!is.na(object$coef))
    return(y)
}

setMethod("getfit",signature(object="coxph"), function(object, ...)
{
	summ = summary(object)
	summ1 = summ$coefficients[,c(1,3)]
	if (length(dim(summ1))==0) {
		didi = dimnames(summ$coefficients)
		summ1=matrix(summ1, nrow=1, ncol=2, dimnames=list(didi[[1]],didi[[2]][c(1,3)]))
	}
	df = object$n-attr(logLik(object),"df")
	return(cbind(summ1, data.frame(df=rep(df, length(summ$coefficients[,1])))))
})

setMethod("getfit",signature(object="coxph.null"), function(object, ...)
{
	return(NULL)
})





# information criteria
setGeneric("aicc", function(object, ...) standardGeneric("aicc"))
setGeneric("aic", function(object, ...) standardGeneric("aic"))
setGeneric("bic", function(object, ...) standardGeneric("bic"))
setGeneric("qaic", function(object, ...) standardGeneric("qaic"))
setGeneric("qaicc", function(object, ...) standardGeneric("qaicc"))
setMethod("aicc", "ANY", function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
	n= nobs(object)
	return(-2*as.numeric(liliac[1]) + 2*k*n/max(n-k-1,0))
})

setMethod("bic", signature(object="ANY"), function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
	n= nobs(object)
	return(-2*as.numeric(liliac[1]) + k*log(n))
})

setMethod("aic", "ANY",  function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
	return(-2*as.numeric(liliac[1])+2*k)
})
setMethod("qaic", "ANY",  function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
 	return(-2*as.numeric(liliac[1]) / getOption("glmulti-cvalue") + 2 * k)

})

setMethod("qaicc", "ANY",  function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
	n= nobs(object)
 	return(-2*as.numeric(liliac[1]) / getOption("glmulti-cvalue") + 2*k*n/max(n-k-1,0))

})


# printing a quick table of IC weights
# kudos to J. Byrnes
setGeneric("weightable", function(object, ...) standardGeneric("weightable"))
setMethod("weightable", "glmulti",  function(object, ...)
{
	#get the summary
	obj_summary<-summary(object)
	#get the functions for output
	funcs<-as.character(object@formulas)
	ret<-data.frame(model=funcs, ic=obj_summary$icvalues, weights=obj_summary$modelweights)
	#a little renaming
	names(ret)[which(names(ret)=="ic")]<-obj_summary$crit
	return(ret)
})






# generic glmulti function
setGeneric("glmulti", function(y, xr, data, exclude=c(), name="glmulti.analysis", intercept=TRUE, marginality=FALSE , bunch = 30, chunk=1, chunks=1, 
		level=2, minsize=0, maxsize=-1, minK=0, maxK=-1, method="h",crit="aic",confsetsize=100,popsize=100,mutrate=10^-3,
		sexrate=0.1,imm=0.3, plotty=TRUE, report=TRUE, deltaM=0.05, deltaB=0.05, conseq=5, fitfunction="glm", resumefile = "id", includeobjects=TRUE,  ...) 
		standardGeneric("glmulti"))


	
# interfaces for formulas/models: calls with missing xr argument
setMethod("glmulti","missing",
function(y, xr, data, exclude, name, intercept, marginality , bunch, chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile, includeobjects,  ...) 
{
	write("This is glmulti 1.0.7, Apr. 2013.",file="")
})

setMethod("glmulti",
def=function(y, xr, data, exclude, name, intercept, marginality ,bunch, chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,includeobjects, ...) 
{
	stop("Improper call of glmulti.")
})


setMethod("glmulti", signature(y="formula", xr= "missing", exclude="missing"),
function(y, xr, data, exclude, name, intercept, marginality , bunch, chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile, includeobjects, ...) {
	if (missing(data))
		tete = terms(y)
	else
		tete = terms(y, data=data)
	oo = attr(tete,"order")
	dep = as.character(attr(tete,"variables"))[2]
	int = attr(tete,"intercept")
	preds = as.character(attr(tete,"variables"))[-(1:2)]
	if (level==2 && max(oo)>1) {
		# get all possible interactions
		interac = attr(tete,"term.labels")[oo==2]
		neotete = terms(as.formula(paste("h~",paste(preds, collapse="*"))))
		neointerac= attr(neotete,"term.labels")[attr(neotete,"order")==2]
		# get exclusions
		for (i in interac)
			neointerac=neointerac[neointerac!=i]
		# same for main effects
		mama = attr(tete,"term.labels")[oo==1]
		exma = preds
		for (j in mama)
			exma = exma[exma!=j]
		exma = c(exma,neointerac)
	} else {
		preds = attr(tete,"term.labels")[oo==1]
		exma=c(1)
	}
	call = match.call()
	call[[match("y", names(call))]] = dep
	call[[length(names(call))+1]] = preds
	names(call)[length(names(call))] ="xr"
		call[[length(names(call))+1]] = exma
	names(call)[length(names(call))] ="exclude"
	
	if (missing(data)) {
		call[[length(names(call))+1]] = environment(y)
		names(call)[length(names(call))] ="data"
	}
	eval(call)
})

setMethod("glmulti", signature(y="character", xr= "missing", exclude="missing"),
function(y, xr, data, exclude, name, intercept, marginality , bunch, chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile, includeobjects, ...) {
	call = match.call()
	call[[match("y", names(call))]] = as.formula(y)
	eval(call)
})

setMethod("glmulti", signature(y="glm", xr= "missing", exclude="missing"),
function(y, xr, data, exclude, name, intercept, marginality ,bunch, chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile, includeobjects, ...) {
	call = match.call()
	call[[match("y", names(call))]] = formula(y)
	call[[length(names(call))+1]] = "glm"
	names(call)[length(names(call))] = "fitfunction"
	if (length(y$data)) {
		call[[length(names(call))+1]] = y$data
		names(call)[length(names(call))] = "data"
	}
	eval(call)
})

setMethod("glmulti", signature(y="lm", xr= "missing", exclude="missing"),
function(y, xr, data, exclude, name, intercept, marginality , bunch, chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile, includeobjects, ...) {		
	call = match.call()
	call[[match("y", names(call))]] = formula(y)
	call[[length(names(call))+1]] = "lm"
	names(call)[length(names(call))] = "fitfunction"
	if (length(summary(y)$call$data)) {
		call[[length(names(call))+1]] = summary(y)$call$data
		names(call)[length(names(call))] = "data"
	}
	eval(call)
})

# convenience function : pass a list of fitted models
setMethod("glmulti", signature(y="list"),
function(y, xr, data, exclude, name, intercept, marginality ,bunch, chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,includeobjects, ...) 
{
	# go for exhaustive screening 
	if (report) write("TASK: Exhaustive screening of candidate models.",file="")
	# init output object
	resulto <- new("glmulti")
	if (chunks>1) {
		resulto@name = paste(name,"_",chunk,".",chunks,sep="")
	} else 
		resulto@name = name
	resulto@params = list(name=name, intercept=intercept, marginality=marginality ,bunch=bunch, chunk=chunk, chunks=chunks, 
			level=level, minsize=minsize, maxsize=maxsize, minK=minK, maxK=maxK, method="h",crit=crit, confsetsize=confsetsize, fitfunction=fitfunction)
	resulto@call=match.call()
	resulto@adi=list(...)

	if (report) write("Evaluating models...",file="")
	flush.console()
	nbmodels = length(y)
	confsetsize=nbmodels
	lesCrit<-numeric(confsetsize)
	lesK<-vector('integer',confsetsize)
	lesForms<-vector('character',confsetsize)
	
	# functions for support and fit
	support = match.fun(crit)
	if (is.function(crit))
		crit = deparse(substitute(crit))
	fitfunc = match.fun(fitfunction)
	if (is.function(fitfunction)) 
  		fitfunction <- deparse(substitute(fitfunction))
	
	# proceed with model evaluation
	modeval=function(momo) {
			cricriT<-support(y[[momo]])
			lesFormsT=as.character(as.expression(terms(formula(y[[momo]]))))
			lesKT=attr(logLik(y[[momo]]),"df")
			list(cricriT,lesFormsT,lesKT)
	}
	sel=nbmodels
	lapply(1:nbmodels, modeval)-> temporat
	lesCrit=sapply(temporat, function(x) x[[1]])
	lesForms=sapply(temporat, function(x) x[[2]])
	lesK=sapply(temporat, function(x) x[[3]])
	
	
	# over !
	if (report) write("Completed.",file="")
	# prepare and return glmulti object
	reglo <- order(lesCrit[1:sel])
	resulto@crits = lesCrit[1:sel][reglo]
	resulto@formulas = lapply(lesForms[1:sel][reglo], as.formula)
	resulto@K = as.integer(lesK[1:sel][reglo])
	resulto@nbmods = as.integer(sel)
	if (includeobjects) resulto@objects = y else resulto@objects = list()
	return(resulto)

})



	

# workhorse function (call with building blocks and exclusions, as in earlier versions 0.5-x)
setMethod("glmulti", signature(y="character", xr= "character"),
function(y, xr, data, exclude, name, intercept, marginality , bunch, chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile, includeobjects, ...) 
{
	if (report) write("Initialization...",file="")

	# some general constants
	if (method=="h") DELTAD=50
	else DELTAD=10


	# handle terms
	# remove duplicates if any
	x<-unique(xr)
	excluzt<-unique(exclude)
	# "qazxcww" is added to avoid zero length arrays that confuse R. It is always ignored.
	databis = model.frame(as.formula(paste(y,"~", paste(x, sep="", collapse="+"),sep="")), data = data)
	ssize<-length(databis[,1])
	xc<- c("qazxcww")
	xq<- c("qazxcww")
	nc<-0
	nq<-0
	nana=names(databis)
	for (i in x) 
		if (is.factor(databis[[which(nana==i)]])) {
			nc<-nc+1
			xc<-c(xc,i)
		} else {
			nq<-nq+1
			xq<-c(xq,i)
		}
	excluz = c("qazxcww",excluzt)
	
	# functions for support and fit
	support = match.fun(crit)
	if (is.function(crit))
		crit = deparse(substitute(crit))
	fitfunc = match.fun(fitfunction)
	if (is.function(fitfunction)) 
  		fitfunction <- deparse(substitute(fitfunction))
	
	
	# instanciate java object
	if (method=="r") {
		# resume a genetic algorithm
		write("Restoring from files...",file="")
		if (resumefile=="id") 
			molly <- .jcall("glmulti/Resumator","Lglmulti/ModelGenerator;","resto",name)
		else molly <- .jcall("glmulti/Resumator","Lglmulti/ModelGenerator;","resto",resumefile)
		if (!.jfield(molly,"Z","ok")) {
			warning("!Could not restore analysis from files.")
			return(-1)
		}

	} else {
		# brand new analysis
		molly <- .jnew("glmulti/ModelGenerator",y,.jarray(xc),.jarray(xq),.jarray(excluz),as.integer(level),as.integer(minsize),
				as.integer(maxsize),as.integer(minK),as.integer(maxK),intercept,marginality)
		if (!.jfield(molly,"Z","ok")) {
			warning("!Oversized candidate set.")
			return(-1)
		}
		# check for the number of predictors
		okidok=T
		if (method!="l") 
		if (level==1 ) {
			okidok = (nc<=32 && nq<=32)
		} else {
			okidok =  (nc<16 && nq<16 && nc*nq<128)
		}
		if (!okidok) {
			warning("!Too many predictors.")
			return(-1)
		}

		# informs java object of the number of levels of factors, if any
		if (nc>0) {
			flevs = integer(nc)
			for (i in 1:nc) 
				flevs[i]<-as.integer(nlevels(databis[,which(nana==xc[i+1])]))
			.jcall(molly,"V","supplyNbLev",.jarray(flevs))
		} else .jcall(molly,"V","supplyNbLev",.jarray(integer(2)))

		# informs it about the number of df consumed by the error distrib
		options(warn=-1)
		.jcall(molly,"V","supplyErrorDF",as.integer(attr(logLik(fitfunc(as.formula(paste(y,"~1")),data=data, ...)),"df")-1))
		options(warn=0)
	}


	#
	   #
		  #


	if (method == "d") {
		# go for simple diagnostic !
		write("TASK: Diagnostic of candidate set.",file="")
		write(paste("Sample size:",ssize),file="")
		write(paste(nc,"factor(s)."),file="")
		write(paste(nq,"covariate(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxc"),"f exclusion(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxq"),"c exclusion(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxcxc"),"f:f exclusion(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxqxq"),"c:c exclusion(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxcxq"),"f:c exclusion(s)."),file="")
		write(paste("Size constraints: min = ",minsize,"max =",maxsize),file="")
		write(paste("Complexity constraints: min = ",minK,"max =",maxK),file="")
		if (marginality) 
			write("Marginality rule.",file="")
		nbcand = .jcall(molly,"I","diagnose")
		if (nbcand==-1) 		
			write("Your candidate set contains more than 1 billion (1e9) models.",file="")
		else 
			write(paste("Your candidate set contains",nbcand,"models."),file="")

	return(nbcand)
	}


	#
	   #
		  #



	if (method=="h") {
	# go for exhaustive screening 
	if (report) write("TASK: Exhaustive screening of candidate set.",file="")
	# init output object
	resulto <- new("glmulti")
	if (chunks>1) {
		resulto@name = paste(name,"_",chunk,".",chunks,sep="")
	} else 
		resulto@name = name
	resulto@params = list(name=name, intercept=intercept, marginality=marginality ,bunch=bunch, chunk=chunk, chunks=chunks, 
			level=level, minsize=minsize, maxsize=maxsize, minK=minK, maxK=maxK, method=method,crit=crit, confsetsize=confsetsize, fitfunction=fitfunction)
	resulto@call=match.call()
	resulto@adi=list(...)

	if (report) write("Fitting...",file="")
	flush.console()
	if(!.jcall(molly,"Z","produceModels",as.integer(chunk),as.integer(chunks), as.integer(bunch))) {
		warning("!Failed to start Java thread.")
		return(-1)
	}
	lesCrit<-numeric(confsetsize)
	lesK<-vector('integer',confsetsize)
	lesForms<-vector('character',confsetsize)
	
	lesObjects = list()

	# proceed with fitting
	curr<-0
	sel<-0

	while(.jcall(molly,"Z","nextModel")) {
		formula<-.jcall(molly,"[S","getCurrModel")

		beber<-lapply(formula, function(x) if (!is.na(x)) fitfunc(as.formula(x), data=data, ...) else NA)
		for (momo in 1:bunch) {
		curr = curr+1
		# convergence ?
		proceed=!is.na(formula[momo])
		if (proceed && fitfunction=="glm" && !beber[[momo]]$converged) {
			proceed<-FALSE
			warning(paste("!fitting function failed to converge. Model skipped: ", formula, sep=""))
		}
		if(proceed) {
			cricri<-support(beber[[momo]])
			if (sel<confsetsize) {
				sel=sel+1
				lesForms[sel]=formula[momo]
				lesCrit[sel]=cricri
				lesK[sel]=attr(logLik(beber[[momo]]),"df")
				lesObjects=c(lesObjects, list(beber[[momo]]))
			} else {
				mini=max(lesCrit)
				if (cricri<mini) {
					ou=which(lesCrit==mini)[1]
					lesForms[ou]=formula[momo]
					lesCrit[ou]=cricri
					lesK[ou]=attr(logLik(beber[[momo]]),"df")
					lesObjects[[ou]] = beber[[momo]]
				}
			}
		}
		if (curr%%(DELTAD) == 0) {
				if (report) {
					write(paste("\nAfter ",curr," models:",sep=""),file="")
					bestofsex=min(lesCrit[1:sel])
					write(paste("Best model: ",gsub(" ","",lesForms[which(lesCrit == bestofsex)]),sep=""),file="")
					write(paste("Crit= ",bestofsex,sep=""),file="")
					write(paste("Mean crit= ",mean(lesCrit[1:sel]),sep=""),file="")
					flush.console()
				} 
				if (plotty) {
					plot(sort(lesCrit[1:sel]),xlab="Best models",ylab=paste("Support (",crit,")"),pch=19,main="IC profile")
					abline(h=bestofsex+2,col="red")
				}
			}

		}
		
	}
	# over !
	if (report) write("Completed.",file="")
	# prepare and return glmulti object
	reglo <- order(lesCrit[1:sel])
	resulto@crits = lesCrit[1:sel][reglo]
	resulto@formulas = lapply(lesForms[1:sel][reglo],as.formula)
	resulto@K = as.integer(lesK[1:sel][reglo])
	resulto@nbmods = as.integer(sel)
	if (includeobjects) resulto@objects = lesObjects[reglo] else resulto@objects = list()
	return(resulto)
	
	#
	   #
		  #
		  #
		#
	#



	} else if (method=="l") {
	# go for exhaustive screening using leaps
	if (report) write("TASK: Exhaustive screening of candidate set, branch-and-bound algorithm.",file="")
	write("[ Be sure to have package leaps installed ]", file="")
	if (level==2 || nc>0) stop("Method l cannot be used with factors or interactions")
	fitfunction="lm"
	# init output object
	resulto <- new("glmulti")
	resulto@name = name
	resulto@params = list(name=name, intercept=intercept, marginality=marginality ,chunk=chunk, chunks=chunks, 
			level=level, minsize=minsize, maxsize=maxsize, minK=minK, maxK=maxK, method=method,crit=crit, confsetsize=confsetsize, fitfunction=fitfunction)
	resulto@call=match.call()
	resulto@adi=list(...)
	if (report) write("Fitting...",file="")
	flush.console()
	# proceed with leaps model selection
	f <- as.formula(paste(y, paste(xq[-1], collapse="+"),sep="~"))
	try(regsubsets(f, data = databis, nbest = confsetsize, really.big=T, intercept=intercept), silent=F)-> lilo				
	if (class(lilo)== "try-error") 
					stop("!call to regsubsets failed.")
	summary(lilo, matrix=T, matrix.logical=T)-> mama
	mama[[7]] -> coda
	mama[[3]] -> lesrss
	nbmods= length(coda[,1])
	as.numeric(sapply(dimnames(coda)[[1]], function(x) strsplit(x, split=" ")[[1]][1]))-> KK
	pena = 0
	if (intercept) pena=1
	length(resid(lm(as.formula(paste(y,"1",sep="~")),data=databis)))-> sssize
	if (crit=="aicc") 
	as.numeric(lapply(1:nbmods, function(x) sssize*log(lesrss[x]/sssize) + 2*(KK[x]+1+pena)*sssize/max(0,sssize-KK[x]-2-pena)))-> lesic
	else if (crit=="aic")
	as.numeric(lapply(1:nbmods, function(x) sssize*log(lesrss[x]/sssize) + 2*(KK[x]+1+pena)))-> lesic
	else as.numeric(lapply(1:nbmods, function(x) sssize*log(lesrss[x]/sssize) + log(sssize)*(KK[x]+1+pena)))-> lesic
				
	# we want to include the null model if appropriate
	if (intercept) {
	sum(lm(as.formula(paste(y, "1",sep="~")), data = databis)$residuals^2)-> nullos
	if (crit=="aicc")  sssize*log(nullos/sssize) + 2*(1+pena)*sssize/max(0,sssize-2-pena)-> icnull
	else if (crit=="aic") sssize*log(nullos/sssize) + 2*(1+pena)-> icnull
	else sssize*log(nullos/sssize) + log(sssize)*(1+pena)-> icnull
	KK=c(KK,0)
	lesic=c(lesic, icnull)
	coda=rbind(coda, rep(F, length(coda[1,])))
	nbmods=nbmods+1
	}
	if (report) write("Completed.",file="")
	
	# extract the best models
	order(lesic)-> oo
	lesic[oo]-> lesic
	KK[oo]-> KK
	coda[oo,]-> coda	
	nbbyk = tapply(KK, KK, length)
	worro = tapply(lesic, KK, max)
	if (max(nbbyk)==confsetsize) thresh = min(worro[nbbyk == confsetsize]) else thresh=max(worro)
	KK[lesic<=thresh] -> KK
	coda[lesic<=thresh,]-> coda
	lesic = lesic[lesic<=thresh]
	nbmodels = length(KK)
	if (report) write(paste(nbmodels, "first best models identified."),file="")
		
	#
	# pack the object
	resulto@crits = lesic
	resulto@K  = as.integer(KK+1+pena)
	# prepare formulas
	apply(coda, 1, function(x) {
			if (intercept) as.formula(paste(y, paste(c("1", xq[-1][x]), collapse="+"),sep="~"))
			else as.formula(paste(y, paste(c("-1", xq[-1][x]), collapse="+"),sep="~"))
	}) -> formo
	resulto@formulas = formo			
	resulto@nbmods = as.integer(nbmodels)
	# if objects are requested to be included, refit models
	if (includeobjects) {
	fitfunc = match.fun(fitfunction)
	lesObjects= lapply(formo, function(x) fitfunc(x, data=data, ...))
	resulto@objects = lesObjects
	}	

	return(resulto)
	
	#
	   #
		  #
		  #
		#
	#


	} else {
		#Go for genetic algorithm approach
		write("TASK: Genetic algorithm in the candidate set.",file="")
		# init output object
		resulto <- new("glmulti")
		resulto@name = name
		resulto@params = list(name=name, intercept=intercept, marginality=marginality ,bnch=bunch, chunk=chunk, chunks=chunks, 
				level=level, minsize=minsize, maxsize=maxsize, minK=minK, maxK=maxK, method=method,crit=crit, confsetsize=confsetsize, fitfunction=fitfunction, popsize=popsize,mutrate=mutrate,
				sexrate=sexrate,imm=imm, plotty=plotty, deltaM=deltaM, deltaB=deltaB, conseq=conseq, resumefile = resumefile)
		resulto@call=match.call()
		resulto@adi=list(...)
		lesObjects=list()
		
		write("Initialization...",file="")
		currgen = 0
		consoude = 0
		gogo=TRUE
		bestofsex = 10000
		minou = 10000
		bestofsexN = 1000
		minouN = 1000
		if (method=="r") 
			popul = .jcall(molly,"[S","initPopAgain", mutrate, imm, sexrate, name)
		else popul = .jcall(molly,"[S","initPop",as.integer(popsize), mutrate, imm, sexrate, as.integer(confsetsize), name)

		tini = Sys.time()
		dyniT = numeric(0)
		dyniB = numeric(0)
		dyniM = numeric(0)
		write("Algorithm started...",file="")
		while (gogo) {
			for (i in 1:DELTAD) { 
				nbtofit = length(popul)
				lesic = numeric(nbtofit)
				if (nbtofit>0) 
					for (m in 1:nbtofit) {
						# fit models
						formula=popul[m]
						beber<- fitfunc(as.formula(formula), data=data, ...)
						liliac<- logLik(beber)
						K<-attr(liliac,"df")
						# convergence ?
						if (fitfunction=="glm" && !beber$converged) 
							lesic[m]=10000
						else lesic[m]=support(beber)
						}
				currgen = currgen+1
				popul = .jcall(molly,"[S","nextGeneration",.jarray(lesic))
			}
			# report about the current state and store things
			lesCrit = .jcall(molly,"[D","reportConfIC")
			bestofsexN=.jcall(molly,"D","reportBestIC")
			minouN=.jcall(molly,"D","reportMeanIC")
			bestform = .jcall(molly,"S","reportbestModel")
			dyniT = c(dyniT,as.numeric(Sys.time()-tini))
			dyniB = c(dyniB,bestofsexN)
			dyniM = c(dyniM,minouN)
			if (report) {
				write(paste("\nAfter ",currgen," generations:",sep=""),file="")
				write(paste("Best model: ",gsub(" ","",bestform),sep=""),file="")
				write(paste("Crit= ",bestofsexN,sep=""),file="")
				write(paste("Mean crit= ",minouN,sep=""),file="")
				flush.console()
				}
			.jcall(molly,"V","printImage")
			if (plotty) {
				plot(sort(lesCrit),xlab="Best models",ylab=paste("Support (",crit,")"),pch=19,main="IC profile")
				abline(h=bestofsexN+2,col="red")
			}
			# Shall we continue ?
			if (length(lesCrit)==confsetsize && minouN - minou >= -deltaM && bestofsexN - bestofsex >= -deltaB) 
				consoude = consoude+1
			else consoude = 0
			if (consoude == conseq) {
				write("Improvements in best and average IC have bebingo en below the specified goals.",file="")
				write("Algorithm is declared to have converged.",file="")
				gogo=FALSE
			} else {
				if (report) write(paste("Change in best IC:",bestofsexN - bestofsex, "/ Change in mean IC:",minouN - minou),file="")
				minou = minouN
				bestofsex = bestofsexN
			} 
			}
		# END OF GA
		write("Completed.",file="")
		lesForms= .jcall(molly,"[S","reportConfMods")
		lesKK = .jcall(molly,"[I","reportConfKs")
		sel = length(lesCrit)
		# prepare and return glmulti object
		reglo <- order(lesCrit)
		resulto@crits = lesCrit[reglo]
		resulto@formulas = lapply(lesForms[reglo],as.formula)
		resulto@K = as.integer(lesKK)[reglo]
		resulto@nbmods = as.integer(sel)
		# if objects are requested to be included, refit models
		if (includeobjects) {
			lesObjects= lapply(lesForms[reglo], function(x) if (!is.na(x)) fitfunc(as.formula(x), data=data, ...) else NA)
			resulto@objects = lesObjects
		}	
		resulto@params = c(resulto@params, list(generations=currgen, elapsed=as.numeric(Sys.time()-tini), dynat=dyniT, dynab=dyniB, dynam=dyniM))
		return(resulto)

	}


	#
	   #
		  #
		  

	# end of function.
})

