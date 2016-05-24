fitLMER.fnc <- function(model,
    item=FALSE, # can be an item identifier such as "Item" or "Word"
	backfit.on=c("F","t"),
	method=c("F","t","z","llrt","AIC", "BIC","relLik.AIC","relLik.BIC"),
	threshold=NULL,
	t.threshold=NULL,
    ran.effects=list(ran.intercepts=as.character(),
             slopes=as.character(),
		     corr=as.character(),
             by.vars=as.character()),
    alpha=NULL,
	alphaitem=NULL,	
    if.warn.not.add=TRUE,
    prune.ranefs=TRUE,
    p.value="upper",
    set.REML.FALSE=TRUE,
	keep.single.factors=FALSE,	
    reset.REML.TRUE=TRUE,
    log.file.name=NULL # or other path and file name or FALSE
){
	if(backfit.on[1]=="F"&&method[1]=="t"){
		warning("resetting argument \"method\" to \"F\"\n")
		method<-"F"
	}

	if(backfit.on[1]=="t"&&method[1]=="F"){
		warning("resetting argument \"method\" to \"t\"\n")
		method<-"t"
	}

	if(as.vector(model@call[1])=="glmer()"){
		if(backfit.on=="F"){
			backfit.on<-"t"
		}
	}

	if(is.null(log.file.name)){
		log.file.name<-file.path(tempdir(),paste("fitLMER_log_",gsub(":","-",
			gsub(" ","_",date())),".txt",sep=""))
	}

  	current.dir=getwd()
  	temp.dir=tempdir()
  	tempdir()

  	if(log.file.name!=FALSE)sink(file=log.file.name,split=TRUE)  
  	cat("======================================================\n")
  	cat("===              backfitting fixed effects         ===\n")
  	cat("======================================================\n")
  	if(backfit.on[1]=="F"){
		mod=bfFixefLMER_F.fnc(model=model,
			item=item,method=method,threshold=threshold,alpha=alpha,
			alphaitem=alphaitem,prune.ranefs=prune.ranefs,
			p.value=p.value,set.REML.FALSE=set.REML.FALSE,
			keep.single.factors=keep.single.factors,reset.REML.TRUE=FALSE,
			log.file=FALSE) 
  	}else{
    	mod=bfFixefLMER_t.fnc(model=model,item=item,method=method,
			threshold=threshold,t.threshold=t.threshold,alphaitem=alphaitem,
			prune.ranefs=prune.ranefs,set.REML.FALSE=set.REML.FALSE,
			keep.single.factors=keep.single.factors,
			reset.REML.TRUE=FALSE,log.file=FALSE)
  	}


	cat("======================================================\n")
	cat("===            forwardfitting random effects       ===\n")
	cat("======================================================\n")
	mod=ffRanefLMER.fnc(model=mod,ran.effects=ran.effects,
		alpha=ifelse(is.null(alpha),0.05,alpha),
		if.warn.not.add=if.warn.not.add,log.file=FALSE)


  	cat("======================================================\n")
  	cat("===            re-backfitting fixed effects        ===\n")
  	cat("======================================================\n")
  	if(backfit.on[1]=="F"){
		mod=bfFixefLMER_F.fnc(model=mod,
			item=item,method=method,threshold=threshold,alpha=alpha,
			alphaitem=alphaitem,prune.ranefs=prune.ranefs,
			p.value=p.value,set.REML.FALSE=set.REML.FALSE,
			keep.single.factors=keep.single.factors,
			reset.REML.TRUE=reset.REML.TRUE,log.file=FALSE) 
  	}else{
		mod=bfFixefLMER_t.fnc(model=mod,item=item,method=method,
			threshold=threshold,t.threshold=t.threshold,alphaitem=alphaitem,
			prune.ranefs=prune.ranefs,set.REML.FALSE=set.REML.FALSE,
			keep.single.factors=keep.single.factors,
			reset.REML.TRUE=reset.REML.TRUE,log.file=FALSE)
  	}

  	options(warn=0)
  	sink(file=NULL,type="message")        

  	if(log.file.name!=FALSE){
    	cat("log file is",log.file.name,"\n")
    	sink(file=NULL)
  	}

  	return(model=mod)
}

