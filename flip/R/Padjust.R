flip.adjust <- function (permTP, method = flip.npc.methods, maxalpha=1, weights=NULL, stdSpace=FALSE, ...) {
    
	##TODO: here is case sensitive, while in npc it is not. here is so for compatibility with p.adjust. shall we change?
  if(length(method)>1) method <- "maxT" else 
    method <- match.arg(method,c(p.adjust.methods,flip.npc.methods))
  
	
	##TODO: add some check for names(weights) compatibility. 
	##TODO: exclude missing names (eg permTP=permTP[names(weights)]
	if(!is.null(weights)) if(is.null(names(weights))) names(weights)<- colnames(permTP)
	
	if(is(permTP,"flip.object")){
	  if(is.null(list(...)$tail)) tail=NULL else permTP@tail=list(...)$tail
		if((method%in%p.adjust.methods) & (is.null(weights))){ 
			#run the standard function
			adjs=p.adjust(p.value(permTP), method = method)
		} else if((method%in%p.adjust.methods) & (!is.null(weights))){
			#standard function but weigthed
			if(method=="hommel") {print("weighted \"hommel\" method not allowed."); return()}
			adjs=p.adjust.w(p.value(permTP), method = method, w=weights)
		} else {		#perform permutation specific procedures
			if(method=="minP") {
				adjs=.maxt.adjust(if(is.null(permTP@permP)) -t2p(permTP@permT,obs.only=FALSE,tail=permTP@tail) else -permTP@permP,  maxalpha,weights=weights)
			} else if(method=="maxT") {
				if(stdSpace) {permTP@permT = .t2stdt(permTP@permT,FALSE)}
				adjs=.maxt.adjust(.setTail(.fixPermT(permTP@permT), tail=permTP@tail), maxalpha,weights=weights)
			} else if(method=="kfwer") {
			  if(is.null(list(...)$k)) k=11 else 
          k=list(...)$k
			  if(stdSpace) {permTP = .t2stdt(permTP@permT,FALSE)}
			  adjs=.p.adjust.kfwer(.setTail(.fixPermT(permTP@permT), tail=tail), k=k)
			} else {
			#otherwise perform closed testing
			if(stdSpace & (method %in% c("sumT", "sumT2"))) {permTP@permT = .t2stdt(permTP@permT,FALSE)}
			if(is.null(permTP@permP) & (method %in% c("Fisher", "Liptak", "minP"))) {
			  permTP@permP=t2p(.fixPermT(permTP@permT),obs.only=FALSE,tail=permTP@tail)
			}
      
      # Define the local test to be used in the closed testing procedure
      mytest <- function(hyps) {p.value(npc(permTP[hyps],method,weights=weights[hyps]))}
      cl <- closed(mytest, names(permTP),alpha=NA)
      adjs=sapply(names(permTP),function(id) adjusted(cl,id))
			}
		}
		#fit it to the flip-object
		permTP@res=cbind(permTP@res,adjs)
    pastemethod <- function(){
      if(method=="kfwer")
        method=paste("k",k,"fwer",sep="")
      method
    }
	  if(paste("Adjust:",pastemethod(),sep="")%in%colnames(permTP@res)){
      i=2
      while(paste("Adjust:",pastemethod(),".",i,sep="")%in%colnames(permTP@res)) i=i+1
        newname=paste("Adjust:",pastemethod(),".",i,sep="")
	  } else newname=paste("Adjust:",pastemethod(),sep="")
    colnames(permTP@res)[length(colnames(permTP@res))]=newname
		return(permTP)
			
	} else { # not a flip.object, return a vector
    if(is.null(list(...)$tail)) tail=NULL else tail=list(...)$tail
    
		if((method%in%p.adjust.methods) & (is.null(weights))){ 
			return(p.adjust(permTP, method = method))
		} else if((method%in%p.adjust.methods) & (!is.null(weights))){
			#standard function but weigthed
			if(method=="hommel") {print("weighted \"hommel\" method not allowed."); return()}
			return(p.adjust.w(permTP, method = method, w=weights))
		} else{		#perform permutation specific procedures
			if(method=="minP") {
				return(.maxt.adjust(-t2p(permTP,obs.only=FALSE,tail=tail), maxalpha,weights=weights))
			} else if(method=="maxT") {
			    if(stdSpace) {permTP = .t2stdt(permTP,FALSE)}
			  return(.maxt.adjust(.setTail(.fixPermT(permTP), tail=tail), maxalpha,weights=weights))
			} else if(method=="kfwer") {
        if(is.null(k)) k=11
			  if(stdSpace) {permTP = .t2stdt(permTP,FALSE)}
			  return(.p.adjust.kfwer(.setTail(.fixPermT(permTP), tail=tail), k=k))
			} else {
			#then perform closed testing
			# Define the local test to be used in the closed testing procedure
			mytest <- function(hyps) {p.value(npc(permTP[,hyps],method,weights=weights[hyps]))}
			cl <- closed(mytest, colnames(permTP),adjust=TRUE)
			adjs=cl@adjusted[1:ncol(permTP)]
			}
		}
	}
	adjs
}


################
.maxt.adjust <- function(permT,maxalpha=1,weights=NULL,m=ncol(permT)) {
  
	#get colnames to 
	if(is.null(colnames(permT))) colnames(permT)=1:m
	
	#define the order of testing hypotheses
	steps=names(sort(permT[1,],decreasing = TRUE))
	
	if(!is.null(weights)) permT=t(weights*t(permT))
	
	#set of NOT yet rejected p-values
	notrejs=rep(TRUE,m)
	names(notrejs)=colnames(permT)
	
	#set of adjusted p-values
	Padjs=rep(1,m)
	names(Padjs)=colnames(permT)
	
	i <- 1
	while((i<=m) & ifelse(i>1,Padjs[steps[i-1]] <= maxalpha,TRUE)){
		Padjs[steps[i]]=max( 
      t2p(c(permT[1,steps[i]], apply(permT[-1,notrejs,drop=FALSE],1,max) )) ,
      Padjs[steps[i-1]] ) #first max ensures monotonicity
		notrejs[steps[i]]=FALSE
		
		#avoid to compute max if test statistic are equal (specially useful in minp)
		while(permT[1,steps[i]]==ifelse(i==m,Inf,permT[1,steps[i+1]])){  
			Padjs[steps[i+1]]=Padjs[steps[i]]
			notrejs[steps[i+1]]=FALSE
			i=i+1
		}
		
		i=i+1
	}

	return(Padjs)
}
