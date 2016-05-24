############################
flip.npc.methods <-
    c("Fisher", "Liptak", "Tippett", "MahalanobisT", "MahalanobisP", "minP", "maxT", "maxTstd", "sumT", "Direct", "sumTstd", "sumT2", "kfwer", "data.sum","data.linComb","data.pc","data.trace")
############################

npc <- function(permTP, comb.funct = c(flip.npc.methods, p.adjust.methods) ,subsets=NULL,weights=NULL, stdSpace=FALSE, ...){
#	on.exit(browser())
	### just in analogy with gt(). to be implemented as flip-options
	trace=TRUE
	
	#### arrange comb.funct
	comb.funct <- match.arg(tolower(comb.funct[1]),tolower(flip.npc.methods))
	comb.funct <- flip.npc.methods[which(tolower(flip.npc.methods)==comb.funct)]
	if(comb.funct == "Tippett") comb.funct ="minP"
	if(comb.funct == "Direct") comb.funct ="sumT"
	
	if(is(permTP,"flip.object")){
    if(!is.null(list(...)$tail)) {  permTP@tail=tail  }
		nperms = permTP@call$B
    
    ##########flipMix type of combinations
    if(comb.funct %in% c("data.sum","data.pc","data.linComb","data.trace")) {
      if(comb.funct %in% c("data.trace")) test <- .trace.between.nptest
      if(comb.funct %in% c("data.sum","data.linComb","data.pc")) {
        if(!is.null(list(...)$statTest)&&(list(...)$statTest == "F"))
        test <-.F.between.nptest else
          test <-.t.between.nptest         
      }
      environment(test) <- permTP@call.env
      environment(test)$otherParams=list(...)
      environment(test)$otherParams$subsets=subsets
      if(!is.null(environment(test)$otherParams$perms))   environment(test)$perms=environment(test)$otherParams$perms
      if(comb.funct %in% c("data.sum","data.linComb","data.pc")) environment(test)$otherParams$onlyMANOVA=TRUE
      if(comb.funct %in% c("data.sum")) environment(test)$otherParams$linComb=1
      res=test()
#      browser()
      if(!exists("flipReturn") || is.null(flipReturn)) 
        flipReturn=list(permT=TRUE,permP=FALSE,call.env=TRUE)
      
      out=.getOut(type="npc",res=list(permT=res$permT,extraInfoPre=list(comb.funct=comb.funct,nVar=res$extraInfoPre$nVar)),data=NULL,tail=list(...)$tail, call=match.call(), 
                flipReturn=flipReturn,call.env=environment(test))
      return(out)
      
    } else	
      if(comb.funct %in% c("Fisher", "Liptak", "minP")) {
			if(!is.null(permTP@permP)){ 
				permTP=permTP@permP			
			} else { 
				  if(!is.null(permTP@permT)) { 
            permTP=t2p(.fixPermT(permTP@permT),obs.only=FALSE,tail=permTP@tail)} else 
				    {print("Joint distribution of p-values not provided. Nothing done."); return()}
			} 
		} else 
      if(comb.funct %in% c("maxT", "maxTstd", "sumT","sumTstd", 
                           "sumT2","MahalanobisT","MahalanobisP")) {
			if(is.null(permTP@permT)) {print("Joint distribution of p-values not provided. Nothing done."); return()}
			if(comb.funct %in% c("MahalanobisT","MahalanobisP")) permTP =permTP@permT else
         permTP=.setTail(.fixPermT(permTP@permT),tail=.fitTail(permTP@permT,permTP@tail) ) 
        if(comb.funct %in% c("maxTstd","sumTstd"))
           permTP=scale(permTP)
            permTP
      }
	} else if(!is.null(list(...)$tail)) {  
	  permTP=.setTail(.fixPermT(permTP),tail=tail)
	}
	
	
	#if(!exists("nperm"))  nperms = list(number=nrow(permTP),seed=NA)
	if(!is.matrix(permTP)) permTP=as.matrix(permTP)
	if(stdSpace & (comb.funct %in% c("maxT", "sumT", "sumT2"))) {permTP = .t2stdt(permTP,FALSE)}
	if(comb.funct=="Fisher"){permTP = -log(permTP)} else
	if(comb.funct=="Liptak"){permTP = -qnorm(permTP)} else
	if(comb.funct=="minP"){permTP = 1/permTP} else
	  if(comb.funct=="sumT2"){permTP = permTP^2}	
	
	
	############
	temp=.getSubsetWeights(weights,subsets,colnames(permTP))
	weights=temp$weights
	subsets=temp$subsets
	many.weights=temp$many.weights
	many.subsets=temp$many.subsets
	one.weight=temp$one.weight
	
	rm(temp)
	
	# prepare progress info
   ###  if (missing(trace)) trace <- gt.options()$trace && (many.weights || many.subsets)
  if (trace && (many.subsets || many.weights)) {
    if (many.subsets) 
      K <- length(subsets) 
    else 
      K <- length(weights) 
    digitsK <- trunc(log10(K))+1
  }

  
  # weight
    if (one.weight) {
      if (length(weights) != ncol(permTP)) 
        stop("length of \"weights\" does not match column count of \"permTP\"")
      all.weights <- weights
    } else {
      if(comb.funct%in%c("sumT","sumTstd")) 
		all.weights <- rep(1/sqrt(ncol(permTP)), ncol(permTP)) else 
		all.weights <- rep(1, ncol(permTP)) 	 
	}
	names(all.weights) = colnames(permTP)
  
	if(comb.funct %in% c("Fisher", "Liptak", "sumT", "sumT2", "sumTstd"))
		  test= function(subset=NULL,weights=NULL){ 
		  permT = matrix(if(is.null(subset)) permTP%*%all.weights else permTP[,subset,drop=FALSE]%*%all.weights[subset]) ;
      permT} else 	
  if(comb.funct %in% c("minP", "maxT", "maxTstd"))
		   test= function(subset=NULL,weights=NULL){ #browser()
					permT = matrix(apply(if(is.null(subset)) { if(one.weight) t(all.weights*t(permTP)) else permTP } else 
					t(all.weights[subset]*t(permTP[,subset,drop=FALSE])) , 1, max))  ; 
					permT
          } else 
  if(comb.funct %in% c("MahalanobisT","MahalanobisP")) 
        test= function(subset=NULL,weights=NULL){
         if(is.null(subset)) 
           if(one.weight) pseudoT<-permTP%*%diag(all.weights) else 
             pseudoT<-permTP  else
             if(one.weight) pseudoT<-permTP[,subset,drop=FALSE]%*%diag(all.weights[subset])  else 
               pseudoT<-permTP[,subset,drop=FALSE]
         if(ncol(pseudoT)==1)
           return(scale(pseudoT)^2)
           
         if(comb.funct %in% "MahalanobisP")          
           pseudoT=qnorm( t2p(pseudoT,tail=1,obs.only=FALSE)*.99999999999)  else
             pseudoT=scale(pseudoT,scale=FALSE)

         ei=eigen(t(pseudoT)%*%pseudoT,symmetric=TRUE)
         ei$values[ei$values<0]=0
         invhalfEV=ei$values^-.5
         keep= is.finite(invhalfEV)
					  permT=matrix(rowSums((pseudoT %*% ei$vectors[,keep] %*% diag(invhalfEV[keep]))^2)*sum(keep))
         permT
					}
            

					
	# Do the test
  if ((!many.subsets) && (!many.weights)) {           # single weighting; single subset
    permT <- test()
	nVar=ncol(permTP)
  } else {     
    L <- if (many.subsets) length(subsets) else length(weights)
    permT <- sapply(1:L, function (i) { 
      if (trace && L>1) {
        cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
        flush.console()
      }
      if (!many.weights) {                                           # single weighting; many subsets
        uit <- test(subset=subsets[[i]]) 
      } else if (!many.subsets) {                                    # many weightings; single subset
        uit <- test(weights=weights[[i]]) 
      } else {                                                      # many weightings; many subsets
        uit <- test(subset=subsets[[i]], weights=weights[[i]])
      } 
	  uit
    })
    if (many.subsets && !is.null(names(subsets))){
      colnames(permT) <- names(subsets)
	}
    else if (many.weights && !is.null(names(weights))){
      colnames(permT) <- names(weights)
	}
  nVar={if(is.null(subsets)) 1 else sapply(subsets,length)}*ifelse(is.null(weights), 1, sapply(weights,length))
  }
  if (trace && (many.subsets || many.weights) && L>1) cat("\n")
	
	
	
	  if(!exists("flipReturn") || is.null(flipReturn)) 
			flipReturn=list(permT=TRUE,permP=FALSE)
  #build the flip-object
	out=.getOut(type="npc",res=list(permT=permT,extraInfoPre=list(comb.funct=comb.funct,nVar=nVar)),data=list(),tail=list(...)$tail, call=match.call(), flipReturn=flipReturn)
	return(out)
}