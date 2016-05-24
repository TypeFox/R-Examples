# TODO: Add comment
# 
# Author: nmuller
###############################################################################


seqerulesdisc <- function(fsubseq, datadiscr, tsef, pvalue=0.1, supvars=NULL, adjust=TRUE, topt=FALSE, link="cloglog", dep=NULL) {
	## supvars is the strata
	nvars <- length(datadiscr[["data"]])
	res <- list("rule"=list(), "discrule"=list(), "model"=list(), "support"=list(), "confidence"=list(), "lift"=list(), "slift"=list())
	eventshr <- list("events"=list(), "hr"=list(), "hrindep"=list(), "pvalues"=list())
	nrule <- 1
	## Total number of sequences
	ntotal <- length(fsubseq$seqe)
	
	
	if(topt==TRUE) {
		subseqs.tmp <- as.character(fsubseq$subseq)
		subseqs <- NULL
		matcountr <- seqeapplysub(fsubseq, method="presence", rules=TRUE)
		
		
		for(i in 1:length(subseqs.tmp)) {
			if(sum(matcountr[,as.character(subseqs.tmp[i])])==1) {
				subseqs <- c(subseqs, subseqs.tmp[i])
			}
		}
	}
	else {
		subseqs <- fsubseq$subseq
	}
	
	if(adjust) {
		pvalue <- pvalue/length(subseqs)
		message("Adjusted p-value threshold : ", pvalue)
	}
	
	for(i in 1:length(subseqs)) {
		sseq <- subseqs[i]
		splitted2 <- strsplit(as.character(sseq), "-")[[1]]
		sseq <- gsub(",", "\\)-\\(", sseq)
		splitted <- strsplit(as.character(sseq), "-")[[1]]
		if((!is.null(dep) && splitted[length(splitted)]==paste("(",dep, ")", sep="")) | is.null(dep)) {
			
			
			hrs <- NULL
			hrindep <- NULL
			pvalues <- NULL
			discrule <- NULL 
			## subsequence must have more than one element
			if(length(splitted)>1) {
				left <- sapply(splitted, function(x) gsub("[\\(\\)]", "", x), USE.NAMES=FALSE)
				for(j in 2:(length(left))) {
					# sélection des données correspondant à la partie 
					data.cur <- datadiscr[["data"]][[left[j]]]
					## when we arrive at the third variable, we control, from the first variable  
					if(j>2) {
						for(k in 1:(j-2)) {
							## selecting ids with events order (a -> b -> c, then a$age < b$age)
							iddep <- tsef[(tsef[,left[k]] < tsef[,left[k+1]]) & tsef[,paste(left[k], "ST", sep="")]==1,1]
							
							data.cur <- data.cur[data.cur$IDPERS%in%iddep,]
							data.cur <- subset(data.cur,get(left[k])==1)
						}
					}
					## we keep as the independent variable the preceding event
					if(!is.null(supvars)) {
						data.varindeps <- subset(data.cur, ,c("AGE", left[j-1]))
					}
					else {
          
						  data.varindeps <- subset(data.cur, ,c(left[j-1]))
           
          }
					## if there is a list two individuals who experienced the event
					if(sum(data.cur[,3]) > 1) {
						if(!is.null(supvars)) {
							try(coxmodel <- coxph(Surv(data.cur[,2], data.cur[,3], data.cur[,(nvars+3)]) ~ data.varindeps[,left[j-1]] + strata(data.varindeps[,supvars])))
						}
						else {
							model <- glm(data.cur[,3] ~ data.cur[,left[j-1]] + data.cur[,"AGE"], family=binomial(link=link))
						}
						z <- summary(model)
						pv <- z$coefficient[2,4]
						coef <- z$coefficient[2,1]
						pvalues<-c(pvalues, pv)
						hrs <- c(hrs, exp(coef))
						
						discrule <- paste(discrule, left[j-1], "(HR=", round(exp(coef), digits=2), ",p=",round(pv,digits=4),") =>", sep="")
						
						
					}
					else { discrule <- "No occurences of " }
				}
				discrule <- paste(discrule, left[length(left)])
				anyp <- !any(pvalues>pvalue)				
				## Saving results in the list
				if(!is.na(anyp) && anyp && !(discrule%in%res[["discrule"]])) {
					res[["discrule"]][[nrule]] <- discrule
					eventshr[["events"]][[nrule]] <- left
					eventshr[["hr"]][[nrule]] <- hrs
					eventshr[["pvalues"]][[nrule]] <- pvalues
					if(!is.null(supvars)) {
						eventshr[["hrindep"]][[nrule]] <- hrindep
					}
					nrule<-nrule+1
				}
				## end of if(length(splitted)>1) {
			}
		}
		## end of main 'for'  
	}
	
	results.frame <- data.frame(
			"Hazards"=unlist(res$discrule)
	)  
	class(results.frame) <- c("rules", "data.frame")
	attr(results.frame, "graph") <- eventshr
	print("END")
	return(results.frame)       
}

