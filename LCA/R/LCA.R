LCA <- function(x,PTLmodel,clique,seed.row,combine.method="Fisher",adjust.method="BH",comparison.alpha=0.05){

	D <- as.matrix(dist(t(x)))

	cat("evaluating comparison probabilities \n")	
	all.combinations <- combn(clique,m=2)
	comb.ps <- mapply(function(a,b){evaluateDiffSignificance(d=D[a,b],diff=x[seed.row,b]-x[seed.row,a],PTLmodel=PTLmodel)},a=all.combinations[1,],b=all.combinations[2,])
	comb.ps <- p.adjust(comb.ps,method=adjust.method)

	sig.combinations <- all.combinations[,which(comb.ps<comparison.alpha)]
	if(sum(comb.ps<comparison.alpha)==0) warning("no significant comparisons")
	out.frame <- NULL
	rownames(out.frame) <- rownames(x)

	if(sum(comb.ps<comparison.alpha)==1){
		sig.combinations <- as.matrix(sig.combinations)
	}

	if(sum(comb.ps<comparison.alpha)>1){
		out.frame <- data.frame(p.value=rep(NA,nrow(x)),adj.P.value=rep(NA,nrow(x)))
		cat(paste("evaluating LCD estimates for",sum(comb.ps<comparison.alpha),"comparisons \n"))
		sig.comb.ps <- comb.ps[which(comb.ps<comparison.alpha)]
		LCD.ps <- mapply(function(a,b){LCD(x1=x[,a],x2=x[,b],seed.row=seed.row,PTLmodel=PTLmodel)},a=sig.combinations[1,],b=sig.combinations[2,])

		# combine p-values for each row (weighted combination on sig.comb.ps)
		cat(paste("combining",nrow(LCD.ps),"p-values for",ncol(LCD.ps),"comparisons \n"))
		if(combine.method=="inverseProd"){
			observation.p <- 1-pchisq(-2*sum(log(sig.comb.ps)),df=2*length(sig.comb.ps))
			invProdFun <- function(x,w){
                                logprobs <- log(1-x)
                                weighted_logprob_sum <- (sum(logprobs*w))/sum(w)
                                out <- exp(weighted_logprob_sum)
                                max(c(2.2e-16,1-out))
                        }
			out.pvals <- 1-((1-apply(LCD.ps,MARGIN=1,invProdFun,w=sig.comb.ps))*(1-observation.p))
		}
		if(combine.method=="Fisher"){
			out.pvals <- apply(LCD.ps,MARGIN=1,function(x,w){max(c(2.2e-16,1-pchisq(-2*sum(log(1-((1-x)*(1-w)))),df=2*length(x))))},w=sig.comb.ps)
		}
		
		adj.pvals <- p.adjust(out.pvals,method=adjust.method)
		# output row-wise p-value estimates
		out.frame <- data.frame(p.value=out.pvals,adj.P.value=adj.pvals)
		rownames(out.frame) <- rownames(x)
	}
	list(LCD=out.frame[order(out.frame$p.value,decreasing=FALSE),],comparisons=sig.combinations)
}


