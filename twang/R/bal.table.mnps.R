bal.table.mnps <- function(x, digits = 3, collapse.to = c("pair","covariate","stop.method")[1], subset.var = NULL, subset.treat = NULL, subset.stop.method = NULL, es.cutoff = 0, ks.cutoff = 0, p.cutoff = 1, ks.p.cutoff = 1, ...) {
   	tmt1 <- tmt2 <- stop.method <- std.eff.sz <- ks <- p <- ks.pval <- max.std.eff.sz <- max.ks <- min.p <- NULL
   	min.ks.pval <- control <- NULL
   	if(x$estimand == "ATE"){
   		pwc <- pairwiseComparison(x, collapse.to = collapse.to)
   		if(!is.null(subset.var) & !(collapse.to == "stop.method")) pwc <- subset(pwc, var %in% subset.var | var %in% paste(subset.var, "<NA>", sep = ":"))
   		if(!is.null(subset.treat)){ 
   			if(length(subset.treat) == 1) pwc <- subset(pwc, tmt1 == subset.treat | tmt2 == subset.treat)
   			if(length(subset.treat) > 1) pwc <- subset(pwc, tmt1 %in% subset.treat & tmt2 %in% subset.treat)   
   		}
   		if(!is.null(subset.stop.method)) pwc <- subset(pwc, stop.method %in% subset.stop.method)
   		if(collapse.to == "pair"){
   			pwc <- subset(pwc, abs(std.eff.sz) >= es.cutoff | is.na(std.eff.sz) | is.nan(std.eff.sz))
   			pwc <- subset(pwc, ks >= ks.cutoff | is.na(ks) | is.nan(ks))
   			pwc <- subset(pwc, p <= p.cutoff | is.na(p) | is.nan(p))
   			pwc <- subset(pwc, ks.pval <= ks.p.cutoff | is.na(ks.pval) | is.nan(ks.pval))   
   		}
   		else {
   			pwc <- subset(pwc, abs(max.std.eff.sz) >= es.cutoff | is.na(max.std.eff.sz) | is.nan(max.std.eff.sz))
   			pwc <- subset(pwc, max.ks >= ks.cutoff | is.na(max.ks) | is.nan(max.ks))
   			pwc <- subset(pwc, min.p <= p.cutoff | is.na(min.p) | is.nan(min.p))
   			pwc <- subset(pwc, min.ks.pval <= ks.p.cutoff | is.na(min.ks.pval) | is.nan(min.ks.pval))   
   		}
   		hldNum <- NULL
   		for(i in 1:length(pwc[1,])) hldNum <- c(hldNum, is.numeric(pwc[1,i]))
   		pwc[,hldNum] <- round(pwc[,hldNum], digits = digits)		
   		return(pwc)
   		}
   	nFits <- x$nFits
   	balTabList <- vector(mode = "list", length = nFits)
   	#if(x$estimand == "ATT")
   	if(collapse.to == "pair") cat(paste("Note that `tx' refers to the category specified as the treatATT, ", x$treatATT, ".\n\n", sep = ""))
   	for(i in 1:nFits) balTabList[[i]] <- bal.table(x$psList[[i]], digits = digits)
   	#if(x$estimand == "ATT") 
   	names(balTabList) <- x$levExceptTreatATT
   	for(i in 1:length(balTabList)){
   		for(j in 1:length(balTabList[[i]])){
   			balTabList[[i]][[j]] <- data.frame(var = row.names(balTabList[[i]][[j]]), balTabList[[i]][[j]], control = names(balTabList)[i], stop.method = names(balTabList[[i]])[j])
   		}
   	}
   	
   	nonTreatATT <- x$levExceptTreatATT
   	newBalTabList <- vector(mode = "list", length = (length(x$stopMethods) + 1) * length(nonTreatATT))
   	cnt <- 1
   	for(i in 1:length(nonTreatATT)){
   		newBalTabList[[cnt]] <- balTabList[[i]][["unw"]]
   		cnt <- cnt + 1
   	}
   	for(i in 1:length(x$stopMethods)){
   		for(j in 1:length(nonTreatATT)){
   			newBalTabList[[cnt]] <- balTabList[[nonTreatATT[j]]][[paste(x$stopMethods[i], ".ATT", sep = "")]]
   			cnt <- cnt + 1
   		}
   	}
   	balTabList <- do.call(rbind, newBalTabList)
   	row.names(balTabList) <- NULL
   	
   	balTabList$stop.method <- gsub(".ATT", "", balTabList$stop.method)

   	if(collapse.to == "pair") {
   		if(!is.null(subset.var)) balTabList <- subset(balTabList, var %in% subset.var | var %in% paste(subset.var, "<NA>", sep = ":"))
   		if(!is.null(subset.treat)) {
   			if(x$estimand == "ATE") balTabList <- subset(balTabList, tmt1 %in% subset.treat | tmt2 %in% subset.treat)
   			else balTabList <- subset(control %in% subset.treat)
   			}
   		if(!is.null(subset.stop.method)) balTabList <- subset(balTabList, stop.method %in% subset.stop.method)
   		balTabList <- subset(balTabList, abs(std.eff.sz) >= es.cutoff | is.na(std.eff.sz) | is.nan(std.eff.sz))
   		balTabList <- subset(balTabList, ks >= ks.cutoff | is.na(ks) | is.nan(ks))
   		balTabList <- subset(balTabList, p <= p.cutoff | is.na(p) | is.nan(p))
   		balTabList <- subset(balTabList, ks.pval <= ks.p.cutoff | is.na(ks.pval) | is.nan(ks.pval))
   		hldNum <- NULL
   		for(i in 1:length(balTabList[1,])) hldNum <- c(hldNum, is.numeric(balTabList[1,i]))
   		balTabList[,hldNum] <- round(balTabList[,hldNum], digits = digits)	   		
   		return(balTabList)
   		}
   	if(collapse.to == "covariate") {
   		colCov <- data.frame(var = rep(unique(balTabList$var), length(unique(balTabList$stop.method))), max.std.eff.sz = NA, min.p = NA, max.ks = NA, min.ks.pval = NA, stop.method = rep(c("unw",x$stopMethods), each = length(unique(balTabList$var))))
   		
   		for(i in 1:nrow(colCov)){
   			colCov$max.std.eff.sz[i] <- with(subset(balTabList, var == colCov$var[i] & stop.method == colCov$stop.method[i]), max(abs(std.eff.sz), ...))
   			colCov$min.p[i] <- with(subset(balTabList, var == colCov$var[i] & stop.method == colCov$stop.method[i]), min(p, ...))
   			colCov$max.ks[i] <- with(subset(balTabList, var == colCov$var[i] & stop.method == colCov$stop.method[i]), max(ks, ...))
   			colCov$min.ks.pval[i] <- with(subset(balTabList, var == colCov$var[i] & stop.method == colCov$stop.method[i]), min(ks.pval, ...))
   		}
   		colCov <- subset(colCov, abs(max.std.eff.sz) >= es.cutoff | is.na(max.std.eff.sz) | is.nan(max.std.eff.sz))
   		colCov <- subset(colCov, max.ks >= ks.cutoff | is.na(max.ks) | is.nan(max.ks))
   		colCov <- subset(colCov, min.p <= p.cutoff | is.na(min.p) | is.nan(min.p))
   		colCov <- subset(colCov, min.ks.pval <= ks.p.cutoff | is.na(min.ks.pval) | is.nan(min.ks.pval))   
   		hldNum <- NULL
   		for(i in 1:length(colCov[1,])) hldNum <- c(hldNum, is.numeric(colCov[1,i]))
   		colCov[,hldNum] <- round(colCov[,hldNum], digits = digits)	   		   		
   		return(colCov)
   	}
   	if(collapse.to == "stop.method"){
   		colStop <- data.frame(max.std.eff.sz = NA, min.p = NA, max.ks = NA, min.ks.pval = NA, stop.method = unique(balTabList$stop.method))
   		for(i in 1:nrow(colStop)){
   			colStop$max.std.eff.sz[i] <- with(subset(balTabList, stop.method == colStop$stop.method[i]), max(abs(std.eff.sz), ...))
   			colStop$min.p[i] <- with(subset(balTabList, stop.method == colStop$stop.method[i]), min(p, ...))
   			colStop$max.ks[i] <- with(subset(balTabList, stop.method == colStop$stop.method[i]), max(ks, ...))
   			colStop$min.ks.pval[i] <- with(subset(balTabList, stop.method == colStop$stop.method[i]), min(ks.pval, ...))
   		}
   		colStop <- subset(colStop, abs(max.std.eff.sz) >= es.cutoff | is.na(max.std.eff.sz) | is.nan(max.std.eff.sz))
   		colStop <- subset(colStop, max.ks >= ks.cutoff | is.na(max.ks) | is.nan(max.ks))
   		colStop <- subset(colStop, min.p <= p.cutoff | is.na(min.p) | is.nan(min.p))
   		colStop <- subset(colStop, min.ks.pval <= ks.p.cutoff | is.na(min.ks.pval) | is.nan(min.ks.pval))     	
   		hldNum <- NULL
   		for(i in 1:length(colStop[1,])) hldNum <- c(hldNum, is.numeric(colStop[1,i]))
   		colStop[,hldNum] <- round(colStop[,hldNum], digits = digits)	 	
   		return(colStop)
   	}
   }
   	
