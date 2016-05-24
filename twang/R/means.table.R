means.table <- function(mnps, stop.method = 1, includeSD = FALSE, digits = NULL){

	if(is.numeric(stop.method)){
		if(!(stop.method %in% 1:length(mnps$stopMethods))){
			stop("If numeric, 'stop.method' must be an integer between 1 and the number of stop.method's used to fit the mnps object.")
		}
		
	}
	
	if(mnps$estimand == "ATE") meansTab <- with(mnps$psList[[1]]$desc$unw, data.frame(popMean = bal.tab$results$ct.mn))
	else {
		meansTab <- with(mnps$psList[[1]]$desc$unw, data.frame(theTreatedMean = bal.tab$results$tx.mn))
		if(includeSD) meansTab <- data.frame(meansTab, with(mnps$psList[[1]]$desc$unw, bal.tab$results$tx.sd))
		}


	if(mnps$estimand == "ATE") namesVec <- "pop.mean"
	else {
		namesVec <- paste(mnps$treatATT, "mean", sep = ".")
		if(includeSD) namesVec <- c(namesVec, paste(mnps$treatATT, "SD", sep = "."))
		}
	
	if(!is.numeric(stop.method)) {
		stop.method.long <- paste(stop.method, mnps$psList[[1]]$estimand, sep = ".")
		j <- match(stop.method.long, names(mnps$psList[[1]]$w))
		j <- j+1
		}
	else j <- stop.method + 1

	if(mnps$estimand == "ATE"){	
		for(i in 1:length(mnps$psList)){
			meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc$unw$bal.tab$results[,"tx.mn"])
			namesVec <- c(namesVec, paste("unwt", names(mnps$psList)[i], "mean", sep = "."))
			meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc[[j]]$bal.tab$results[,"tx.mn"])
			namesVec <- c(namesVec, paste("wt", names(mnps$psList)[i],"mean", sep = "."))
			meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc$unw$bal.tab$results[,"std.eff.sz"])
			namesVec <- c(namesVec, paste("unwt", names(mnps$psList)[i], "smd",sep = "."))
			meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc[[j]]$bal.tab$results[,"std.eff.sz"])
			namesVec <- c(namesVec, paste("wt", names(mnps$psList)[i], "smd", sep = "."))		
		
			if(includeSD){
				meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc$unw$bal.tab$results[, "tx.sd"])
				meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc[[j]]$bal.tab$results[,"tx.sd"])		
				namesVec <- c(namesVec, paste("unwt", names(mnps$psList)[i], "SD",  sep = "."))
				namesVec <- c(namesVec, paste("wght", names(mnps$psList)[i], "SD", sep = "."))
		
			}
			
		}	
		
	}
	
	if(mnps$estimand == "ATT"){
		for(i in 1:length(mnps$psList)){
			meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc$unw$bal.tab$results[,"ct.mn"])
			namesVec <- c(namesVec, paste("unwt", names(mnps$psList)[i], "mean", sep = "."))
			meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc[[j]]$bal.tab$results[,"ct.mn"])
			namesVec <- c(namesVec, paste("wt", names(mnps$psList)[i], "mean", sep = "."))			
			meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc$unw$bal.tab$results[,"std.eff.sz"])
			namesVec <- c(namesVec, paste("unwt", names(mnps$psList)[i], "smd", sep = "."))
			meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc[[j]]$bal.tab$results[,"std.eff.sz"])
			namesVec <- c(namesVec, paste("wt", names(mnps$psList)[i], "smd", sep = "."))
			
			if(includeSD){
				meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc$unw$bal.tab$results[, "ct.sd"])
				meansTab <- data.frame(meansTab, mnps$psList[[i]]$desc[[j]]$bal.tab$results[, "ct.sd"])
				namesVec <- c(namesVec, paste("unwt", names(mnps$psList)[i], "SD",  sep = "."))
				namesVec <- c(namesVec, paste("wght", names(mnps$psList)[i], "SD", sep = "."))
			}					
		}
		
	}
	
	names(meansTab) <- namesVec
	row.names(meansTab) <- row.names(mnps$psList[[1]]$desc$unw$bal.tab$results)
	
	if(is.null(digits))
	return(meansTab)
	else return(round(meansTab, digits = digits))
		
	
}