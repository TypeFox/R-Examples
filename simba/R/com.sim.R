"com.sim" <-
function(veg, subs, simil="soerensen", binary=TRUE, permutations=1000, alpha=0.05, bonfc=TRUE, ...) {   
    strata <- unique(subs)
	stratan <- length(strata)
	comb <- outer(strata, strata, "paste")
	comb <- comb[row(comb)<col(comb)]
	comb <- t(sapply(strsplit(comb, split=" "), function(x) x[1:2]))
	if(binary){
	   tmp <- lapply(c(1:nrow(comb)), function(x) diffmean(as.numeric(sim(veg[subs==(comb[x,1]),], method=simil)), as.numeric(sim(veg[subs==(comb[x,2]),], method=simil,))))
	}
	else{
	   tmp <- lapply(c(1:nrow(comb)), function(x) diffmean(as.numeric(vegdist(veg[subs==(comb[x,1]),], method=simil)), as.numeric(vegdist(veg[subs==(comb[x,2]),], method=simil,))))
	}
	comb <- data.frame(comb)
	names(comb) <- c("X", "Y")
	##bonferroni correction
	if (bonfc) {
	alpha <- alpha/stratan
	}
	comb$mean.x <- as.numeric(sapply(sapply(tmp, function(x) x[4]), function(x) x[1]))
	comb$mean.y <- as.numeric(sapply(sapply(tmp, function(x) x[5]), function(x) x[1]))
	comb$diff <- as.numeric(sapply(sapply(tmp, function(x) x[2]), function(x) x[1]))
	comb$sig <- as.numeric(sapply(sapply(tmp, function(x) x[7]), function(x) x[1]))
	comb$sigs <- ifelse(comb$sig >= alpha, "ns", "*")
    comb$F <- as.numeric(sapply(sapply(tmp, function(x) x[6]), function(x) x[1]))
	comb$sigF <- as.numeric(sapply(sapply(tmp, function(x) x[8]), function(x) x[1]))
    comb$sigsF <- ifelse(comb$sigF >= alpha, "ns", "*")
    con.tab <- mama(comb[,c(1:2,7)])
	res <- list(call=match.call(), method=simil, out=comb, cons=con.tab, strata=stratan, permutations=permutations)
	class(res) <- "cslist"
	res
    }