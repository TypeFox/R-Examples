"metascore" <-
function(...) {
	names <- as.character(substitute(list(...)))[-1]
# generate global table
	first <- 1
	for (i in names) {
		print(i)
		if (first) {
			gtab <- get(i)
			cnams <- names(gtab)
			if (checkcls(gtab)) stop(paste("error evaluating object",i))
			gtab <- data.frame(name=gtab$name,beta=gtab$beta,seb=gtab$seb,P=gtab$P,Pgc=gtab$Pgc,coding=gtab$coding,strand=gtab$strand,study=rep(i,length(gtab$name)),stringsAsFactors=F)
			first <- 0
		} else {
			tmp <- get(i)
			if (any(names(tmp) != cnams)) 
				stop(paste("column names in obejct",i,"(",names(tmp),") do not match to (",cnams,")"))
			if (checkcls(tmp)) stop(paste("error evaluating object",i))
			tmp <- data.frame(name=tmp$name,beta=tmp$beta,seb=tmp$seb,P=tmp$P,Pgc=tmp$Pgc,coding=tmp$coding,strand=tmp$strand,study=rep(i,length(tmp$name)),stringsAsFactors=F)
			gtab <- rbind(gtab,tmp)
			
		}
		
	}
# analyse global table
	gtab$passed <- 0
	usnps <- unique(gtab$name)
	for (i in usnps) {
		pres <- gtab[gtab$name==i,]
		cod <- pres$coding[1]
		str <- pres$strand[1]
		nstu <- dim(pres)[1]
		if (nstu>1) {
			for (j in c(2:nstu)) {
				print(pres)
				if (pres$strand[j] != str) {
		  			revs <- alleleID.raw2char()[as.character(alleleID.revstrand())]
					names(revs) <- alleleID.raw2char()[names(alleleID.revstrand())]
					print(revs)
					pres$coding[j] <- revs[pres$coding[j]]
				}
				print(pres)
			}
		}
	}
}

"checkcls" <-
function(data) {
	err <- 0
	nams <- names(data)
	reqnams <- c("Pgc","name","P","beta","seb","coding","strand")
	m <- match(reqnams,nams)
	if (any(is.na(m))) {
		cat("arguments do not match:\n")
		cat("required arguments are:\n")
		print(reqnams)
		cat("observed arguments are:\n")
		print(nams)
		err <- 1
	} 
	err
}
