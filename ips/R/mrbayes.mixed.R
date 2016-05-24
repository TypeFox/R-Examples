mrbayes.mixed <- function(x, file, nst = 6, 
                          rates = "invgamma", ngammacat = 4, 
                          nruns = 2, ngen = 1000000, 
                          printfreq = 100, samplefreq = 10,  
                          nchains = 4, savebrlens = "yes", 
                          temp = 0.2, burnin = 10, 
                          contype = "allcompat", 
                          run = TRUE){
	
	# get numbers of taxa and characters
	ntax <- dim(x)[1]
	ncha <- dim(x)[2]
	p2 <- min(which(x[1, ] %in% c("0", "1")))
	
	# datatype
	# --------
	datatype <- paste("datatype=mixed(DNA:1-", p2-1, ",Standard=", 		p2, "-", ncha, ")", sep = "")
	
	# assemble nexus file
	nexus <- vector(length = ntax + 11)
	nexus[1] <- "#nexus"
	nexus[c(2, 6, 8, ntax + 11)] <- ""
	nexus[3] <- "begin data;"
	nexus[4] <- paste("\tdimensions ntax=", ntax, " nchar=", 		ncha, 	";", sep = "")
	nexus[5] <- paste("\tformat", datatype, "missing=N gap=-;")
	nexus[7] <- "matrix"
	for (i in 1:ntax){
		s <- paste(x[i,], collapse = "")
		s <- paste(rownames(x)[i], toupper(s))
		nexus[i + 8] <- s
	}
	nexus[ntax + 9] <- ";"
	nexus[ntax + 10] <-"end;"
	
	
	# write MrBayes block
	bayes <- vector(length = 4)
	bayes[1] <- "begin mrbayes;"
	bayes[2] <- paste("\tlset nst=", nst, " rates=", 		rates, " ngammacat=", ngammacat,";", sep = "")
	bayes[3] <- paste("\tmcmc nruns=", nruns, " ngen=", 		as.integer(ngen), " printfreq=", printfreq, 		" samplefreq=", samplefreq, " nchains=", 		nchains, " savebrlens=", savebrlens, 			" temp=", temp, ";", sep = "")
	bayes[4] <- "end;"
	
	nexus <- c(nexus, bayes)
	write(nexus, file)
	
	# start mrbayes
	if ( run ){
		if (.Platform$OS.type == "unix"){
		  system(paste("./mb > execute", file))
		} else {
      system(paste("mrbayes ", file, ".bayes", sep = ""))
		}
		tr <- read.nexus(paste(file, ".con", sep = ""))
		return(tr)
	}
}
