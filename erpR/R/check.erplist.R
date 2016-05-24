check.erplist<-function(erplist = NULL){
	# preliminary checks
	if (is.null(erplist)){
	stop("an erplist object containing ERP data frames must be specified!", call.=F)
	}
	
	#### SET A SERIES OF OBJECTS COLLECTING THE CHECKS MADE
	errors=list(
	dupl.error = FALSE
	, no.number.error = FALSE
	, incons.number.error = FALSE
	
	)
	
	## CHECK DUPLICATED
	dupl=duplicated(names(erplist))
	if (any(dupl)){
		errors$dupl.error=TRUE
		dupl.names=names(erplist)[dupl]
		dupl.names=paste(dupl.names, sep=" ", collapse="\n")
		cat("WARNING!\nThe erplist contains the following duplicated objects:\n", dupl.names, "\n\n")
	}
	
	## CHECK NO NUMBER IN NAMES
	erplist.len=1:length(erplist)
	nonumber=erplist.len[!erplist.len%in%grep("[0-9]", names(erplist))]
	
	if (length(nonumber>0)){
		errors$nonumber.error=TRUE
		nonumber.names=names(erplist)[nonumber]
		nonumber.names=paste(nonumber.names, sep=" ", collapse="\n")
		cat("WARNING!\nThe erplist contains the following objects with incorrect naming:\n", nonumber.names, "\n\n")
		cat("An erp object name should be composed by a \"base\" and a \"number\" ( see help(erpR) )\n\n")
	}
	
	# CHECK UNBALANCED SUBJECTS
	numbers= as.numeric(gsub("[^[:digit:]]","", names(erplist))) #notice the double square parantheses. It is the key 
	numbers.dat=data.frame(table(numbers))
	numbers.dat=apply(numbers.dat, 2, as.character)
	numbers.dat=apply(numbers.dat, 2, function(x){sprintf("%3s", x)})	
	
	if (length(unique(numbers.dat[,"Freq"]))!=1){
		errors$incons.number.error = TRUE
		# the line abovecheck if there is consistency. If there is, then all numbers have the same Freq (see numbers.dat columns)
		# and then the length of unique is equal to 1. Other wise it is not equal to 1.
		cat("WARNING!\nThe erplist contains some inconsistencies in the objects:\n")
		cat("num", "Freq", "\n")
		for (i in 1:dim(numbers.dat)[1]){
			cat(numbers.dat[i,], "\n")
		}
		cat("\n")
	}

	# RETURN THAT EVERYTHING IS OK
	if(all(!unlist(errors))){
		cat("No problem found in the erplist\n\n")
	}
		
}

	