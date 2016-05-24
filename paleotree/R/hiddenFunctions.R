#hidden functions, all having to do with model fitting junk

make_paleotreeFunc<-function(f,parnames,parbounds){
	#this should be a hidden function!!
	np<-length(parnames)	#number of parameters 
	#just automates making a paleotreeFunc
	attr(f,"class")<-c("paleotreeFunc",class(f))
	#parnames
	value<-as.character(parnames)
	if(any(is.na(value))){stop("NA values in parnames replacement")}
	if(any(duplicated(value))){stop("Duplicated names in parnames replacement")}
	attr(f,"parnames")<-parnames
	#parbounds
	if(!is.list(parbounds) | !length(parbounds)==2){stop("parbounds needs to be a list composed of two vectors")}
	lower<-as.numeric(parbounds[[1]])
	if(length(lower)!=np){stop("length of new lower parbounds not equal to number of parameters")}
	if(any(is.na(lower))){stop("NA values in lower parbounds replacement")}
	upper<-as.numeric(parbounds[[2]])
	if(length(upper)!=np){stop("length of new upper parbounds not equal to number of parameters")}
	if(any(is.na(upper))){stop("NA values in upper parbounds replacement")}
	names(parbounds[[1]])<-names(parbounds[[2]])<-parnames
	attr(f,"parbounds")<-parbounds
	#
	attr(f,"np")<-length(parnames)
	return(f)
	}


constrainParsePaleo<-function(formula, names.lhs, names.rhs,extra=NULL){
	#this should be a hidden function!!
	#based on Rich FitzJohn's functions for diversitree 10-22-13
		#all comments with double ## are Rich's original comments, with some editting
	#
	## Parsing constraints:
	## The LHS of a formula must be a single variable name that exists in "names.lhs"
	##
	## The RHS can be one of numeric value or an expression
	## If it is an expression, then all variable names must be found in
	## names.rhs (or perhaps in the containing environment - check in the future?)
	# 
	formula <- as.formula(formula)
	if ( length(formula) != 3L ) {stop("Invalid formula")}
	lhs <- formula[[2]]	#of type symbol
	rhs <- formula[[3]]	#of type language
	 ## Checking the lhs is easy: is the lhs in the list of allowable
	 ## names and of length 1? Anything that does not match this is invalid.
	if ( !is.name(lhs) ) {stop("Invalid target on LHS of formula" )}	#If one result term in formula, is.name = true
	lhs.is.target <- is.na(match(as.character(lhs), names.lhs))
	 ## Checking the rhs is more difficult. We are OK if any of the
	 ## following is met:
	 ## Numeric values (checked at the end)
	 ## If all vars in rhs are in names.rhs
	 ## There is a single var and it is in extra
	 ## Failing that, if the rhs is a single variable that does exist in
	 ## the calling environment.
	if( is.language(rhs) ) {
		vars <- all.vars(rhs)
		ok <- (all(vars %in% names.rhs) || length(vars) == 1 && vars %in% extra)
		if( !ok && length(vars) == 1 ) {
			e <- parent.frame()
			if( exists(vars, e) ) {
				rhs <- get(vars, e)
				ok <- TRUE
				}
   			}
		if( !ok ){stop("Invalid RHS of formula:\n\t", as.character(rhs))}
		if( as.character(lhs) %in% vars ){stop("LHS cannot appear in RHS")}
	}else{
		if( !is.numeric(rhs) ) {stop("RHS must be expression, variable or number")}
		}
	res <- list(lhs, rhs)
  	attr(res, "lhs.is.target") <- lhs.is.target
	return(res)
	}

expandConstrainForm<-function(formula,breakNames,nparcat){
	#another function to be hidden at all costs
	#take an expression like p.all~q.all 
	formula <- as.formula(formula)
	if ( length(formula) != 3L ) {stop("Invalid formula")}
	lhs <- formula[[2]]	#of type symbol
	rhs <- formula[[3]]	#of type language
	if ( !is.name(lhs) ) {stop("Invalid target on LHS of formula" )}	#If one result term in formula, is.name = true
	#for now, we can only have one rhs term
	if ( !is.name(rhs) ) {stop("Can't have more than single, non-modified RHS term with match/all" )}
	#break 'em
	lhs<-as.vector(sapply(all.vars(lhs),function(x) unlist(strsplit(x,".",fixed=TRUE))))
	rhs<-as.vector(sapply(all.vars(rhs),function(x) unlist(strsplit(x,".",fixed=TRUE))))
	#right length?
	if(!length(lhs)==nparcat){stop(paste("formula LHS doesn't have right number of categories in label in",format(formula)))}
	if(!length(rhs)==nparcat){stop(paste("formula RHS doesn't have right number of categories in label in",format(formula)))}
	#find alls and matches
	lhsAll<-lhs=="all"; rhsAll<-rhs=="all"; lhsMatch<-lhs=="match"; rhsMatch<-rhs=="match"
	#expand all 'all' statements
	lhs<-expand.grid(lapply(1:nparcat,function(x) if(lhsAll[x]){unique(breakNames[,x])}else{lhs[x]}))
	rhs<-expand.grid(lapply(1:nparcat,function(x) if(rhsAll[x]){unique(breakNames[,x])}else{rhs[x]}))
	#and like that, we no longer need both lhs and rhs
	lhs<-list(rbind(rhs,lhs))
	#expand all match statements
	if(any(lhsMatch) | any(rhsMatch)){
		#first, test that the matches are symmetric
		if(!which(lhsMatch)==which(rhsMatch)){stop(paste("'match' statements aren't symmetric in",format(formula)))}
		#we only need to worry about lhs now, but need to account for each lhs matrix made
		replaceList<-list()
		for(j in which(lhsMatch)){	
			for(k in unique(breakNames[,j])){
				for(i in 1:length(lhs)){
					replaceMatch<-lhs[[i]]
					replaceMatch[,j]<-k
					replaceList<-c(replaceList,list(replaceMatch))
					}
				}
			}
		lhs<-replaceList
		}
	#EXPANSION DONE
	#get rid of any double entries within each set of equivalancies
	lhs<-lapply(lhs,unique)
	#OKAY now to turn these into formulas, first make strings
	#turning these freaking data.frames back into strings is more of a nightmare than expected
	rhs<-lapply(lhs,function(x) x[-1,])
	rhs<-lapply(rhs,function(x) apply(x,1,function(y) paste0(c(sapply(y,as.character)),collapse=".")))
	lhs<-lapply(lhs,function(x) x[1,])
	lhs<-lapply(lhs,function(x) paste0(c(sapply(x,as.character)),collapse="."))
	#now make a formula for each unit of rhs
	res<-as.vector(sapply(1:length(lhs),function(x) sapply(lhs[[x]],function(y) paste(rhs[[x]],y,sep=" ~ "))))
	return(res)
	}