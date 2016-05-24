
callList2terms=function(cl)
{	lst=cl
	if(!is.list(lst)) return(deparse(lst))
	
	first=lst[[1L]]
	if(length(lst)==1L) return(as.character(lst[[1L]]))
	
	if(identical(first, colon) && length(lst)==3L ){
		return(c(Recall(as.list(lst[[2L]])), Recall(as.list(lst[[3L]]))))
	}
	tmp=sapply(lst[-1], as.character)
	tmpn = names(lst[-1])
	if(!is.null(tmpn)) for(i in which(nchar(tmpn)>0L)) tmp[i] = paste(tmpn[i], '=', tmp[i])
	paste0(as.character(first), "(", paste(tmp, collapse=', '), ')')
}

splitTerm=function(term)
{
	tcall=as.list(as.call(parse(text=term))[[1]])
	unique(callList2terms(tcall))
}

sortTerm=function(term, priority)
{
	singleTerms = splitTerm(term)
	if(missing(priority)) singleTerms=sort(singleTerms) else{
		idx=singleTerms%in%priority
		singleTerms=c(sort(singleTerms[idx]), sort(singleTerms[!idx]))
	}
	paste0(singleTerms, collapse=':')
}
