sort.data.frame<-function(x,decreasing,by,...){
	signs<-c()
	get.signs<-function(f){
		l<-length(signs)
		lf<-length(f)
		if(as.character(f[[1]]) %in% c("+","-")){
			signs[l+1]<<-as.character(f[[1]])
			names(signs)[l+1]<<-if(is.call(f[[lf]])) format(f[[lf]]) else as.character(f[[lf]])
			if(lf>2){
				if(!is.symbol(f[[2]]))
					get.signs(f[[2]])
				if(lf>=3 && !is.symbol(f[[3]]))
					get.signs(f[[3]])	
			}			
		}
	}
	bycall<-as.list(match.call()[-1])$by
	bycall<-gsub(" ","",if(is.call(bycall)) format(bycall) else as.character(bycall))
	if(substr(bycall,1,1) != "~")
		by<-as.formula(paste("~",bycall))
	
	if(length(gregexpr("~",format(by))[[1]])>1.5)
		stop("formula must be one-sided")	
	
	# If the first character is not + or -, add +
	formc <- as.character(by[2]) 
	if(!is.element(substring(formc, 1, 1), c("+", "-")))
		formc <- paste("+", formc, sep = "")
	
	f<-as.formula(paste("~1 ",formc))
	
	
	get.signs(f[[2]])
	signs<-rev(signs)
	mf<-as.data.frame(lapply(names(signs),function(nm) eval(parse(text=nm),x, parent.frame())))
	names(mf)<-names(signs)
	if(nrow(mf)!=nrow(x))
		stop("error: 'by' specification implies incorrect numbers of rows")
	
	# Build a list of arguments to pass to "order" function
	calllist <- list()
	for(i in 1:length(signs)){
		varsign <- signs[i]
		v<-mf[[names(signs)[i]]]
		if(is.factor(v)){
			if(varsign == "-") {
				calllist[[i]] <- -rank(v)
			} else {
				calllist[[i]] <- rank(v)
			}
		} else {
			if(varsign == "-") {
				calllist[[i]] <- -v
			} else {
				calllist[[i]] <- v
			}
		}
	}
	return(x[do.call("order", calllist), , drop=FALSE])
}
