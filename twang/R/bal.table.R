## Generic function for extracting balance tables from ps and other objects
bal.table <- function(x, digits = 3, collapse.to = c("pair","covariate","stop.method")[1], subset.var = NULL, subset.treat = NULL, subset.stop.method = NULL, es.cutoff = 0, ks.cutoff = 0, p.cutoff = 1, ks.p.cutoff = 1, timePeriods = NULL, ...){
	if(!(class(x) %in% c("mnps", "iptw", "mniptw"))){
   bal.tab <- bal.table.ps(x, digits = digits)
   return(bal.tab)
   }
   else if(class(x) == "iptw"){
   	if(is.null(timePeriods)) timePeriods <- 1:length(x$psList)
   	for(i in timePeriods){
   		cat("Balance at time ", x$uniqueTimes[i], ":\n")
   		print(bal.table.ps(x$psList[[i]], digits = digits))
   		cat("\n")
   	}
   }
   else if(class(x) == "mnps"){
   	bal.table.mnps(x=x, digits = digits, collapse.to = collapse.to, subset.var = subset.var, subset.treat = subset.treat, subset.stop.method = subset.stop.method, es.cutoff = es.cutoff, p.cutoff = p.cutoff, ks.p.cutoff = ks.p.cutoff, ...)
   }
   else if(class(x) == "mniptw"){
   	if(is.null(timePeriods)) timePeriods <- 1:length(x$psList)
   	for(i in timePeriods){
   		cat("Balance at time ", x$uniqueTimes[i], ":\n")
   		print(bal.table.mnps(x$psList[[i]], digits = digits, collapse.to = collapse.to, subset.var = subset.var, subset.treat = subset.treat, subset.stop.method = subset.stop.method, es.cutoff = es.cutoff, p.cutoff = p.cutoff, ks.p.cutoff = ks.p.cutoff, ...))
   		cat("\n")
   		
   	}
   }
   	
}


