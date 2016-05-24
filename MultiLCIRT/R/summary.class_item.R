summary.class_item <-function(object, ...){
	
# preliminaries
	r = length(out$lk)+1
	out = object
# print output	
	cat("\nCall:\n")
    print(out$call)
    cat("\nClustering output:\n")
	lg = rep(0,r-1); for(j in 1:(r-1)) lg[j] = length(out$groups[[j]])
	np = out$np0-out$np; pv = pchisq(out$height,np,lower.tail=FALSE)
	table = cbind(out$merge,out$height,np,pv,matrix(0,r-1,max(lg)))
	for(i in 1:(r-1)) table[i,6:(5+lg[i])] = sort(out$groups[[i]])
	names = c("items","","deviance","df","p-value","newgroup")
	for(i in 7:ncol(table)) names = c(names,"")
	colnames(table) = names
	names = NULL; for(i in 1:(r-1)) names = c(names,paste("Step",i))
	rownames(table) = names
	table[,c(3,5)] = round(table[,c(3,5)],3)
    print(table)
    cat("\n")

}