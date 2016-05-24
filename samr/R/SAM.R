SAM = function(x, y = NULL, censoring.status = NULL, 
	resp.type = c("Quantitative", "Two class unpaired", "Survival", 
		"Multiclass", "One class", "Two class paired", "Two class unpaired timecourse", 
		"One class timecourse", "Two class paired timecourse", 
		"Pattern discovery"),  
	geneid = NULL, genenames = NULL, s0 = NULL, s0.perc = NULL, 
	nperms = 100, center.arrays = FALSE, testStatistic = c("standard", 
		"wilcoxon"), time.summary.type = c("slope", "signed.area"), 
	regression.method = c("standard", "ranks"), return.x = TRUE, 
	knn.neighbors = 10, random.seed = NULL,
	logged2 = FALSE, fdr.output = 0.20, eigengene.number = 1)
{
	this.call <- match.call()
	xl.mode = "regular"
	xl.time = NULL
	xl.prevfit = NULL
	
       if(fdr.output <0 | fdr.output>1){
        stop("Error: fdr.output must be between 0 and 1")
       }

	## set gene id and gene names
	if (is.null(geneid))
	{
		geneid = as.character(1:nrow(x))
	}
	if (is.null(genenames))
	{
		genenames = paste("g", as.character(1:nrow(x)), sep = "")
	}
	
	## construct the data list
	data = list(x = x, y = y, censoring.status = censoring.status, 
		geneid = geneid, genenames = genenames, logged2 = logged2, 
		eigengene.number = eigengene.number)
	
	## call samr with approprate parameters
	samr.obj = samr(data, resp.type = resp.type, assay.type = "array",
		s0 = s0, s0.perc = s0.perc, nperms = nperms, center.arrays = center.arrays, 
		testStatistic = testStatistic, time.summary.type = time.summary.type, 
		regression.method = regression.method, return.x = return.x, 
		knn.neighbors = knn.neighbors, random.seed = random.seed)
		
	## construct a full delta table with no fold change constraints
	delta.table <- samr.compute.delta.table(samr.obj)
	
	## get the significant gene table
	siggenes.table <- del <- NULL
	delta.table <- delta.table[delta.table[, "# called"] > 0, , drop=FALSE]
	if (nrow(delta.table) > 0)
	{
		## find the appropriate delta
		oo <- which(delta.table[, "median FDR"] >= fdr.output)
		if (length(oo) > 0)
		{
			oo <- oo[length(oo)]
		}
		else
		{
			oo <- 1
		}
		
		## grab the right part of the the delta.table and get the significant gene table
		delta.table <- delta.table[oo : nrow(delta.table), , drop=FALSE]
		del <- delta.table[1, "delta"]
		siggenes.table <- samr.compute.siggenes.table(samr.obj, del, data, delta.table)
		
		## get the right range of table
		rang = 4 : 8
		if (resp.type == "Multiclass")
		{
			nclass = length(table(y))
			rang = 3 : (ncol(siggenes.table$genes.up))
		}
		if (resp.type == "Quantitative" | resp.type == "Pattern discovery" | resp.type == "Survival")
		{
			rang = 4 : 7
		}
		
		siggenes.table$genes.up[, rang] = round(as.numeric(siggenes.table$genes.up[, rang]), 3)
		siggenes.table$genes.lo[, rang] = round(as.numeric(siggenes.table$genes.lo[, rang]), 3)
		siggenes.table$genes.up = siggenes.table$genes.up[, -1]
		siggenes.table$genes.lo = siggenes.table$genes.lo[, -1]
	}

	## construct the return values
	out = list(samr.obj = samr.obj, del = del, delta.table = delta.table, siggenes.table = siggenes.table)
	out$call = this.call
	class(out) = "SAMoutput"
	
	return(out)
}

print.SAMoutput = function(x, ...) {
	cat("Call:\n")
	dput(x$call)
	mat1 = x$siggenes.table$genes.up
	cat("", fill = T)
	cat("Genes up", fill = T)
	print(mat1, quote = FALSE)
	mat2 = x$siggenes.table$genes.lo
	cat("", fill = T)
	cat("Genes down", fill = T)
	print(mat2, quote = FALSE)
	invisible()
}

plot.SAMoutput = function(x, ...) {
	samr.plot(x$samr.obj, x$del)
	invisible()
} 
