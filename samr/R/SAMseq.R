######################################################################
#		A concise version of samr for sequencing data
######################################################################
SAMseq <- function(x, y, censoring.status = NULL, 
	resp.type = c("Quantitative", "Two class unpaired", "Survival", "Multiclass", "Two class paired"), 
	geneid = NULL, genenames = NULL, nperms = 100, random.seed = NULL, nresamp = 20, fdr.output = 0.20)
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
	data = list(x = x, y = y, censoring.status = censoring.status, geneid = geneid, genenames = genenames)
	
	## call samr with approprate parameters
	samr.obj = samr(data, resp.type = resp.type, assay.type = "seq", nperms = nperms, return.x = TRUE, 
		random.seed = random.seed, nresamp = nresamp, nresamp.perm = nresamp)
		
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
		} else
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
		if (resp.type == "Quantitative" | resp.type == "Survival")
		{
			rang = 4 : 7
		}
		
		siggenes.table$genes.up[, rang] <- round(as.numeric(siggenes.table$genes.up[, rang]), 3)
		siggenes.table$genes.lo[, rang] <- round(as.numeric(siggenes.table$genes.lo[, rang]), 3)
		
		tname <- colnames(siggenes.table$genes.up)
		siggenes.table$genes.up <- siggenes.table$genes.up[, 
			(tname != "Row") & (tname != "Numerator(r)") 
			& (tname != "Denominator(s+s0)")]
		tname <- colnames(siggenes.table$genes.lo)
		siggenes.table$genes.lo <- siggenes.table$genes.lo[, 
			(tname != "Row") & (tname != "Numerator(r)") 
			& (tname != "Denominator(s+s0)")]
	}

	## construct the return values
	out = list(samr.obj = samr.obj, del = del, delta.table = delta.table, siggenes.table = siggenes.table)
	out$call = this.call
	class(out) = "SAMoutput"
	
	return(out)
}

######################################################################
#		implementation of the generic function print
######################################################################
print.SAMoutput = function(x, ...)
{
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

######################################################################
#		implementation of the generic funciton plot
######################################################################
plot.SAMoutput = function(x, ...)
{
	samr.plot(x$samr.obj, x$del)
	invisible()
} 
