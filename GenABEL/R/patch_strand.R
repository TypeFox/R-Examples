
"patch_strand" <- function(data,snpid,strand,based_on="snpnames", quiet = TRUE)
{
	
	if (class(data) != "gwaa.data" && class(data) != "snp.data") 
		stop("data must be of gwaa.data or snp.data-class")
	possible_bases <- c("map","snpnames")
	if (!any(based_on == possible_bases)) 
		stop(paste("based_on must be one of: ",possible_bases))
	if (length(snpid) != length(strand)) 
		stop("length of 'snpid' is npot equal to length of 'strand'")
	if (!is.character(strand) && !is.factor(strand)) 
		stop("strand must be character or factor")
	strand <- as.character(strand)
	if (length(levels(as.factor(strand)))>3)
		stop("no more than three levels of strand (+,-,u) are allowed")
	if (!all(levels(as.factor(strand)) %in% c("+","-","u")))
		stop (paste("only three levels of strand (+,-,u) are allowed; now",levels(as.factor(strand))))
	if (length(snpid) != unique(length(snpid))) 
		stop("non-unique snpid's")
	
	if (class(data) == "gwaa.data")
		wdata <- data@gtdata
	else
		wdata <- data
	
	if (based_on == "map")
		ga_snpid <- as.integer(wdata@map)
	else if (based_on == "snpnames")
		ga_snpid <- wdata@snpnames
	else stop("wrong based_on argument")
	
	if (class(ga_snpid) != class(snpid))
		warning("classes of snpid and internal data snpid do not match")
	
	ga_matched_snps <- which(ga_snpid %in% snpid)
	matched_snps <- which(snpid %in% ga_snpid)
	
	if (!quiet) cat("identified",length(ga_matched_snps),"SNPs to be patched\n")
	new_strand <- strand[matched_snps]
	old_strand <- as.character(wdata@strand)
	if (!quiet) cat("Changes table:\n")
	if (!quiet) print(table(old_strand[ga_matched_snps],new_strand))
	if (!quiet) cat("changing strand for",sum(old_strand[ga_matched_snps] != new_strand),"SNPs\n")
	old_strand[ga_matched_snps] <- new_strand
	raw_strand <- rep(0,length(old_strand))
	raw_strand[old_strand == "+"] <- 1
	raw_strand[old_strand == "-"] <- 2
	
	rawVal <- as.raw(raw_strand)
	names(rawVal) <- snpnames(wdata)
	wdata@strand <- new("snp.strand",rawVal)
#	print(table(as.character(wdata@strand[ga_matched_snps]),new_strand))
	
	if (class(data) == "gwaa.data")
		wdata <- new("gwaa.data",phdata=data@phdata,gtdata=wdata)

	if (!quiet) cat("... done\n")
	return(wdata)
}