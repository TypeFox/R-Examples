findprobe <-
function (dna, probe, tol = 3) 
{
    xxx <- 
	function(x) 
	{sum((aa[x + seq_along(bb)])!= bb)}
	
	aa <- strsplit(dna, "")[[1]]
    bb <- strsplit(probe, "")[[1]]
    cc <- sapply(0:(length(aa) - length(bb)), xxx)
	
    pb = which(cc <= tol)
    
	if (table(cc <= tol)["FALSE"] == length(cc)) 
	{ pb = 0 }
	return(pb)
}

