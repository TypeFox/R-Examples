frag.dat <-
function (fil, enznames, enzdata) 
{
	if(length(fil) == 1)
	{stop("Only the object in class \"fasta\" could be used.")}
	if(length(fil) >= 2){
	   if(!inherits(fil, "fasta")){
	      fil <- as.fasta(fil)
	   }
	}
	
	for (i in 1:(length(fil)/2)) {
        res.cut = enzCut(enznames = enznames, DNAsq = fil[2 * i], enzdata = enzdata)	
        cutSite = paste(res.cut$RFLP.site, collapse = ",")
        fragLength = paste(res.cut$RFLP.frag, collapse = ",")
        T5 = res.cut$TRFLP["T5"]
        T3 = res.cut$TRFLP["T3"]
        dft = data.frame(enznames = res.cut$enz$nam, recogSite = res.cut$enz$site, 
            cutting_Site = cutSite, fragment_Length = fragLength, 
            T5 = T5, T3 = T3)
        rownames(dft) = fil[2 * i - 1]
        if (i == 1) 
            dat = dft
        if (i > 1) 
            dat = rbind(dat, dft)
    }
    return(dat)
}

