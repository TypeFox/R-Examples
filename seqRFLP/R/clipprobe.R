clipprobe <-
function (fil, forProbe, bacProbe, tol = 3, clipped.only = TRUE) 
{
    if(length(fil) == 1){
	fil <- c(">inputseq",fil)
	fil <- as.fasta(fil)
	}
	if(length(fil) >= 2){
	    if(!inherits(fil, "fasta")){
	       stop("The input sequences must be in class \"fasta\".")
	    }
	}
	
	for (i in 1:(length(fil)/2)) {
        stri = toupper(fil[i * 2])
        re = ""
        forw = min(findprobe(stri, forProbe, tol))
        back = max(findprobe(stri, bacProbe, tol))
        
		if (back == 0 & forw == 0) {
		
            revstri = revComp(stri)
            forw = min(findprobe(revstri, forProbe, tol))
            back = max(findprobe(revstri, bacProbe, tol))
            
			if (forw > 0 | back > 0){ 
                stri = revstri
			   }
        }
        if (back > 0 & forw > 0) {
            stri = substr(stri, forw, (back + nchar(bacProbe) - 1))
            re = "_cliped"
        }
        nam = paste(fil[(2 * i) - 1], re, sep = "")
		
        if (i == 1){ 
            DNA = c(nam, stri)
			}
			
        if (i > 1) {
            DNA = c(DNA, nam, stri)
			}
    }
    if (clipped.only){ 
        DNA = DNA[sort(c(grep("_cliped", DNA), grep("cliped", DNA) + 1))]
	   }
	DNA = gsub("_cliped$", "", DNA)
	class(DNA) <- "fasta"
	return(DNA)
}

