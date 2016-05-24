### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### sample.textmatrix.R
### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### created 2006-07-31

sample.textmatrix <- function (textmatrix, samplesize, index.return=FALSE) {

	if (class(textmatrix) != "textmatrix") {
		stop("[sample.textmatrix] ERROR: first argument not a textmatrix."); 
	}

	filelist = colnames(textmatrix)
	rnd_sample = sample(1:length(filelist), samplesize)
	if (index.return) {
		return(list(x=filelist[rnd_sample],ix=rnd_sample))
	} else {
		return(filelist[rnd_sample])
	}
	
} # sample.textmatrix()
