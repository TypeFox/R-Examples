#### Function edgesub as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

edgesub <-
function(x, pattern = "-", replacement = "?"){
	loc <- fmatch(dna = x, pattern = pattern)
	newdna <- paste(paste(rep(x = replacement, times = loc), collapse = ""), 
			  substring(x, loc+1, nchar(x)), sep = "")		  
	revdna <- reverse(newdna)		  
	locrev <- fmatch(revdna, pattern = pattern)		  
	revdnaresult <- paste(paste(rep(x = replacement, times = locrev), 
	       collapse = ""), substring(revdna, locrev+1, nchar(newdna)), sep = "")
	dnaresult <- reverse(revdnaresult)
	return(dnaresult)
}

