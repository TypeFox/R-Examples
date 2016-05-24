### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### textmatrix
### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### dependencies: library("RStem")

### HISTORY
### 
### 2009-09-14
###    * reduced data load to compare only input vector with rest of
###      the matrix -- thanks to input from Yue Shan, National Cheng Kung University, Taiwan 
### 2009-08-19
###    * created from input from Robert Koblischke / Reinhard Dietl and the
###      folks from the Idiom project
###

associate <- function ( textmatrix, term, measure="cosine", threshold=0.7 ) {

	# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
	# calc similarity table of all terms vectors to each other
	
	if (measure == "cosine" ) { 
		term2term=cosine(textmatrix[term,],t(textmatrix))
	} else if (measure == "pearson") {
		term2term=cor(textmatrix[term,], t(textmatrix), method="pearson")
	} else if (measure == "spearman") {
		term2term=cor(textmatrix[term,], t(textmatrix), method="spearman")
	} else {
		stop("[associate] ERROR: This closeness measure is not supported.")
	}
	
	# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
	# get index of the terms that correlate with computer higher than .7
	# options: (all) or (>threshold)
	
	# only the ones above threshold
	ix = which(term2term>threshold)
	ranks = sort(term2term[ix], decreasing=TRUE)

	# remove 'term' from list
	ranks = ranks[-which(names(ranks)==term)]
	
	if(length(ranks)>0) return (ranks) else return(NULL)

} # associate ()

