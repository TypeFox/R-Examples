# Author: ianfellows
###############################################################################

#d: the corpus for which term frequencies are calculated
#
#topN: If specified, only the 'topN' most frequent terms are returned. If more 
#      terms are requested than available, all terms are returned.
#
#percent: If specified, only the top 'percent' % of terms are returned. If more
#	  terms are requested than available, all terms are returned.
#
#sorted: A string specifying how to sort the terms.  'none' for no sorting,
#	'alpha' for alphanumeric sorting, and 'freq' for sorting by frequency.
#
#decreasing: If TRUE, terms are sorted in decreasing order, if FALSE, sorted
#	     ascending order.
#
#useDocFreq: If TRUE, the returned frequencies are for the total number of
#	documents in which the term occurs.  If false, they are the total 
#	number of occurrences.
#
#minFreq: Terms with *TOTAL* frequencies below this threshold will not be 
#	  included in the output.
#

term.freq <- function(
		d, 
		topN=0, #only include the first 'topN' terms
		percent=0,
		sorted=c("none", "alpha", "freq"),
		decreasing=FALSE,
		useDocFreq=FALSE,		
		minFreq=1) {
	
	dtm <- DocumentTermMatrix(d, control = list(tolower=FALSE, minWordLength=1));
	
	#remove terms whose total ocurrences < minFreq
	dtm <- dtm[ ,findFreqTerms(dtm, minFreq)];
	
  #clamp freqs of dtm between 0 and 1 for doc. frequency
	if (useDocFreq == TRUE) 
	{
		dtm <- dtm > 0
	}
	x <- apply(dtm,2,sum);

  sorted <- match.arg(sorted);

  #negate x instead of using decreasing = true, so that ties are sorted alphabetically increasing without reversing 'names'
  if (topN > 0) { #use absolute number of terms
	  o <- order(-x, names(x))[1:min(topN, length(x))];
	  x <- x[o];
  } else if(percent > 0) { #use percentage of terms
    o <- order(-x, names(x))[1:(length(x) * percent / 100)];
    x <- x[o];
  } 

  if(sorted == "alpha") {
    x <- x[order(names(x), decreasing=decreasing)];
  }
  else if (sorted == "freq") {
	  #x <- sort(x, decreasing=decreasing);
		
		#want to sort ties alphabetically
		#using the decreasing parameter
		#would cause these to be sorted reverse alphabetically in some cases.
		if (decreasing)
		{
			x <- x[order(-x, names(x))];
		}
		else
		{
			x <- x[order(x, names(x))];
		}
  }

  return(x);
}

make.color.scale<- function(aColor, bColor, steps, gradientExp=.5){
	len <- steps-1
	ret <- 0
	for (i in 0:len)
	{
		alpha = (i/len)^gradientExp # doing this yields a better spread of color
		ct = (1 - alpha) * aColor + alpha * bColor
		ret[i+1] <- rgb(ct[1], ct[2], ct[3])
	}
	return(ret);
}

#pretty.print.freq.totals <- function()
#{
#	#Pretty print freq totals
#	
#	#TODO collapse too long words
#	
#	terms <- c('loooong', 'w', 'hello', 'a', 'notaword', 'kotobanai')
#	freqs <- c('1', '2', '3', '4' , '0', '-1')
#	
#	nCol <- 3
#	colWidth <- 23
#	
#	for ( ip1 in 1:ceiling(length(terms)/nCol) )
#	{
#		i <- ip1 - 1 #CURSE YOU 1-based indexing!
#		rowStart <- (nCol*i + 1)
#		rowEnd <- min(rowStart + nCol - 1, length(terms))
#		
#		termRow = terms[rowStart:rowEnd]
#		freqRow = freqs[rowStart:rowEnd]
#		
#		cat(format(termRow,width=colWidth, justify='centre'), sep='|');
#		cat('\n');
#		cat(format(freqRow,width=colWidth, justify='centre'), sep='|');
#		cat('\n')
#		cat(paste(rep("-",nCol*colWidth),collapse=''))
#		cat('\n');
#	}
#}
#pretty.print.freq.totals()

