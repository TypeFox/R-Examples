########################################################################

### Summarization Results Exploration

### Collection of code for examining results of the textreg summarizer

### L. Miratrix

########################################################################




##########################################################################
## Functions to play with output from textreg summarizer
##
## by Luke Miratrix
##########################################################################



#' Call gregexpr on the content of a tm Corpus.
#'
#' Pull out content of a tm corpus and call gregexpr on that content represented
#' as a list of character strings.
#'
#' If 'corpus' is already a character vector, it just calls
#' gregexpr with no fuss (or warning).
#'
#' @return This method gives results exactly as if \code{\link{gregexpr}} were called on the Corpus 
#' represented as a list of strings.
#'
#' @param pattern See gregexpr
#' @param corpus Either a character vector or tm Corpus object.
#' @param ignore.case See gregexpr
#' @param perl See gregexpr
#' @param fixed See gregexpr
#' @param useBytes See gregexpr
#' @import NLP
#' @return See gregexpr.
#' @seealso gregexpr
tm_gregexpr = function( pattern, corpus, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE ) {
	
	if ( is.character( corpus ) ) {
		rs = gregexpr( pattern, corpus, ignore.case, perl, fixed, useBytes )
	} else {
		# TO DO: Is this the fast way, or should it be cast to 
		# character and then called without sapply?
		rs = sapply( corpus, function( minicorp ) {
			gregexpr(  pattern, content(minicorp), ignore.case, perl, fixed, useBytes )
		} )
	}
	rs
}



#' Convert phrases to appropriate search string.
#"
#' @description
#' Will change, e.g., "test * pig+" to appropriate regular 
#' expression to find in the text.
#'
#' @param phrases List of strings denoting the phrases to be
#' searched for.  
make_search_phrases = function( phrases ) {
	if ( length( phrases ) == 0 ) {
		return( c() ) 
	}
	
	x = gsub( "*", "\\w+", phrases, fixed=TRUE )
	x = gsub( "+", "\\w*", x, fixed=TRUE )
	x = gsub( "\\w\\w*", "\\w+", x, fixed=TRUE )  # sad hack
	
	x = paste( "\\b", x, "\\b", sep="" )
	x
}



#' Count number of times documents have a given phrase.
#'
#' Given a list of phrases, count how many documents they appear in
#' and subdivide by positive and negative appearance.
#'
#' This method does not consider multiple counts of phrases within documents.
#' Phrases can have wildcards and stemming notation.  See \code{\link{grab.fragments}}.
#'
#' @seealso grab.fragments
#'
#' @return a dataframe of statistics.  per.pos is the percent of the 
#' documents with the phrase that are positively labeled.  per.tag is
#' the percent of the positively labeled documents that have the 
#' phrase.
#' 
#' @param phrases List of strings
#' @param labeling Vector of +1/0/-1 labels
#' @param corpus A corpus object from tm package
#' @family textregCounting
#' @import NLP
#' @export
#' @examples
#' library( tm )
#' data( bathtub )
#' lbl = meta( bathtub )$meth.chl
#' make.count.table( c("bathtub","strip+", "vapor *"), lbl, bathtub )
make.count.table = function( phrases, labeling, corpus ) {
	
	mcorp = corpus[labeling != 0]
	#print(mcorp)
	#	names(labeling) = unlist( meta( corpus, type="local", "ID" ) )
	
	cnts = sapply( phrases, function(x) { 
		#x = gsub( "*", "\\w*", x, fixed=TRUE )
	    x = make_search_phrases( x )
		#cat( "searching '", x, "'\n", sep="" )
		#pos = tm_index( mcorp, x )
		pos = tm_index(mcorp, FUN = function(doc) any(grep(x, content(doc))))
		#nab = unlist( meta( pos, type="local", "ID" ) )
		hits = labeling[pos]
		c( sum(hits > 0), length(hits) )
	} )
	
	rownames(cnts) = c("n.pos","n")
	res = data.frame( n.pos = cnts[1,], n=cnts[2,] )
	res$per.pos = round(100*res$n.pos/res$n)
	
	ntag = sum( labeling==1 )
	res$per.tag = round( 100 * res$n.pos / ntag )
	res	
}


#' @title Make a table of where phrases appear in a corpus
#'
#' @description
#' Generate a n by p phrase count matrix, with n being number of documents
#' and p being number of phrases:
#'  \\tabular{rrrrr}{
#'   0 \\tab 0 \\tab 0 \\tab 0 \\tab 0 \\cr
#'   1 \\tab 6 \\tab 2 \\tab 0 \\tab 0 \\cr
#'   8 \\tab 0 \\tab 0 \\tab 0 \\tab 0}

#' This is the phrase equivilent of a document-term matrix.
#'
#' @param phrase_list List of strings
#' @param corpus A corpus object from tm package
#' @return a n X p matrix, n being number of documents, p being number of phrases.
#' @family textregCounting
#' @export
#' @examples
#' library( tm )
#' data( bathtub )
#' lbl = meta( bathtub )$meth.chl
#' head( make.phrase.matrix( c("bathtub","strip+", "vapor *"), bathtub ) )
make.phrase.matrix = function( phrase_list, corpus ) {
	
	# TODO Avoid this conversion, but keep vectorization
	# ability of gregexpr?
	# strings = convert.tm.to.character( corpus )
	
	countmatches <- function(pattern,strings){
		#cat( "counting '", pattern, "'\n" )
		matchlist <-  tm_gregexpr(pattern,strings)			
	#	matchlist <-  tm_index(mcorp, FUN = function(doc) any(grep(x, content(doc))))

		matchsimple <- sapply( matchlist, as.vector)
		if ( length(corpus) > 1 ) {
			matchnumbers <- sapply(matchsimple,function(x) if(any(is.na(x))|all(x == -1)) 0 else length(x))
		} else {
			matchnumbers <- if(any(is.na(matchsimple))|all(matchsimple == -1)) 0 else length(matchsimple)
		}
		return(matchnumbers)
	}

	# the following makes the term-doc matrix for the phrases of interest
	phrases = make_search_phrases( phrase_list )
	phrasecount <-sapply(phrases,countmatches,corpus)

	if ( length(corpus) == 1 ) {
		phrasecount = matrix( phrasecount, nrow=1 )
	} 
	colnames(phrasecount) = phrase_list
	phrasecount
}



#' Count phrase appearance.
#'
#' Count number of times a _single_ phrase appears in the corpus
#'
#' @param phrase A string
#' @param corp A corpus object from tm package
#' @family textregCounting
#' @export
#' @examples
#' library( tm )
#' data( bathtub )
#' phrase.count( "bathtub", bathtub )
phrase.count = function( phrase, corp ) {
	stopifnot( length(phrase) == 1 )
	cnts = make.phrase.matrix( c( phrase ), corp )
	cnts[,1]
}







#' Grab all fragments in a corpus with given phrase.
#'
#' @description
#' Search corpus for passed phrase, using some wildcard notation.  Return snippits
#' of text containing this phrase, with a specified number of characters before and
#' after.  This gives context for phrases in documents.
#'
#' Use like this \code{frags = grab.fragments( "israel", bigcorp )}
#'
#' Can take phrases such as 'appl+' which means any word starting with "appl."  
#' Can also take phrases such as "big * city" which consist of any three-word phrase with "big" 
#' as the first word and "city" as the third word.
#'
#' If a pattern matches overlapping phrases, it will return the first but not 
#' the second.
#'
#' @param phrase Phrase to find in corpus
#' @param corp is a tm corpus
#' @param char.before Number of characters of document to pull before phrase to give 
#'        context.
#' @param char.after As above, but trailing characters.  Defaults to char.before value.
#' @param cap.phrase TRUE if the phrase should be put in ALL CAPS.  False if left alone.
#' @param clean True means drop all documents without phrase from list. False means leave
#'    NULLs in the list.
#' @return fragments in corp that have given phrase.List of lists.  First list is 
#'       len(corp) long
#'  with NULL values for documents without phrase, and lists
#'   of phrases for those documents with the phrase
#' @examples
#' library( tm )
#' docs = c( "987654321 test 123456789", "987654321 test test word 123456789", 
#'        "test at start", "a test b", "this is a test", "without the t-word",
#'        "a test for you and a test for me" )
#' corpus <- Corpus(VectorSource(docs))
#' grab.fragments( "test *", corpus, char.before=4, char.after=4 )
#' @export
grab.fragments = function( phrase, corp, char.before = 80,
 			 char.after=char.before, cap.phrase=TRUE, clean=FALSE ) {

	stopifnot( length( phrase ) == 1 )
	stopifnot( nchar(phrase) > 0 )
	
	if ( length(corp) == 0 ) {
		return( c() )	
	}
	phrase.exp= make_search_phrases(phrase)
	phrase.exp2 = paste( "(", phrase.exp, ")", sep="" )		

	#corp = convert.tm.to.character( corp )
	a = tm_gregexpr( phrase.exp, corp )
	res = lapply( 1:length(corp), function( k ) {
		spos = a[[k]]
		topk = attr(a[[k]], "match.length" ) + char.after - 1
		if ( spos[[1]] != -1 ) {
			strs = substring( corp[[ k ]], spos - char.before, spos + topk )
			if ( cap.phrase ) {
				strs = gsub( phrase.exp2, "\\U\\1", strs, perl=TRUE )
			} else {
				strs
			}
		} else {
			c()
		}
	} )
	names(res) = 1:length(corp)
	#npos = res[ sapply(res, length) > 0 ]
	#npos
	
	
	if ( clean ) {
		res.null = sapply( res, is.null )
		res = res[ !res.null ]
	}
	res
}




#' @title Sample fragments of text to contextualize a phrase.
#' 
#' @description
#' Take a phrase, a labeling and a corpus and return text fragments
#' containing that phrase.
#'
#' Grab all phrases and then give sample of N from positive class
#' and N from negative class.  Sampling is to first sample from documents
#' and then sample a random phrase from each of those documents.
#'
#' @export
#' @import tm
#' @param phrases Phrases to examine (a list of strings)
#' @param labeling -- a vector of the same length as the corpus
#' @param corp Corpus object (tm package Corpus object)
#' @param N size of sample to make.
#' @param char.before Number of characters of document to pull before phrase to give 
#'        context.
#' @param char.after As above, but trailing characters.  Defaults to char.before value.
#' @param metainfo -- extra string to add to the printout for clarity if many such
#'     printouts are being generated.
#' @family sample.fragments
#' @examples
#' library( tm )
#' data( bathtub )
#' sample.fragments( "bathtub", meta(bathtub)$meth.chl, bathtub )
sample.fragments = function( phrases, labeling, corp, N=10, char.before=80, char.after=char.before, metainfo=NULL ) {
	
	stopifnot( "Corpus" %in% class(corp) )

	lb = labeling
	stopifnot( !is.null( lb ) )
	stopifnot( length(lb) == length(corp) )
	nP = sum( lb > 0 )
	nN = sum( lb < 0 )

	clusters = sapply( phrases, function( phrase ) {
		cnts = phrase.count( phrase, corp )
		
		selP = which( lb  > 0 & cnts > 0 )
		nfP = length(selP)	
		if ( nfP > N ) {
			selP = sort( sample( selP, N ) )
		}
		selN = which( lb  < 0 & cnts > 0 )
		nfN = length(selN)
		if ( nfN > N ) {
			selN = sort( sample( selN, N ) )
		}
		crpP = corp[selP]
		crpN = corp[selN]
		
		resP = grab.fragments( phrase, crpP, char.before=char.before, char.after=char.after )
		resN = grab.fragments( phrase, crpN, char.before=char.before, char.after=char.after )
		
		resP = sapply( resP, sample, 1 )
		resN = sapply( resN, sample, 1 )	
	
		#resP = sapply( resP, function(x) { gsub( ph, toupper(phrase), x ) } )
		#resN = sapply( resN, function(x) { gsub( ph, toupper(phrase), x ) } )
		
		fragment.sample =  list( phrase=phrase, resP=resP, resN=resN, metainfo=metainfo, nfP=nfP, nP=nP, nfN=nfN, nN=nN, selP = selP, selN=selN )
		class( fragment.sample ) = "fragment.sample"
		fragment.sample
	}, simplify=FALSE )
	
	clusters
}



#' Is object a fragment.sample object?
#' @export
#' @aliases fragment.sample
#' @param x the object to check.
#' @family sample.fragments
is.fragment.sample = function( x ) {
	inherits(x, "fragment.sample")
}






#' Pretty print results of phrase sampling object.
#' @export
#' @param x A fragment.sample object.
#' @param ... No extra options passed.
#' @family sample.fragments
print.fragment.sample = function( x, ... ) {
	
	stopifnot( is.fragment.sample( x ) )
	
	with( x, {
		#ph = make_search_phrases( phrase )
	
		pp = 100 * nfP/nP
		pn = 100 * nfN/nN
		
		if ( !is.null(metainfo) ) {
			cat( sprintf( "\nProfile of Summary Phrase: '%s' - %s", phrase, metainfo ) )
		} else {
			cat( sprintf( "\nProfile of Summary Phrase: '%s'", phrase ) )
		}
		cat( sprintf( "\nPositive: %d/%d = %.2f\nNegative: %d/%d = %.2f\nAppearance of '%s' in positively marked documents:\n", 
					nfP,nP,pp, nfN,nN, pn, phrase ) )
		
		sapply( resP, function(x) { cat( "* ", x, "\n", sep="" ) } )
		
		if ( !is.null(metainfo) ) {
			cat( sprintf( "\nAppearance of '%s' in baseline documents for %s:\n", phrase, metainfo ) )
		} else {
			cat( sprintf( "\nAppearance of '%s' in baseline documents.\n", phrase ) )
		}
		sapply( resN, function(x) { cat( "* ", x, "\n", sep="" ) } )
	} )
	invisible( x )	
}









