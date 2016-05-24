#  Get raw data ready to be 'summarized' by textreg C++ function
#
#  Based on Robin's version based on Luke's earlier code
#  12-9-2012



#' @title Clean text and get it ready for textreg.
#'
#' @description 
#' Changes multiline documents to single line.  Strips extra whitespace and punctuation.
#' Changes digits to 'X's.  Non-alpha characters converted to spaces.
#'
#' @import tm
#' @export
#' @param bigcorp  A tm Corpus object.
#' @examples
#' library( tm )
#' txt = c( "thhis s! and bonkus  4:33pm and Jan 3, 2015. ", 
#'          "   big    space\n     dawg-ness?")
#' a <- clean.text( Corpus( VectorSource( txt ) ) )
#' a[[1]]
clean.text = function( bigcorp ) {

  convert.numbers.and.simplify <- function(x) {
	         x2 = paste(x, collapse = " ")
	         x2 = gsub("[0-9]", "X", x2)
	         x2 = gsub("[^a-zA-Z]", " ", x2)
	         PlainTextDocument(x2, id = meta(x,"id"), language = meta( x, "language") )
	}
	
	
	bigcorp <- tm_map(bigcorp, content_transformer (tolower))
	bigcorp <- tm_map(bigcorp, convert.numbers.and.simplify)
	bigcorp <- tm_map(bigcorp, stripWhitespace)
	#bigcorp <- tm_map(bigcorp, removePunctuation)
	#bigcorp <- tm_map(bigcorp, removeNumbers)
	
	bigcorp
}


#' @title Save corpus to text (and RData) file.
#'
#' @description
#' Small utility to save a corpus to a text file (and RData file) for ease of use.
#'
#' It is possibly recommended to pass a filename to the C++ function \code{\link{textreg}}
#' rather than the entire corpus for
#' large text since I believe it will otherwise copy over everything due to the coder's (my) poor
#' understanding of how RCpp converts objects.
#'
#' @import tm
#' @param bigcorp A tm Corpus object.
#' @param filename The first part of the filename.  A rda and txt extension will be appended to the two generated files.
#'
#' @export
save.corpus.to.files = function( bigcorp, filename="corpus" ) {

	# save our data for later use.  What is disk space, really?
	save(bigcorp, file=paste( filename, ".rda", sep="" ) )

	# first the documents
	out_file = file( paste( filename, ".txt", sep="" ) )
	open( out_file, open="w" )
	ig = sapply( bigcorp, function( x ) { writeLines( content(x), out_file ) } )
	close(out_file)
}


