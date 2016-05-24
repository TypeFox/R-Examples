##
## Support code for manipulating the stuff handed back from the C++ program
##


#' Make matrix of where phrases appear in corpus.
#'
#' Construct a $n X p$ matrix of appearances for selected phrases out of textreg object.
#' $n$ is the number of documents, $p$ is the number of phrases selected in the result object `rules.'
#' 
#' @export
#' @param rules Either a textreg.result object or the rules list from such an object.
#' @param n  (Optional) If giving a rules list, the number of documents in corpus.
phrase.matrix = function( rules, n ) {

	if ( is.textreg.result( rules ) ) {
		n = rules$notes$n
		rules = rules$rules
	}
	
	mat = matrix( 0, nrow=n, ncol=length(rules) )
	
	colnames(mat) = sapply( rules, function(x) { x$ngram } )
	
	for ( i in 1:length(rules) ) {
		rl = rules[[i]]$support + 1
		#rl = subset( rl, rl < 0 )
		mat[ rl, i ] = rules[[i]]$weight
	}
	
	mat[ , "*intercept*"] = 1
	mat
}



#' Predict labeling with the selected phrases.
#' 
#' Given raw text and a textreg model, predict the labeling by counting appearance of 
#' relevant phrases in text and then multiplying these counts by the beta vector
#' associated with the textreg object.  Just like linear regression.
#'
#' @export
#' @param object A textreg.result object
#' @param new.text If you want to predict for new text, pass it along.
#' @param return.matrix TRUE means hand back the phrase appearance pattern matrix.
#' @param ... Nothing can be passed extra.
#' @return Vector of predictions (numbers).
#'
#' @examples
#' res = textreg( c( "", "", "A", "A" ), c( -1, -1, 1, 1 ), 
#'       C=1, Lq=1, convergence.threshold=0.00000001, verbosity=0 )
#' predict( res )
#' predict( res, new.text=c("A B C A") )
predict.textreg.result = function( object, new.text= NULL, return.matrix=FALSE, ... ) {
	stopifnot( is.textreg.result( object ) )
	model = object
	
	if ( !is.null( new.text ) ) {
		keyphrase.mat = make.phrase.matrix( model$model$ngram, new.text )
		keyphrase.mat[ , "*intercept*" ] = 1
	} else {
		keyphrase.mat = phrase.matrix( model )
	}
	model = model$model
	
	kp = sweep( keyphrase.mat, 2, model$Z, FUN="/" )
	
	rsp = as.numeric( kp %*% model$beta )
	
	if ( return.matrix ) {
		attr( rsp, "keyphrase.matrix" ) <- keyphrase.mat
	}
	rsp
}


#' Calculate total loss of model (Squared hinge loss).
#' 
#' Calculate the loss for a model in predicting the -1/+1 labeling.
#' If new text and labeling given, then calc loss on the new text and labeling.
#' This can be useful for cross validation and train-test splits.
#'
#' @export
#' @param model.blob The model returned from \code{\link{textreg}}
#' @param new.text New text (string or tm Corpus) to predict labeling for
#' @param new.labeling Labeling to go with new text.
#' @param loss Type of loss to calc for.
#' @return Three numbers: total loss, loss from prediction, loss from penalty term
#' @examples
#' data( testCorpora )
#' res = textreg( c( "", "", "A", "A" ), c( -1, -1, 1, 1 ), C=1, Lq=1,
#'           convergence.threshold=0.00000001, verbosity=0 )
#' calc.loss( res )
#' calc.loss( res, new.text=c("A B C A"), new.labeling=c(1) )
calc.loss = function( model.blob, new.text=NULL, new.labeling=NULL, loss=c( "square.hinge", "square", "hinge") ) {

	loss = match.arg( loss )

	model = model.blob$model
	#	k.mat = phrase.matrix( model.blob$rules, model.blob$notes$n )
	
	pd = predict( model.blob, new.text )
	
	if ( is.null( new.labeling ) ) {
		if ( !is.null( new.text ) ) {
			stop( "New text without new labeling" )
		}
		new.labeling = model.blob$labeling
	}
	
	if ( loss =="square.hinge" ) {
		loss = sum( pmax( (1 - new.labeling*pd), 0 )^2 )
	} else if ( loss == "square" ) {
		loss = sum( (pd-new.labeling)^2 )
	} else if ( loss == "hinge" ) {
		loss = sum( pmax( (1 - pd* new.labeling), 0 ) )
	}
		
	pen = model.blob$notes$C * sum( abs( model$beta[ -1 ] ) )
	
	c( tot.loss=loss+pen, loss=loss, penalty=pen )
}


#' Clean up output from textreg.
#' 
#' Calculate some useful statistics (percents, etc) and return as dataframe.
#'
#' @export
#' @param model The model returned from \code{\link{textreg}}
#' @param short True if the output should be abbrviated for easy consumption.
#' @return Dataframe with statistics on the terms in the model
#' @family textreg.result
reformat.textreg.model = function( model, short=TRUE ) {
	stopifnot( is.textreg.result( model ) )
	
	npos = model$model$posCount[[1]]
	nneg = model$model$negCount[[1]]
	
	mod = model$model
	mod$per = mod$posCount / mod$totalDocs
	mod$perPos = mod$posCount / npos
	mod$perNeg = mod$negCount / nneg
	
	#mod$perWord = mod$posWordCount / mod$support
	
	if ( !short ) {
		mod[ c( "ngram", "beta", "Z", "support", "totalDocs", "posCount", "negCount", "per", "perPos", "perNeg" ) ]
		#, "posWordCount", "negWordCount", "perWord" ) ]
	} else {
		mod = mod[ c("ngram", "support", "totalDocs", "posCount", "per", "perPos") ]
		names(mod) = c("phrase", "num.phrase", "num.reports", "num.tag", "per.tag", "per.phrase" )
		mod$per.tag = round( 100 * mod$per.tag )
		mod$per.phrase = round( 100 * mod$per.phrase )
		mod = mod[ order( mod$per.phrase, decreasing=TRUE ), ]
		mod
	}
	
}

