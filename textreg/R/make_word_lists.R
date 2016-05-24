

# Add phrase counts, etc., to a word table
annotate.table = function( a, dat ) {
	
	#words = do.call( rbind, sapply( dat, function(d) { d$model[-1,] }, simplify=FALSE ) )
	# NOTE FOR ABOVE: ngram rather than phrase is var name.  Fix
	words = do.call( rbind, sapply( dat, function(d) { reformat.textreg.model( d ) }, simplify=FALSE ) )
	words = words[ !duplicated( words$phrase ), ]
	words
	words$beta = words$Z = NULL
	a = as.data.frame(a)
	a$phrase = rownames(a)
	#ords = a$ngram
	a = merge( a, words, by="phrase", all.x=TRUE, sort=FALSE )
	a
	
}



# Utility function: Make small dataframe for making word lists
textreg.to.df = function( textreg.res, RunID = 1 ) {   
    data.frame( RunID = RunID, 
                word = textreg.res$model$ngram,
                weight = textreg.res$model$beta / textreg.res$model$Z )
}





#' Collate multiple regression runs.
#'
#' This method makes a table of several regression runs side by side.
#' The table has rows being phrases and the
#' columns being the regression runs.  A number is usually the 
#' weight found for that word at that window.
#' If multiple runs have the same phrase, row will
#' have multiple entries. 
#' 
#' Method will also order rows based on calculated importance of phrases.
#' Multiple ways of ordering are possible, via the \code{method} argument.
#'
#' Finally, the table can be annotated with descriptive statistics of the phrases.
#'
#' Warning: this method DOES NOT flip negative weight words (so negative weight usually look
#' less imporant in the ordering).
#'
#' See the bathtub vignette for an example of this method.
#'
#' @return If annotate = true, a dataframe with each column corresponding to an textreg.result 
#'   object (and possibly extra columns about phrases).  Otherwise a matrix of the word scores.
#' @export
#' @param result.list  List of mix of textreg.result objects and dataframes with two columns of
#'        "word" and "weight".  (The latter is for merging lists from other regression packages.)
#' @param model.names Names of the textreg.result objects
#' @param M maximum number of words to keep
#' @param topic String A name for the topic
#' @param method Different ways to sort the phrases.  'word' means make a list of words.
#' @param annotate  Add summary statistics to table such as phrase counts, etc.
make.list.table = function (result.list, model.names = names(result.list), M = 100, 
                            topic = "Summary Collection", method = c("rank", "weight", "count", "word"),
                            annotate=TRUE) 
{
    if (!is.list(result.list)) {
        stop("Need a list of result objects")
    }
    stopifnot(length(result.list) > 1)
    df.list = result.list
    for ( i in seq_along(df.list) ) {
        if ( is.textreg.result( df.list[[i]] ) ) {
            df.list[[i]] = textreg.to.df( df.list[[i]] )
        } 
        stopifnot( nrow(df.list[[i]] ) > 0 )
        df.list[[i]]$RunID = i
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
        stop( "need to install plyr package to run make.list.table or other list generation methods" )
    }
    dat = plyr::rbind.fill( df.list )
    
    Xs = unique(dat$word)
    first = tapply(dat$RunID, dat$word, min)
    first = tapply(dat$RunID, dat$word, min)
    method = match.arg(method)
    if (method == "rank") {
        dat$rnk = dat$weight
        for (i in unique(dat$RunID)) {
            runi = dat$RunID == i
            dat$rnk[runi] = rank(abs(dat$weight[runi]))
        }
        a = tapply(dat$rnk * sign(dat$weight), list(dat$word, 
                                                    dat$RunID), sum)
    }
    else if (method == "weight" || method == "word") {
        a = tapply(abs(dat$weight), list(dat$word, dat$RunID), 
                   sum)
    }
    else if (method == "count") {
        a = tapply(sign(dat$weight), list(dat$word, dat$RunID), 
                   length)
    }
    else {
        stop(gettextf("Method %s not understood as valid mode in make.list.table", 
                      method))
    }
    stopifnot(all(names(first) == rownames(a)))
    keeps = rownames(a) != "*intercept*"
    
    if ( sum( keeps ) == 1 ) {
        names = rownames(a)
        b = matrix( a[keeps,], nrow=1 )
        rownames(b) = names[keeps]
    } else {
        a = a[keeps, ]

        # sort the list
        a[is.na(a)] = 0
        b2 = as.list(data.frame(abs(a)))
        ord = do.call(order, c(b2, decreasing = TRUE))
        b = a[ord, ]
        b[b == 0] = NA
    }
    
    if (method == "word") {
        words = rownames(b)
        b = matrix(as.character(b), nrow = nrow(b))
        for (i in 1:ncol(b)) {
            b[, i] = ifelse(is.na(b[, i]), "", words)
        }
        rownames(b) = words
    }
    colnames(b) = model.names

    if ( annotate ) {
        ngo = sapply( result.list, is.textreg.result )
        if ( sum( ngo ) > 0 ) {
            b = annotate.table( b, result.list[ngo] )
        }
        
        attr(b, "num.models" ) = length( result.list )	
    }
    
    attr(b, "topic" ) = topic
    attr(b, "method" ) = method
    #class( b ) = c( class(b), "list.table" )
    b
}


old.make.list.table = function( result.list, model.names = names(result.list), M = 100, topic="Summary Collection",
				method=c("rank","weight","count","word"),
				annotate=TRUE ) {
	
	if ( !is.list( result.list ) ) {
		stop( "Need a list of result objects" )
	}
	stopifnot( length( result.list ) > 1 )
	
	if ( length( result.list ) == 1 ) {
		# doesn't work yet.
		runids = rep( 1, nrow(result.list[[1]]$model) )
		words = result.list[[1]]$model$ngram 
		weights = result.list[[1]]$model$beta / result.list[[1]]$model$Z 	
	} else {
		runids = do.call( c, sapply( 1:length(result.list), function( d ) { rep( d, nrow(result.list[[d]]$model) ) } ) )
		words = do.call( c, sapply( result.list, function( d ) { d$model$ngram } ) )
		weights = do.call( c, sapply( result.list, function( d ) { d$model$beta / d$model$Z } ) )
	}
	dat = data.frame( RunID = runids, word=words, weight=weights )
	
	Xs = unique( dat$word )
	
	# in which window did each word first occur?
	first = tapply( dat$RunID, dat$word, min )
	
	#	sapply( 1:length(dat), function( d ) { rep( d, nrow(dat[[d]]$nrow) ) } )
	#wordlists = sapply( dat, function( d ) { d$model$ngram[-1] }, simplify=FALSE 
	#Xs = unique( unlist( wordlists ) )
			
	#cat( "# Unique words found = ", length(Xs), " / ", length(dat$word ), "\n" )
	
	# in which window did each word first occur?
	first = tapply( dat$RunID, dat$word, min )
	
	# How to pick which words to keep
	method = match.arg(method)
	if ( method=="rank" ) {
		dat$rnk = dat$weight
		for ( i in unique(dat$RunID) ) {
			runi = dat$RunID==i
			dat$rnk[runi] = rank(abs(dat$weight[runi]))
		}
		a = tapply( dat$rnk*sign(dat$weight), list(dat$word, dat$RunID), sum )
	} else if ( method == "weight" || method=="word" ) {
		a = tapply( abs(dat$weight), list(dat$word, dat$RunID), sum )
	} else if ( method == "count" ) {
		# make grid of number of times a word appeared across all windows 
		a = tapply( sign(dat$weight), list(dat$word, dat$RunID), length )
	# } else if ( method == "word" ) {
		# a = tapply( dat$word, list(dat$word, dat$RunID), paste )
		
	} else {
		stop( gettextf( "Method %s not understood as valid mode in make.list.table", method ) )
	}
	stopifnot( all( names(first) == rownames(a) ) )

	a[is.na(a)] = 0
	a = a[ rownames(a) != "*intercept*", ]
	#a = a / max(a)
	
	# order words by the columns (time) going from earliest time to last
	# to get words ordered by first appearance, and then by
	# strength of first appearance
	b2 = as.list( data.frame( abs(a) ) ) # a[1:nrow(a),] )
	ord = do.call( order, c( b2, decreasing=TRUE ) )
	b = a[ ord, ] 
		
	b[ b == 0 ] = NA

	if ( method=="word" ) {
		words = rownames(b)
		b = matrix( as.character(b), nrow=nrow(b) )
		for ( i in 1:ncol(b) ) {
			b[,i] = ifelse( is.na( b[,i] ), "", words )
		}
		rownames(b) = words
	}
	colnames(b) = model.names	

	if ( annotate ) {
		b = annotate.table( b, result.list )
		attr(b, "num.models" ) = length( result.list )	
	}

	attr(b, "topic" ) = topic
	attr(b, "method" ) = method
    #class( b ) = c( class(b), "list.table" )
	b
}







#' Graphic showing multiple word lists side-by-side.
#' 
#' This method basically makes a visual plot of a list table (which you call first).
#' 
#' @export
#'
#' @param model.list Matrix (or data.frame) from the make.list.table call.
#' @seealso make.list.table
#' @param M is the max number of words to show in chart
#' @param linespace Where to space
#' @param ytick Put y tick marks
#' @param dates Dates to put on bottom
#' @param main Main title
#' @param xlab Label for x-axis
#' @param xaxt Plot an x-axis (see par)
#' @param mar Margin of plot (see par)
#' @param color.breaks  Cut-points (like on a histogram) defining the different color levels.
#' @param color.ramp List of colors to use from lowest value (potentially negative weights) to highest.  If both color.breaks and color.ramp passed, color.breaks is list one longer than color.ramp. 
#' @param ... Extra arguments for the core image() call that plots the word weights.
list.table.chart = function( model.list, M=100, linespace=4, ytick=NULL, 
		dates=NULL, 
		main = paste( "Word Appearance for ", attr(model.list,"topic"), "\n(Method: ", attr(model.list,"method"),")", sep=""),
		xlab="Model",
		mar = c(3,5,2.5,0.1),
		xaxt="y",
		color.breaks = NULL,
		color.ramp = NULL,
		... ) {
	

	if ( is.data.frame(model.list) ) {
		nms = model.list$phrase
		b = model.list[2:(1+attr(model.list,"num.models"))]
		b = as.matrix(b)
		rownames(b) = nms
	} else {
		b = model.list
	}

	# flip so it plots first row at top
    b = b[ nrow(b):1, ]
	
	# heuristic for grabbing "best" words if only plotting
    # M words total
    if ( nrow( b ) > M ) {
        bm = rank( apply(abs(b), 1, max, na.rm=TRUE ), na.last=FALSE )
	    bm2 = rank( apply(abs(b), 1, mean, na.rm=TRUE ), na.last=FALSE )
	    bmm = pmax( bm, bm2 )
	    b[b==0] = NA
	
	    # take only M words deemed best.
	    b = b[ bmm >= quantile( bmm, 1-M/length(bmm)), ]
    }
    
	 # maximum weight to plot
     max = max(abs(b), na.rm=TRUE)
     if ( min( b, na.rm=TRUE ) < 0 ) {
         min = -max
     } else {
         min = 0
     }
     
     # make breakpoints in weights if needed.
     if ( is.null( color.breaks ) ) {
        if ( is.null( color.ramp ) ) {
            color.breaks <- seq(-max, max, length.out=11 )
        } else {
            color.breaks <- seq(-max, max, length.out=length(color.ramp) + 1)
        }
     } else {
         if ( !any( 0 == color.breaks ) ) {
             warning( "No cut-point at 0 for color.breaks.  Hard to seperate negative and positive words." )
         }
     }
     
     if ( is.null( color.ramp ) ) { # always do color
		 sw = which.min( abs( color.breaks ) )
		 Ln = length( color.breaks ) - 1
		 reds = greens = rep( 0, Ln )
		 
		 if ( sw > 1 ) {
		     reds[1:(sw-1)] = seq( 1.0, 0.5, length=(sw-1) )
         }
		 if ( sw < Ln ) {
		     greens[sw:Ln] = seq( 0.5, 1.0, length=1+Ln-sw )
		 }
         color.ramp <- rgb( reds, greens, 0 ) 
	 }

	 # Plot Data Map
	 par(mar = mar)
	
	image( x=1:ncol(b),
			y = 1:nrow(b), zlim=c(min, max),
			t(b), 
			xlab=xlab, ylab="", 
			xaxt="n", yaxt="n", bty="n",
			col=color.ramp,
			breaks=color.breaks,
			main=main, ... )

	# horizontal lines
	for ( i in seq(1,nrow(b), by=linespace ) ) {
		abline( h=i, lty=3 )
	}
	
	# tickmarks along bottom
	if ( xaxt=="y" ) {
	    ticks = seq(1,ncol(b),by=1)
    	axis(BELOW<-1, at=ticks, labels=colnames(b), 
    					cex.axis=0.7, las = 2)
    	if ( !is.null( dates ) ) {
    		ticks = seq(1,ncol(b)+1,by=1)
    		axis(BELOW<-1, at=ticks - 0.5, labels=dates, col="grey", 
    					cex.axis=0.5, las = 2, tick=FALSE)
    	}
	}
	
	# words along left
	axis(LEFT <-2, at=1:nrow(b), labels=rownames(b),
 						las= HORIZONTAL<-1,
 						cex.axis=0.5)
	
	invisible( b )
}



