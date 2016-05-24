pr.curve<-function( scores.class0, scores.class1=scores.class0, weights.class0=NULL, 
		weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}}, sorted = FALSE, curve = FALSE, 
		minStepSize=min(1,ifelse(is.null(weights.class0),1,sum(weights.class0)/100)),
		max.compute=F, min.compute=F, rand.compute=F){
	if(!sorted){
		o0<-order(scores.class0);
		scores.class0<-scores.class0[o0];
		if(!is.null(weights.class0)){
			weights.class0<-weights.class0[o0];
		}
		o1<-order(scores.class1);
		scores.class1<-scores.class1[o1];
		if(!is.null(weights.class1)){
			weights.class1<-weights.class1[o1];
		}
	}
	compute.pr(scores.class0,scores.class1,weights.class0,weights.class1,curve,minStepSize,max.compute,min.compute,rand.compute);
}


roc.curve<-function( scores.class0, scores.class1=scores.class0, weights.class0=NULL, 
		weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}}, sorted = FALSE, curve = FALSE, 
		max.compute=F, min.compute=F, rand.compute=F){
	if(!sorted){
		o0<-order(scores.class0);
		scores.class0<-scores.class0[o0];
		if(!is.null(weights.class0)){
			weights.class0<-weights.class0[o0];
		}
		o1<-order(scores.class1);
		scores.class1<-scores.class1[o1];
		if(!is.null(weights.class1)){
			weights.class1<-weights.class1[o1];
		}
	}
	compute.roc(scores.class0,scores.class1,weights.class0,weights.class1,curve,max.compute,min.compute,rand.compute);
}


check <- function( n, weights ) {
	if( !is.null( weights ) ) {
		if( n != length(weights) ) {
			stop( "The weights must have the same length as the scores." );
		}
		if( sum( weights < 0 ) != 0 ) {
			stop( "The weights must be non-negative." );
		}
	}
}

compute.pr <- function( sorted.scores.class0, sorted.scores.class1=sorted.scores.class0, weights.class0 = NULL, 
						weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}}, curve = FALSE, 
						minStepSize=min(1,ifelse(is.null(weights.class0),1,sum(weights.class0)/100)),
						max.compute=F, min.compute=F, rand.compute=F ){
	
	check( length(sorted.scores.class0), weights.class0 );
	check( length(sorted.scores.class1), weights.class1 );
	
	if( !is.null(sorted.scores.class1) & ( length(sorted.scores.class0) != length(sorted.scores.class1) | 
			suppressWarnings( sum(sorted.scores.class0 != sorted.scores.class1) > 0 ) 
		) & is.null(weights.class0) & is.null(weights.class1) ){
		weights.class0<-c(rep(1,length(sorted.scores.class0)),rep(0,length(sorted.scores.class1)));
		sorted.scores.class0<-c(sorted.scores.class0,sorted.scores.class1);
		o0<-order(sorted.scores.class0);
		sorted.scores.class0<-sorted.scores.class0[o0];
		weights.class0<-weights.class0[o0];
		weights.class1<-1-weights.class0;
		sorted.scores.class1<-sorted.scores.class0;
	}
		
	davis.and.goadrich <- ( length(sorted.scores.class0) == length(sorted.scores.class1) & 
								suppressWarnings( sum( sorted.scores.class0 != sorted.scores.class1 ) == 0 ) & 
								length(weights.class0) == length(weights.class1) &
								suppressWarnings( sum( weights.class0 != (1 - weights.class1) ) == 0 ) &
								sum(weights.class0 != 0 & weights.class0 != 1)==0);
								
		#( is.null( weights.class0 ) | sum( weights.class0 != 1 ) == 0 ) & ( is.null( weights.class1 ) | sum( weights.class1 != 1 ) == 0 );

	i.old <- 0; j.old <- 0; i <- 0; j <- 0; d <- length( sorted.scores.class1 ); m <- length( sorted.scores.class0 );
	help1 <- 0; help2 <- 0;
	auc.GD <- ifelse(davis.and.goadrich,0,NA); auc.integral <- 0; fn <- 0; tn <- 0;
	
	nw0 <- is.null( weights.class0 );
	nw1 <- is.null( weights.class1 );
	
	pos <- ifelse( nw0, m, sum( weights.class0 ) );
	neg <- ifelse( nw1, d, sum( weights.class1 ) );
	
	while( ( j<d ) & sorted.scores.class0[ i + 1 ] > sorted.scores.class1[ j + 1 ] ){
		tn <- tn + ifelse( nw1, 1, weights.class1[ j + 1 ] );
		j <- j + 1;
	}
	p <- c( ( pos - fn ) / pos, ( pos - fn ) / ( pos - fn + neg - tn ), sorted.scores.class0[ i + 1 ] );
	ci <- 1;
	if( curve ){
		list.curve <- create.curve( length( sorted.scores.class0 ) + length( sorted.scores.class1 ) );
		list.curve <- append.to.curve( list.curve, p, ci );
		ci <- ci + 1;
	}else{
		list.curve <- NULL;
	}
	
	unique <- !( j < d & sorted.scores.class0[ i + 1 ] == sorted.scores.class1[ j + 1 ] );
	from.motif <- unique;
	
	while( i< m & j < d ){
		i.old <- i;
		j.old <- j;
		tn.old <- tn;
		fn.old <- fn;
		
		if( !unique || from.motif ){
			while( i + 1 < m & sorted.scores.class0[ i + 1 ] == sorted.scores.class0[ i + 2 ] ){
				fn <- fn + ifelse( nw0, 1, weights.class0[ i + 1 ] );
				i <- i + 1;
			}
			fn <- fn + ifelse( nw0, 1, weights.class0[ i + 1 ] );
			i <- i + 1;
		}						
		if( !unique || !from.motif ){
			while( j + 1 < d & sorted.scores.class1[ j + 1 ] == sorted.scores.class1[ j + 2 ] ){
				tn <- tn + ifelse( nw1, 1, weights.class1[ j + 1 ] );
				j <- j + 1;
			}
			tn <- tn + ifelse( nw1, 1, weights.class1[ j + 1 ] );
			j <- j + 1;
		}
		score<-0;
		if( i < m & j < d ){
			if( sorted.scores.class0[ i + 1 ] == sorted.scores.class1[ j + 1 ] ){
				unique <- F;
				score <- sorted.scores.class0[ i + 1 ];
			}else{
				unique <- T;
				if( sorted.scores.class0[ i + 1 ] < sorted.scores.class1[ j + 1 ] ){
					from.motif <- T;
					score <- sorted.scores.class0[ i + 1 ];
				}else{
					from.motif <- F;
					score <- sorted.scores.class1[ j + 1 ];
				}
			}
		} else {
			if( i < m ) {
				score <- sorted.scores.class0[ i + 1 ];
			} else if( j < d ) {
				score <- sorted.scores.class1[ j + 1 ];
			} else {
				#i=m, j=d
				max = max(sorted.scores.class0[ m ],sorted.scores.class1[ d ]);
				score = max; #+ 0.01*( max - min(sorted.scores.class0[ 1 ],sorted.scores.class1[ 1 ]) ); #max + arbitrary offset
			}
		}
		
		if( fn == fn.old ) {#i == i.old ){
			old.p<-p;
			p <- c( p[ 1 ], ( pos - fn ) / ( pos - fn + neg - tn ), score );
			if(is.nan(p[2])){
				p<-old.p;
			}
			if( curve ){
				list.curve <- append.to.curve( list.curve, p, ci );
				ci <- ci + 1;
			}
		}else{
			p.b <- p[ 1 ];
			p.a <- ( pos - fn ) / pos;
			
			if( davis.and.goadrich ){
				if( i < m | j < d ){# TODO
					prop.term <- ( tn - tn.old ) / ( fn - fn.old );
					h1 <- p[ 1 ]; h2 <- p[ 2 ];
					c <- fn.old + 1;
					help.j <- tn.old + prop.term;
					while( c <= fn ){
						help1 <- (pos - c) / pos;
						help2 <- (pos - c) / ( pos - c + neg - help.j );
						help.j <- help.j + prop.term;
						auc.GD <- auc.GD + ( h2 + help2 ) / 2 * ( h1 - help1 );
						#print(c(1,auc.GD,i=i.v,m=pos,j=j.v,d=neg,c=c,i.old=i.v.o,j.old=j.v.o))
 						h1 <- help1;
 						h2 <- help2;
						c <- c + 1;
					}
				}else{
					auc.GD <- auc.GD + p[ 2 ] * p[ 1 ];
				}
			}
			
			h <- ( tn - tn.old ) / ( fn - fn.old );
			a <- 1 + h;
			b <- ( neg - tn - h * ( pos - fn ) ) / pos;

			if( !isTRUE(all.equal(b, 0)) ){
				auc.integral <- auc.integral + ( p.b - p.a - b / a * ( log( a * p.b + b ) - log( a * p.a + b ) ) ) / a;
			}else{
				auc.integral <- auc.integral + ( p.b - p.a ) / a;
			}
			
			prop.term <- min( ( fn - fn.old ) / ( i - i.old ), minStepSize );
			h <- h*prop.term;
			help.i <- fn.old + prop.term;
			i.old <- i.old + 1;
			help.j <- tn.old + h;
			k=1;
			while( help.i < fn ){
				p <- c( ( pos - help.i ) / pos, ( pos - help.i ) / ( pos - help.i + neg - help.j ), score );#interpolate score?
				if( curve ){
					list.curve <- append.to.curve( list.curve, p, ci );
					ci <- ci + 1;
				}
				k=k+1;
				help.j <- tn.old + k*h;
				help.i <- fn.old + k*prop.term;
			}
			if( p.a != p[ 1 ] ){
				temp <- ( pos - fn ) / ( pos - fn + neg - tn );
				if(is.nan(temp)){
					temp <- p[2];
				}
				p <- c( p.a, temp, score );
				if( curve ){
					list.curve <- append.to.curve( list.curve, p, ci );
					ci <- ci + 1;
				}
			}
		}		
	}
	
	if( i < m ){
		help1 <- 0;
		if( davis.and.goadrich ){
			auc.GD <- auc.GD + p[ 2 ] * ( p[ 1 ] - help1 );
		}
		
		auc.integral <- auc.integral + p[ 2 ] * ( p[ 1 ] - help1 );
		
		p <- c( help1, p[ 2 ], sorted.scores.class0[ i + 1 ] );
		if( curve ){
			list.curve <- append.to.curve( list.curve, p, ci );
			ci <- ci + 1;
		}
	}
	if(curve){
		list.curve<-shrink.curve( list.curve );
	#	list.curve<-rbind(c(list.curve[1,1],list.curve[1,2],min(sorted.scores.class0,sorted.scores.class1)),
	#					  list.curve,
	#					  c(list.curve[nrow(list.curve),1],list.curve[nrow(list.curve),2],max(sorted.scores.class0,sorted.scores.class1)))
	}
	res<-list( type = "PR", auc.integral = auc.integral, auc.davis.goadrich = auc.GD, curve=list.curve );
	
	if(max.compute){
		scores0<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores0<-rep(1,length(sorted.scores.class0));
		}else{
			scores0<-weights.class0;
		}
		scores1<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores1<-rep(0,length(sorted.scores.class1));
		}else{
			scores1<-weights.class0;
		}
		
		max.res<-pr.curve( scores.class0=scores0, scores.class1=scores1,weights.class0=weights.class0,
					weights.class1=weights.class1,curve=curve,minStepSize=minStepSize);
		res<-c(res,list(max=max.res));
	}
	
	if(min.compute){
		scores0<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores0<-rep(0,length(sorted.scores.class0));
		}else{
			scores0<-(-weights.class0);
		}
		scores1<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores1<-rep(1,length(sorted.scores.class1));
		}else{
			scores1<-(-weights.class0);
		}
		
		min.res<-pr.curve( scores.class0=scores0, scores.class1=scores1,weights.class0=weights.class0,
					weights.class1=weights.class1,curve=curve,minStepSize=minStepSize);
		res<-c(res,list(min=min.res));
	}
	if(rand.compute){
		rand.auc<-NULL;
		if(is.null(weights.class0)){
			rand.auc<-length(sorted.scores.class0)/(length(sorted.scores.class0)+length(sorted.scores.class1));	
		}else{
			rand.auc<-sum(weights.class0)/sum(weights.class0+weights.class1);
		}
		rand.curve<-create.curve( 2 );
		rand.curve<-append.to.curve( rand.curve, c(0,rand.auc,0), 1 );
		rand.curve<-append.to.curve( rand.curve, c(1,rand.auc,0), 2 );
		rand.result<-list( type = "PR", auc.integral = rand.auc, auc.davis.goadrich = rand.auc, curve=rand.curve );
		class(rand.result)<-"PRROC";
		
		res<-c(res,list(rand=rand.result));
	}
	
	class(res)<-"PRROC";
	res
}

compute.roc<-function( sorted.scores.class0, sorted.scores.class1=sorted.scores.class0, weights.class0 = NULL, 
					   weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}}, curve = FALSE,
					   max.compute=F, min.compute=F, rand.compute=F){
	
	if( !is.null(sorted.scores.class1) & ( length(sorted.scores.class0) != length(sorted.scores.class1) | 
		 suppressWarnings( sum(sorted.scores.class0 != sorted.scores.class1) > 0 ) 
		) & is.null(weights.class0) & is.null(weights.class1) ){
			weights.class0<-c(rep(1,length(sorted.scores.class0)),rep(0,length(sorted.scores.class1)));
			sorted.scores.class0<-c(sorted.scores.class0,sorted.scores.class1);
			o0<-order(sorted.scores.class0);
			sorted.scores.class0<-sorted.scores.class0[o0];
			weights.class0<-weights.class0[o0];
			weights.class1<-1-weights.class0;
			sorted.scores.class1<-sorted.scores.class0;
	}
	
	i <- 0; j <- 0; d <- length( sorted.scores.class1 ); m <- length( sorted.scores.class0 );
	fn <- 0; tn <- 0;
	
	nw0 <- is.null( weights.class0 );
	nw1 <- is.null( weights.class1 );
	
	pos <- ifelse( nw0, m, sum( weights.class0 ) );
	neg <- ifelse( nw1, d, sum( weights.class1 ) );
	
	erg <- 0;
	ci <- 1;
	p <- c( 1, 1, min(sorted.scores.class0,sorted.scores.class1) );
	if( curve ){
		list.curve <- create.curve( length( sorted.scores.class0 ) + length( sorted.scores.class1 ) );
		list.curve <- append.to.curve( list.curve, p, ci );
		ci <- ci + 1;
	}else{
		list.curve <- NULL;
	}
	
	unique <- F; from.motif <- F;
	if( sorted.scores.class0[ i + 1 ] == sorted.scores.class1[ j + 1 ] ){
		unique <- F;
	}else{
		unique <- T;
		if( sorted.scores.class0[ i + 1 ] < sorted.scores.class1[ j + 1 ] ){
			from.motif <- T;
		}else{
			from.motif <- F;
		}
	}
	
	while( i < m & j < d ){
		score <- 0;
		if( unique ){
			if( from.motif ){
				while( i < m & sorted.scores.class0[ i + 1 ] < sorted.scores.class1[ j + 1 ] ){
					fn <- fn + ifelse( nw0, 1, weights.class0[ i + 1 ] );
					score <- sorted.scores.class0[ i + 1 ];
					i <- i + 1;
				}
				
			}else{
				while( j < d & sorted.scores.class0[ i + 1 ] > sorted.scores.class1[ j + 1 ] ){
					tn <- tn + ifelse( nw1, 1, weights.class1[ j + 1 ] );
					score <- sorted.scores.class1[ j + 1 ];
					j <- j + 1;
				}
				#score <- sorted.scores.class0[ i + 1 ];
			}
		}else{
			while( i + 1 < m & sorted.scores.class0[ i + 1 ] == sorted.scores.class0[ i + 2 ] ){
				fn <- fn + ifelse( nw0, 1, weights.class0[ i + 1 ] );
				i <- i + 1;
			}
			while( j + 1 < d & sorted.scores.class1[ j + 1 ] == sorted.scores.class1[ j + 2 ]){
				tn <- tn + ifelse( nw1, 1, weights.class1[ j + 1 ] );
				j <- j + 1;
			}
			fn <- fn + ifelse( nw0, 1, weights.class0[ i + 1 ] );
			tn <- tn + ifelse( nw1, 1, weights.class1[ j + 1 ] );
			i <- i + 1;
			j <- j + 1;
			score <- sorted.scores.class0[ i ];
		}
		
		help1 <- ( neg - tn ) / neg;
		help2 <- ( pos - fn ) / pos;
		erg <- erg + ( p[ 2 ] + help2 ) / 2 * ( p[ 1 ] - help1 );
		p <- c( help1, help2, score );
		if(curve){
			list.curve <- append.to.curve( list.curve, p, ci );
			ci <- ci + 1;
		}
		
		if( i < m & j < d ){
			if( sorted.scores.class0[ i + 1 ] == sorted.scores.class1[ j + 1 ] ){
				unique <- F;
			}else{
				unique <- T;
				if( sorted.scores.class0[ i + 1 ] < sorted.scores.class1[ j + 1 ] ){
					from.motif <- T;
				}else{
					from.motif <- F;
				}
			}
		}
	}
	
	if(curve){
		p <- c( 0, 0, max( sorted.scores.class0, sorted.scores.class1 ) );
		list.curve <- append.to.curve( list.curve, p, ci );
		ci <- ci + 1;
		list.curve<-shrink.curve( list.curve );
		list.curve<-rbind(c(list.curve[1,1],list.curve[1,2],min(sorted.scores.class0,sorted.scores.class1)),
						  list.curve,
						  c(list.curve[nrow(list.curve),1],list.curve[nrow(list.curve),2],max(sorted.scores.class0,sorted.scores.class1)))
	}
	res<-list( type = "ROC", auc = erg, curve=list.curve );
	
	
	if(max.compute){
		scores0<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores0<-rep(1,length(sorted.scores.class0));
		}else{
			scores0<-weights.class0;
		}
		scores1<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores1<-rep(0,length(sorted.scores.class1));
		}else{
			scores1<-weights.class0;
		}
		
		max.res<-roc.curve( scores.class0=scores0, scores.class1=scores1,weights.class0=weights.class0,
							 weights.class1=weights.class1,curve=curve);
		res<-c(res,list(max=max.res));
	}
	
	if(min.compute){
		scores0<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores0<-rep(0,length(sorted.scores.class0));
		}else{
			scores0<-(-weights.class0);
		}
		scores1<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores1<-rep(1,length(sorted.scores.class1));
		}else{
			scores1<-(-weights.class0);
		}
		
		min.res<-roc.curve( scores.class0=scores0, scores.class1=scores1,weights.class0=weights.class0,
							 weights.class1=weights.class1,curve=curve);
		res<-c(res,list(min=min.res));
	}
	if(rand.compute){
		rand.auc<-0.5;
		rand.curve<-create.curve( 2 );
		rand.curve<-append.to.curve( rand.curve, c(0,0,0), 1 );
		rand.curve<-append.to.curve( rand.curve, c(1,1,0), 2 );
		rand.result<-list( type = "ROC", auc=rand.auc, curve=rand.curve );
		class(rand.result)<-"PRROC";
		
		res<-c(res,list(rand=rand.result));
	}
	
	class(res)<-"PRROC";
	res	
}

shrink.curve <- function( curve ){
	if( is.null( curve ) ){
		curve;
	}else{
		curve[ !is.na( curve[ , 1 ] ), ];
	}
}

create.curve <- function( n ){
	m <- matrix( NA, nrow=n, ncol=3 );
	m
}

append.to.curve <- function( curve, p, row ){
	if( row>=nrow( curve ) ){
		curve2 <- matrix( NA, nrow=nrow( curve ) * 2, ncol=3 );
		curve2[ 1:nrow( curve ), ] <- curve;
		curve <- curve2;
	}
	curve[ row, ] <- p;
#	print(c(row,p))
#	if(is.nan(p[2])){
#		traceback(0)
#	}
	curve
}

print.PRROC<-function(x,...){
	if(x$type == "PR"){
		cat("\n  Precision-recall curve\n");
		cat("\n    Area under curve (Integral):\n");
		cat("    ",x$auc.integral,"\n");
		if( !is.null(x$max) & !is.null(x$min) ){
			cat("\n    Relative area under curve (Integral):\n");
			cat("    ",(x$auc.integral - x$min$auc.integral)/(x$max$auc.integral-x$min$auc.integral),"\n");
		}
		cat("\n    Area under curve (Davis & Goadrich):\n");
		if(!is.null(x$auc.davis.goadrich) & !is.na(x$auc.davis.goadrich)){
			cat("    ",x$auc.davis.goadrich,"\n");
			if( !is.null(x$max) & !is.null(x$min) ){
				cat("\n    Relative area under curves (Davis & Goadrich):\n");
				cat("    ",(x$auc.davis.goadrich - x$min$auc.davis.goadrich)/(x$max$auc.davis.goadrich-x$min$auc.davis.goadrich),"\n");
			}
		}else{
			cat("    cannot be computed for weighted data\n");
		}
		
	}else{
		cat("\n  ROC curve\n");
		cat("\n    Area under curve:\n");
		cat("    ",x$auc,"\n");
		if( !is.null(x$max) & !is.null(x$min) ){
			cat("\n    Relative area under curve:\n");
			cat("    ",(x$auc - x$min$auc)/(x$max$auc-x$min$auc),"\n");
		}
	}
	
	if(!is.null(x$curve)){
		cat("\n    Curve for scores from ",min(x$curve[,3])," to ",max(x$curve[,3]),"\n");
		cat("    ( can be plotted with plot(x) )\n\n");
	}else{
		cat("\n    Curve not computed ( can be done by using curve=TRUE )\n");
	}
	
	if(!is.null(x$max)){
		cat("\n\n    Maximum AUC:\n");
		if(x$type == "PR"){
			cat("    ",x$max$auc.integral," ",x$max$auc.davis.goadrich,"\n");
		}else{
			cat("    ",x$max$auc,"\n");
		}
	}
	
	if(!is.null(x$min)){
		cat("\n\n    Minimum AUC:\n");
		if(x$type == "PR"){
			cat("    ",x$min$auc.integral," ",x$min$auc.davis.goadrich,"\n");
		}else{
			cat("    ",x$min$auc,"\n");
		}
	}
	
	if(!is.null(x$rand)){
		cat("\n\n    AUC of a random classifier:\n");
		if(x$type == "PR"){
			cat("    ",x$rand$auc.integral," ",x$rand$auc.davis.goadrich,"\n");
		}else{
			cat("    ",x$rand$auc,"\n");
		}
	}
}


plot.PRROC<-function(x, xlim=c(0,1), ylim=c(0,1), auc.main=TRUE, auc.type=c("integral","davis.goadrich"), 
					 legend=ifelse(is.logical(color) & color==TRUE,4,NA), xlab=NULL, ylab=NULL, main=NULL, color=TRUE, lwd=3, 
					 add=FALSE, scale.color=hsv(h=seq(0,1,length=100)*0.8, s=1, v=1), 
					 max.plot = FALSE, min.plot = FALSE, rand.plot = FALSE, fill.area = (max.plot & min.plot),
					 maxminrand.col = grey(0.5), fill.color = grey(0.95),
					 ...){
	auc.type<-match.arg(auc.type);
	if(is.null(x$curve)){
		stop("Curve is NULL. Use curve=T in pr.curve or roc.curve to obtain one.");
	}
	if(ncol(x$curve) != 3){
		stop("Curve has wrong dimension");
	}
	if(is.null(xlab)){
		my.xlab<-ifelse(x$type=="PR","Recall","FPR");
	}else{
		my.xlab<-xlab;
	}
	if(is.null(ylab)){
		my.ylab<-ifelse(x$type=="PR","Precision","Sensitivity");
	}else{
		my.ylab<-ylab;
	}
	
	if(is.null(main)){
		my.main<-paste(x$type," curve",sep="",collapse="");
	}else{
		my.main<-main;
	}
	if(auc.main){
		my.main<-paste(my.main,"\nAUC = ",format(ifelse(x$type=="PR",ifelse(auc.type=="integral",x$auc.integral,x$auc.davis.goadrich),x$auc)),sep="",collapse="");
	}
 
	
	max.curve<-NULL;
	if(!is.null(x$max) & !is.null(x$max$curve)){
		max.curve<-x$max$curve;
	}
	min.curve<-NULL;
	if(!is.null(x$min) & !is.null(x$min$curve)){
		min.curve<-x$min$curve;
	}
	rand.curve<-NULL;
	if(!is.null(x$rand) & !is.null(x$rand$curve)){
		rand.curve<-x$rand$curve;
	}
	
	x<-x$curve;
	
	cols<-1;
	segment=F;
	plotscale.color=F;
	if( is.logical(color) ){
		if(color){
			min<-min(x[,3]);
			max<-max(x[,3]);
	
			cols<-getColor( scale.color, x[,3], min, max );
			plotscale.color=T;
			segment=T;
		}else{
			cols<-1;
			segment<-F;
		}
	}else {
		cols<-color;
		segment<-F;
	}
	
	if(!add & !is.na(legend) & (is.numeric(legend) | suppressWarnings(legend==TRUE)) & plotscale.color ){
		if(is.logical(legend)){
			legend<-4;
		}
		m<-NULL;widths<-rep(1,2);heights<-rep(1,2)
		if(legend == 1){
			m<-matrix(c(1,2),nrow=2);
			heights<-c(4,lcm(2));
		}else if(legend==2){
			m<-matrix(c(2,1),nrow=1);
			widths=c(lcm(2.5),4);
		}else if(legend==3){
			m<-matrix(c(2,1),nrow=2);
			heights=c(lcm(2),4);
		}else{
			m<-matrix(c(1,2),nrow=1);
			widths=c(4,lcm(2.5));
		}
		layout(mat = m,widths = widths,heights = heights);
		
	}#else if(!add){
	#	layout(1);
	#}
	
	if(!add){
		plot(0,xlim=xlim,ylim=ylim,col=0,xlab=my.xlab,ylab=my.ylab,main=my.main,...);
	}
	
	if( !add ){
		if( fill.area & !is.null(max.curve) & !is.null(min.curve)){
			xs<-c(min.curve[,1],max.curve[nrow(max.curve):1,1],min.curve[1,1]);
			ys<-c(min.curve[,2],max.curve[nrow(max.curve):1,2],min.curve[1,2]);
			polygon( x = xs, y = ys, density = -1, border = NA, col = fill.color );
		}
		
		if(max.plot & !is.null(max.curve)){
			lines(max.curve[,1],max.curve[,2],col=maxminrand.col, lty="dashed", ...);
		}
		
		if(min.plot & !is.null(min.curve)){
			lines(min.curve[,1],min.curve[,2],col=maxminrand.col, lty="dotted", ...);
		}
		
		if(rand.plot & !is.null(rand.curve)){
			lines(rand.curve[,1],rand.curve[,2],col=maxminrand.col, lty="dotdash", ...);
		}
	}
	
	d=nrow(x);
	if( segment ) {
		segments( x[1:(d-1),1], x[1:(d-1),2], x[2:d,1], x[2:d,2], col=cols, lwd=lwd, ...);
	} else {
		lines( x[,1], x[,2], col=cols, lwd=lwd, ...);
	}
	
	if(!add & legend & !is.numeric(color) & color == TRUE){
		scale<-seq( min, max, length = 100 );
		cols<-getColor( scale.color, scale, min, max );
		bak<-par("mar");
		on.exit(par(mar=bak));
		if(legend==2 | legend==4){
			if(legend==4){par(mar=c(5,1,4,2)+0.1);}else{par(mar=c(5,2,4,1)+0.1);}
			image(c(1),scale,matrix(scale,nrow=1),col=cols,xlab="",ylab="",axes=F)
		}else{
			if(legend==1){par(mar=c(2,4,0,2)+0.1);}else{par(mar=c(0,4,2,2)+0.1);}
			image(scale,c(1),matrix(scale,ncol=1),col=cols,xlab="",ylab="",axes=F)
		}
		axis(legend)
		layout(1)
	}
	

}

getColor <- function( scale, x, min=min(x), max=max(x) )  {
	return( scale[round(1 + (length(scale)-1) * (x - min)/(max-min))] );
}
