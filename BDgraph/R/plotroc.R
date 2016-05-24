# function for ROC plot
outRoc = function( G, prob, cut.num )
{
	G[ lower.tri( G, diag = TRUE ) ]       <- 0
	prob[ lower.tri( prob, diag = TRUE ) ] <- 0
	p          = nrow(prob)
	sumEdges   = sum(G)
	sumNoEdges = p * ( p - 1 ) / 2 - sumEdges
	
	tp = c( 1, rep( 0, cut.num ) )
	fp = tp

	cutPoint = ( 0:cut.num ) / cut.num
	
	for ( i in 2:cut.num )
	{
		# checking for cut pints
		estG = 0 * G
		estG[prob > cutPoint[i]] = 1

		tp.all <- ( G != 0 ) * ( estG != 0 ) 
		fp.all <- ( G == 0 ) * ( estG != 0 ) 	
		tp[i]  <- sum( tp.all ) / sumEdges
		fp[i]  <- sum( fp.all ) / sumNoEdges
	}
	
	return( list( tp = tp, fp = fp ) )
}
    
# To plot ROC curve
plotroc = function( sim.obj, bdgraph.obj, bdgraph.obj2 = NULL, bdgraph.obj3 = NULL, cut.num = 20, smooth = FALSE, label = TRUE )
{
    if ( class(sim.obj)     == "sim" ) G = as.matrix( sim.obj $ G ) else G = as.matrix( sim.obj )
    if ( class(bdgraph.obj) == "bdgraph" )
    {
		p_links = bdgraph.obj $ p_links
		if( is.null( p_links ) ) p_links = plinks( bdgraph.obj, round = 10 )
		prob = as.matrix( p_links )
	}else{
		prob = as.matrix( bdgraph.obj )
	}
    
    output = outRoc( G = G, prob = prob, cut.num = cut.num )
    x      = output $ fp
    y      = output $ tp
 	
	if ( smooth == TRUE )
	{
		fit = smooth.spline( x = x, y = y )
		x   = c( 0, fit $ x )
		y   = c( 0, fit $ y )
	}
	
	# par( mar = c( 3.8, 4.2, 1.8, 1 ) )
    plot( NA, type = "l", col = "black", cex.lab = 1.3, cex.main = 2, cex.axis = 1.2,
         main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate", ylim = c(0,1), xlim = c(0,1) )
    points( x = x, y = y, type = "l", col = 1, lty = 1, lw = 2 )
  
    if( !is.null( bdgraph.obj2 ) && is.null( bdgraph.obj3 ) )
    {
		if ( class(bdgraph.obj2)  == "bdgraph" )
		{
			p_links2 = bdgraph.obj2 $ p_links
			if( is.null( p_links2 ) ) p_links2 = plinks( bdgraph.obj2, round = 10 )
			prob2 = as.matrix( p_links2 )
		}else{
			prob2 = as.matrix( bdgraph.obj2 )
		}

        output2 = outRoc( G = G, prob = prob2, cut.num = cut.num )
		x2      = output2 $ fp
		y2      = output2 $ tp

		if ( smooth == TRUE )
		{
			fit2 = smooth.spline( x = x2, y = y2 )
			x2   = c( 0, fit2 $ x )
			y2   = c( 0, fit2 $ y )
		}
		
        points( x = x2, y = y2, type = "l", col = 2, lty = 2, lw = 2 )
		if ( label ) 
			legend( "bottomright", c( "bdgraph.obj", "bdgraph.obj2" ), lty = 1:2, col = 1:2, lwd = c( 2, 2 ), cex = 1.5 )
    }
    
    if( !is.null( bdgraph.obj2 ) && !is.null( bdgraph.obj3 ) )
    {
		if ( class( bdgraph.obj2 )  == "bdgraph" )
		{
			p_links2 = bdgraph.obj2 $ p_links
			if( is.null( p_links2 ) ) p_links2 = plinks( bdgraph.obj2, round = 10 )
			prob2 = as.matrix( p_links2 )
		}else{
			prob2 = as.matrix( bdgraph.obj2 )
		}
   
		if ( class( bdgraph.obj3 )  == "bdgraph" )
		{
			p_links3 = bdgraph.obj3 $ p_links
			if( is.null( p_links3 ) ) p_links3 = plinks( bdgraph.obj3, round = 10 )
			prob3 = as.matrix( p_links3 )
		}else{
			prob2 = as.matrix( bdgraph.obj3 )
		}
   
        output2 = outRoc( G = G, prob = prob2, cut.num = cut.num )
		x2      = output2 $ fp
		y2      = output2 $ tp
   
        output3 = outRoc( G = G, prob = prob3, cut.num = cut.num )
		x3      = output3 $ fp
		y3      = output3 $ tp
   
		if ( smooth == TRUE )
		{
			fit2 = smooth.spline( x = x2, y = y2 )
			x2   = c( 0, fit2 $ x )
			y2   = c( 0, fit2 $ y )
			
			fit3 = smooth.spline( x = x3, y = y3 )
			x3   = c( 0, fit3 $ x )
			y3   = c( 0, fit3 $ y )
		}
		
        points( x = x2, y = y2, type = "l", col = 2, lty = 2, lw = 2 )
        points( x = x3, y = y3, type = "l", col = 3, lty = 3, lw = 2 )
        if ( label ) 
			legend( "bottomright", c( "bdgraph.obj", "bdgraph.obj2", "bdgraph.obj3" ), lty = 1:3, col = 1:3, lwd = c( 2, 2 ), cex = 1.5 )
    }    
}
       
