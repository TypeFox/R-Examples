findScorePercentile <-
function( score, section ) {
	
	#data( svq_score_table )
	#data( saw_score_table = NULL )
	
	#create tables that include extreme values for the findInterval function
	svq_tmp = rbind( c(0,0,0), svq_score_table, c(801,171,100) )
	saw_tmp = rbind( c(-1,-1), saw_score_table, c(7,100) )
	
	#lookup percent of scores below 'score'
	if( section=="verbal" ) {
		pct = svq_tmp$v[ findInterval( score, svq_tmp$score) ]
	}
	if( section=="quant" ){
		pct = svq_tmp$q[ findInterval( score, svq_tmp$score) ]
	}
	if( section=="writing" ){
		pct = saw_tmp$aw[ findInterval( score, saw_tmp$score) ]
	}
		
	return( pct )
}
