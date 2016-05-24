translateScore <-
function( old_score, section ) {

	#data( verbal_conc_table )
	#data( quant_conc_table )
	
	#rename concordenance tables to streamline execution
	verbal = verbal_conc_table
	quant = quant_conc_table
	
	#create tables that include extreme values for the findInterval function
	tmp = rbind( c(0,0,0), get(section), c(801,171,100) )
		
	#convert scores to current scale
	if( section=="verbal" ) {
		score = tmp$new[ findInterval( old_score, tmp$old) ]
	}
	
	if( section=="quant" ) {
		score = tmp$new[ findInterval( old_score, tmp$old) ]
	}
		
	return( score )
}
