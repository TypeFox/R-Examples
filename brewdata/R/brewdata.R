brewdata <-
function( years=2015, term="F", degree="phd", focus="statistics", 
	resolution=10, map=FALSE ) {
	
	#error handling
	tryCatch({
		message("Checking user-specified parameters for errors or typos...")

		#user specified variables
		TERM = toupper(term)
		YEARS = unique( substr(years, 3, 4) )
		DEGREE = tolower(degree)
		FOCUS = tolower(focus)
		},
		error=function(cond) {
			message( 
			"One or more brewdata inputs is invalid."
			)
			stop("brewdata stopped.", call.=FALSE )
		}
	)

	#Scape GradeCafe Results Search and fill the 'data' data.frame
	data = getGradCafeData( YEARS, TERM, DEGREE, FOCUS )
	
	#Fill the 'results' data.frame with the parsed info. The results in the
	#'data' data.frame are unusable strings. Call the parseResults method to 
	#gather self-reported personal statistics (gpa, gre, decision date).
	results = data.frame()
	for( results_i in data$results ){
		results = rbind( results, parseResults( results_i ) )
	}
	
	#Translate all pre-2011 scores to the current scale
	V=results$gre_v; Q=results$gre_q
	I_oldv = which( V>170 ); I_oldq = which( Q>170 )
	results$gre_v[ I_oldv ] = translateScore( V[I_oldv], "verbal" )
	results$gre_q[ I_oldq ] = translateScore( Q[I_oldq], "quant" )
	
	#Find the percentiles for the scores
	v_pct = findScorePercentile( results$gre_v, "verbal" )
	q_pct = findScorePercentile( results$gre_q, "quant" )
	aw_pct = findScorePercentile( results$gre_aw, "writing" )
	percentiles = cbind( v_pct, q_pct, aw_pct )
	
	#Replace unusable strings with the parsed data, then parse the names
	school_name = parseSchools( as.matrix( data[,1] ), resolution, map )
	
	dec_stat = data.frame( results[,8], data[,-c(1:2,4)] )
	colnames( dec_stat ) = c( "decision", "status" )
	out = cbind( school_name, dec_stat, results[,1:4], percentiles, 
		results[,5:7] )
	return( out )
}
