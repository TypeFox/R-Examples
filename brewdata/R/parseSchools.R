parseSchools <-
function( original_name, resolution=10, map=FALSE ) {	

	#data( error_key )
	#data( slang_key )
	#data( dict )
	
	#find & remove bad text values that could crash the function
	junk = original_name %in% showNonASCII( original_name )
	original_name[junk] = "NA"
	parsed = as.data.frame( original_name )
	colnames(parsed) = "names"
	parsed = raw = sapply( parsed$names, tolower )

	#do some basic cleaning using regular expressions
	dirty_char = c(" {0,}- {0,}", "[.]", ",", 
		"un[ iv][a-z]* {0,}| {0,}u {1,}| u$",
		"the | at | [oi]f| in ", "texas$",
		"[(][a-z ]*[)]", " {2,}", "^ | $")
	clean_char = c(" ", "", " ", 
		" univ ", 
		"", "texas austin","", " ", "" )
	tmp = 1:length( clean_char )
	for(i in tmp ) parsed = gsub( dirty_char[i], clean_char[i], parsed )
	rm(tmp)
	
	#replace common naming variations for the same parsed	
	tmp=1:dim( error_key )[1]
	for(i in tmp ) parsed[grepl( error_key[i,1], parsed )] = error_key[i,2]
	rm(tmp)
		
	#try replacing abreviations and slang versions of the university names
	s = which( nchar(parsed)<8 )
	slang_nms = parsed[s]
	replace_bool = slang_nms %in% slang_key[,1]
	replace_index = match( slang_nms[ replace_bool ], slang_key[,1] )
	parsed[s][replace_bool] = slang_key[replace_index,2]
		
	#if all else fails, then replace ambiguous names with the closest match
	tmp = parsed[!junk]
	dist_mat = stringdistmatrix( tmp, dict, method='lcs')
	
	min_dist = apply( dist_mat, 1, min )
	min_index = apply( dist_mat, 1, which.min )
	
	#mn => match names assignment function that uses the resolution parameter
	#resolution=0 => only exact matchs are returned <=> sparse names
	#resolution>20 => course matchs are allowed <=> better aggregated stats
	R = resolution
	mn = function(j) ifelse( min_dist[j]<R, dict[ min_index[j],1 ], tmp[j] )
	tmp = sapply(1:length(tmp), mn )
	parsed[!junk] = tmp; rm(tmp)
	
	#which names changed
	#out = unique( cbind( raw, parsed )[order(parsed),] )
	#changed = which( out[,1]!= out[,2] )
	
	map_out = cbind( parsed, raw)
	colnames( map_out ) = c("school_name","original_name")
	if( map ) out = map_out else out = parsed
		
	return( out )
}
