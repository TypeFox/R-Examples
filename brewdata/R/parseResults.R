parseResults <-
function(result) {
	
	#find the admissions decision
	act = grepl( "accepted", tolower(result) )				#accepted
	wtl = grepl( "wait", tolower(result) )					#wait listed
	rej = grepl( "rejected", tolower(result) )				#rejected
	inv = grepl( "interview", tolower(result) )				#interview invite
	def = ifelse( act|wtl|rej|inv, FALSE, TRUE )			#default
	dec_ind = c( act, wtl, rej, inv, def )
	dec_key = c("A", "W", "R", "I", "O" )
	decision = dec_key[ dec_ind ]
	
	#break text feild into date and stats chunks
	chunk = unlist( strsplit( as.character(result), "GPA:") )
	
	#Parse the date. These rules are pretty stable since the dates are 
	#reported consistantly across the dataset.
	date_ck = unlist( strsplit( chunk[1], " on " ) )[2]
	date_str = unlist( strsplit( date_ck, " ") )[1:3]
	day = as.numeric( date_str[1] )
	mon = date_str[2]
	yr = as.numeric( date_str[3] )
	date_act = data.frame( mon, day, yr )
	
	#Parsing the self-reported strings is much trickier than dates. There are
	#some weird entries in this chunk, so the failGracefully method helps make
	#the most of what's available. 
	failGracefully = function(j) 
	{
		trash = as.data.frame( t( vector( mode="numeric", length=4 ) ) )
		colnames( trash ) = c( "gpa", "gre_v", "gre_q", "gre_aw" )
		return( trash )
	}
	
	#Buckle-up. We're about to try parsing the personal stats string...
	#Start by testing whether or not there's any data to available to evaluate.
	#If there is, try to extract the stats. Otherwise, just fill the stats row
	#with zeros. 
	if( !is.na( chunk[2] ) )
	{
		#Take a swing at the stat chunk. There's no guarantee we can actually
		#use the results. We'll do our best with the tryCatch method. 
		stat_ck = chunk[2]
		stat_str = unlist( strsplit( stat_ck, "(V/Q/W)") )[2]
		score = unlist( strsplit( stat_str, " ") )[2]
		gre = unlist( strsplit( score, "/") )
		gre_v = gre[1]
		gre_q = gre[2]
		gre_aw = unlist( strsplit( gre[3], "GRE") )
		gpa = unlist( strsplit( stat_ck[1], "GRE") )[1]
		
		#attempt to parse any stats reported
		stat_act = tryCatch(
			{
				#Try converting the stats to numbers we can use. If it works,
				#we'll be done and return the output. Otherwise, proceed down
				#the line. 
				v = as.numeric( gre_v )
				q = as.numeric( gre_q )
				aw = as.numeric( gre_aw )
				gpa = as.numeric( gpa )	
				success = data.frame( gpa, v, q, aw )
				colnames(success) = c( "gpa", "gre_v", "gre_q", "gre_aw" )
				success
			},
			error=function(cond) {
				#If we get an error, lets...
				failGracefully(1)
			},
			warning=function(cond) {
				#If we get an warning, lets...
				failGracefully(1)
			}
		)
	} else {
		#We end-up here if the stats chunk is completely empty. 
		stat_act = failGracefully(1)
	}
	
	#Return the OUTPUTS to fetchData and exit.
	return( cbind( stat_act, date_act, decision ) )
}
