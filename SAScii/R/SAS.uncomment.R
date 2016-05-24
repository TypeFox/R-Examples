SAS.uncomment <- 
function( SASinput , starting.comment , ending.comment ){

	#remove /* */
	for ( i in 1:length(SASinput) ){
		
		#test if the line contains a slash_asterisk (or any opening comment character)
		slash_asterisk <- unlist( gregexpr( starting.comment , SASinput[ i ] , fixed = T ) )
		
		#test if the line contains a asterisk_slash (or any closing comment character)
		asterisk_slash <- unlist( gregexpr( ending.comment , SASinput[ i ] , fixed = T ) )
		
		#only if the line contains an slash_asterisk..
		if ( !( -1 %in% slash_asterisk ) ){
			
			#if there's a closing asterisk_slash on that line..
			if ( max( asterisk_slash ) > min( slash_asterisk ) ){
				
				SASinput[ i ] <- sub( substr( SASinput[ i ] , slash_asterisk[1] , asterisk_slash[1] + 1 ) , "" , SASinput[ i ] , fixed = T )
				
				#and re-do the line just in case there's more than one comment!!
				i <- i - 1
				
			} else {
			
				#delete the rest of the line
				SASinput[ i ] <- sub( substr( SASinput[ i ] , slash_asterisk[1] , nchar( SASinput[ i ] ) ) , "" , SASinput[ i ] , fixed = T )
				
				#start a new counter
				j <- i
				
				#keep going until you find a asterisk_slash
				while( max( asterisk_slash ) < 0 ){
					j <- j + 1
					
					#look for asterisk_slash again
					asterisk_slash <- unlist( gregexpr( ending.comment , SASinput[ j ] , fixed = T ) )
					
					#if the asterisk_slash doesn't exist, delete the whole line
					if ( max( asterisk_slash ) < 0 ) SASinput[ j ] <- ""
					#otherwise delete until the asterisk_slash
					else SASinput[ j ] <- sub( substr( SASinput[ j ] , 1 , asterisk_slash[1] + 1 ) , "" , SASinput[ j ] , fixed = T )
				}
			}
		}
	}

	SASinput
}

