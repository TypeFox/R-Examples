parse.SAScii <- 
function( sas_ri , beginline = 1 , lrecl = NULL ){

	##and now actually pull the entire file into R from the FTP site, line-by-line
	SASinput <- readLines( sas_ri )
	
	#remove all tab characters
	SASinput <- gsub( "\t" , " " , SASinput )
	
	##if the SAS code includes more than one INPUT, start at the user-specified beginline
	SASinput <- SASinput[ beginline:length(SASinput) ]
	
	SASinput <- toupper( SASinput )

	#remote all /* and */ from the code
	SASinput <- SAS.uncomment( SASinput , "/*" , "*/" )

	#remote all * and ; from the code
	SASinput <- SAS.uncomment( SASinput , "*" , ";" )

	#find the first line with the word INPUT in it, which is where the ASCII variable locations occur.

	#lines that start with input
	firstline <- grep("INPUT",SASinput)[1]

	#find the first semicolon ending that input line
	a<-grep(";",toupper(SASinput))
	lastline<-min(a[a>firstline])

	#isolate the Fixed-Width File (FWF) input lines
	FWFlines<-SASinput[firstline:lastline]

	#remove the word input from the first line
	input_word <- unlist( gregexpr( "INPUT" , FWFlines[1] , fixed = T ) )
	FWFlines[ 1 ] <- substr( FWFlines[ 1 ] , input_word + 5 , nchar(FWFlines[ 1 ]) )

	#remove the semicolon from the last line
	semicolon <- unlist( gregexpr( ";" , FWFlines[length(FWFlines)] , fixed = T ) )
	FWFlines[length(FWFlines)] <- substr( FWFlines[length(FWFlines)] , 1 , semicolon - 1 )

	#put a space in front of all dollar signs
	for ( i in 1:length( FWFlines ) ) FWFlines[ i ] <- gsub( "$" , " $ " , FWFlines[ i ] , fixed = T )
	for ( i in 1:length( FWFlines ) ) FWFlines[ i ] <- gsub( "-" , " - " , FWFlines[ i ] , fixed = T )

	#remove all fully-blank lines
	FWFlines <- FWFlines[ which( gsub( " " , "" , FWFlines ) != "" ) ]

	#break apart all FWF lines
	z<-strsplit(FWFlines," ",perl=T)

	#initiate massive character vector
	SAS.input.lines <- NULL

	for ( i in 1:length( z ) ){
		#throw out all splits that are empty
		z[[i]] <- gsub( "-" , " " , z[[i]] )
		z[[i]] <- z[[i]][ which( gsub( " " , "" , z[[i]] ) != "" ) ]
		
		#and then combine everything into one huge character vector
		SAS.input.lines <- c( SAS.input.lines , z[[i]] )
	}

	##create FWF structure file (x)
	x <- data.frame(NULL)

	i <- j <- 1

	#pull out the second, third, and fourth elements after input line
	elements_2_4 <- SAS.input.lines[ 2:4 ]
	#remove dollar signs for this test, they don't count
	elements_2_4 <- elements_2_4[ elements_2_4 != "$" ]
	
	#figure out from the first line if the numbers are widths of each column
	#or if they're the actual location on the file
	#look at the first line-- how many non $ numerics are there before you hit the second variable name?
	widths_not_places <- 
		#if there was a dollar sign, the length will be two..
		(length( elements_2_4 ) == 2 & 
		is.na( as.numeric( as.character( elements_2_4[2] ) ) ) ) 
	
	
	#look for any @ symbols in the input lines!
	if ( sum(grepl("@",SAS.input.lines))>0 ){
	#if the input lines appear to contain @START VARNAME FORMAT then use this block:


		#cycle through entire character vector
		while ( i < length( SAS.input.lines ) ){

			start.point <- as.numeric( gsub( "@" , "" , SAS.input.lines[ i ] , fixed = T ) )
			
			#skip the first time:
			if ( i > 1 ){
				#if there's room between the current start point and the previous width, add some empty space
				if ( x[ j - 1 , "start" ] + x[ j - 1 , "width" ] < start.point ){
					#this creates a negative width
					x[ j , "width" ] <- ( x[ j - 1 , "start" ] + x[ j - 1 , "width" ] ) - start.point
					j <- j + 1
				}
			}
		
			#set first word to variable name
			x[ j , "start" ] <- start.point
			x[ j , "varname" ] <- SAS.input.lines[ i + 1 ]
			
			#if there's a dollar sign between second word and the format, record that this is of type character
			if ( SAS.input.lines[ i + 2 ] == "$" ){
				x[ j , "char" ] <- T
				i <- i + 1
			} else x[ j , "char" ] <- F
			
			#remove leading f's and char's
			for ( k in c("F","CHAR") ){
				SAS.input.lines[ i + 2 ] <- gsub( k , "" , SAS.input.lines[ i + 2 ] , fixed = T )
			}
				
			#if the length has a period, split it
			if ( grepl( "." , SAS.input.lines[ i + 2 ] , fixed = T ) ){
				period <- unlist( gregexpr( "." , SAS.input.lines[ i + 2 ] , fixed = T ) )
				x[ j , "width" ] <- as.numeric( substr( SAS.input.lines[ i + 2 ] , 1 , period - 1 ) )
				divisor <- substr( SAS.input.lines[ i + 2 ] , period + 1 , nchar( SAS.input.lines[ i + 2 ] ) )
			} else {
				x[ j , "width" ] <- as.numeric( SAS.input.lines[ i + 2 ] )
				divisor <- ""
			}
			
			if ( divisor != "" ) {
				x[ j , "divisor" ] <- 1 / 10^as.numeric(divisor)
			} else x[ j , "divisor" ] <- 1
			
			i <- i + 3
			j <- j + 1
		}

	} else if ( widths_not_places ) {
	#if the input lines appear to contain VARNAME LENGTH then use this block:		

		#cycle through entire character vector
		while ( i < length( SAS.input.lines ) ){

			#set first word to variable name
			x[ j , "varname" ] <- SAS.input.lines[ i ]
			
			#if there's a dollar sign between first word and the first number, record that this is of type character
			if ( SAS.input.lines[ i + 1 ] == "$" ){
				x[ j , "width" ] <- as.numeric( SAS.input.lines[ i + 2 ] )
				x[ j , "char" ] <- T
				i <- i + 3
			
				#otherwise record that it's type numeric
			} else {
				x[ j , "width" ] <- as.numeric( SAS.input.lines[ i + 1 ] )
				x[ j , "char" ] <- F
				i <- i + 2
			}
			
			#search for a divisor
			if ( grepl( "." , SAS.input.lines[ i ] , fixed = T ) ){
					
				period <- unlist( gregexpr( "." , SAS.input.lines[ i ] , fixed = T ) )
				
				divisor <- substr( SAS.input.lines[ i ] , period + 1 , nchar( SAS.input.lines[ i ] ) )
				
				x[ j , "divisor" ] <- 1 / 10^as.numeric(divisor)
				i <- i + 1
			} else x[ j , "divisor" ] <- 1
			
			#jump to the next row of x
			j <- j + 1
		}
		
	#if the input lines appear to contain VARNAME #START - #END then use this block:	
	} else {
		#cycle through entire character vector
		while ( i < length( SAS.input.lines ) ){

			#set first word to variable name
			x[ j , "varname" ] <- SAS.input.lines[ i ]
			
			#if there's a dollar sign between first word and the first number, record that this is of type character
			if ( SAS.input.lines[ i + 1 ] == "$" ){
				x[ j , "start" ] <- SAS.input.lines[ i + 2 ]
				
				#check if the width was one number or two..
				if ( 
					#if it isn't numeric..
					is.na( as.numeric( SAS.input.lines[ i + 3 ] ) ) |
					#or if it contains a period..
					grepl( "." , SAS.input.lines[ i + 3 ] , fixed = T )
					){
						#then it's moved too far because the width was a single digit..
						x[ j , "end" ] <- x[ j , "start" ]
						
						#and it should move back one
						i <- i - 1
					} else {
						#otherwise, if it's a character string, 
						x[ j , "end" ] <- SAS.input.lines[ i + 3 ]
					}
				
				x[ j , "char" ] <- T
				i <- i + 4
			
			#otherwise record that it's type numeric
			} else {
				x[ j , "start" ] <- SAS.input.lines[ i + 1 ]

				#check if the width was one number or two..
				if ( 
					#if it isn't numeric..
					is.na( as.numeric( SAS.input.lines[ i + 2 ] ) )  |
					#or if it contains a period..
					grepl( "." , SAS.input.lines[ i + 2 ] , fixed = T )
					){
						#then it's moved too far because the width was a single digit..
						x[ j , "end" ] <- x[ j , "start" ]
						
						#and it should move back one
						i <- i - 1
					} else {
						#otherwise, if it's a character string, 
						x[ j , "end" ] <- SAS.input.lines[ i + 2 ]
					}

				x[ j , "char" ] <- F
				i <- i + 3
			}
			
			#search for a divisor
			if ( grepl( "." , SAS.input.lines[ i ] , fixed = T ) ){
				
				period <- unlist( gregexpr( "." , SAS.input.lines[ i ] , fixed = T ) )
				
				divisor <- substr( SAS.input.lines[ i ] , period + 1 , nchar( SAS.input.lines[ i ] ) )
				
				x[ j , "divisor" ] <- 1 / 10^as.numeric(divisor)
				
				i <- i + 1
			} else x[ j , "divisor" ] <- 1
			
			#BUT if we're on the second row already..
			if ( j > 1 ) {
			
				#IF current row's start > previous row's end + 1
				if ( as.numeric( x[ j , 'start' ] ) > as.numeric( x[ j - 1 , 'end' ] ) + 1 ){
					#then you need to add in some blank space!
					x <-
						rbind( 
							x[ 1:(j-1) , ] ,
							NA ,
							x[ j , ]
						)
					
					#add one to j, since you've added a row
					j <- j + 1
					
					#and add a negative
					x[ j - 1 , 'start' ] <- as.numeric( x[ j - 2 , 'end' ] ) + 1
					x[ j - 1 , 'end' ] <- as.numeric( x[ j , 'start' ] ) - 1
				}
			}
			
			#jump to the next row of x
			j <- j + 1
		}
		
		#the width should be the end position minus the beginning position, plus one
		x <- 
			transform(
				x ,
				width = as.numeric(end) - as.numeric(start) + 1 
			)
		
		#if there's no variable name, it should be a negative.
		x[ is.na( x[ , 'varname' ] ) , 'width' ] <- ( -1 * x[ is.na( x[ , 'varname' ] ) , 'width' ] )
		
	}

	#limit to only four columns
	x <- x[ , c("varname","width","char","divisor") ]

	#finally, if the final logical record length is specified by the user..
	if ( !is.null( lrecl ) ){
		
		#if it's the same as the sum of the widths already in x, do nothing (specifying it was unnecessary)
		
		#if it's less than the sum of the absolute values of current widths..
		if ( lrecl < sum( abs( x$width ) ) ) stop ( "specified logical record length (lrecl) parameter is shorter than the SAS columns constructed" )
		
		#if it's more than the sum of the absolute value of the current widths..
		if ( lrecl > sum( abs( x$width ) ) ){
		
			#blank space containing the difference should be added onto the tail of x
			length.of.blank.record.to.add.to.end <- ( lrecl - sum( abs( x$width ) ) )
		
			x[ nrow( x ) + 1 , 'width' ] <- -length.of.blank.record.to.add.to.end
		}
		

	}
	
	x
}

