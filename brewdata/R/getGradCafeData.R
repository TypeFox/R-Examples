getGradCafeData <-
function( years, term, degree, focus ) {
	
	#create usable grep string, e.g. "F10|F11|F14"
	app_cycles = paste( term, years,"|", sep="", collapse="")
	tmp = nchar(app_cycles)-1
	app_cycles = substring(app_cycles,0,tmp); rm(tmp)
	
	#create usable grep string, e.g. "2011|2013|2015"
	app_years = paste( "20", years,"|", sep="", collapse="")
	tmp = nchar(app_years)-1
	app_years = substring(app_years,0,tmp); rm(tmp)

	#set the url and find the total number of pages in the dataset
	url = paste("http://www.thegradcafe.com/survey/index.php?q=",focus,
		"&t=a&pp=250&o=&p=",sep="")
	max_i = getMaxPages( url )
	
	#initialize loop paramaters
	keep_scrolling = fetch = TRUE
	df = data.frame(); i=1; record_ct=0
	
	#build dataset
	while( fetch )
	{
		cat( "Downloading page", i, "from GradCafe Results Search.\n" )
		page_i = paste(url,i, sep="")
		raw_i = readHTMLTable( page_i )[[1]]
		df_i = data.frame( 
			raw_i$Institution,  
			raw_i$"Decision & Date", 
			raw_i$St1,
			raw_i$"Program (Season)"
			)[-1,]
	
		colnames(df_i) = c("school","results","status","program")
	
		#Next couple lines force the loop to advance until we visit the app year.
		#The keep_scrolling variable evaluates to TRUE until we see the first
		#records. Once we see that first record, the app_year_term variable
		#takes over. The loop will continue until we get to the end of the
		#dataset we want. Note that initalizing record_ct to 0 and using max
		#record_ct guarantees the loop terminates after we've created the
		#dataset we want--and not before.
		check_cycle = grepl( min(years), df_i$program ) 
		check_year = grepl( min(app_years), df_i$results )
		record_ct = max( sum( check_cycle, check_year ), record_ct )
		keep_scrolling = ifelse( record_ct==0, TRUE, FALSE )
		
		#test whether or not we should break out of the loop
		keep_data = check_cycle | check_year
		if( keep_scrolling | sum(keep_data)>0 & i<min( 500, max_i+1 ) ) 
		{
			df = rbind( df, df_i[ keep_data,] )
			i=i+1
		} else {
			fetch=FALSE
			data = df[ grep( degree, tolower( df$program ) ), ]     
		}
	}
	return( data )
}