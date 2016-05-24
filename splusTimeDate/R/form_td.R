format.timeDate <- 
  list("1/3/1998" = list( output = "%m/%d/%Y", input = "%m/%d/%Y"),
       "3/1/1998" = list( output = "%d/%m/%Y", input = "%d/%m/%Y"),
       "01/03/1998" = list( output = "%02m/%02d/%Y", input = "%m/%d/%Y"),
       "03/01/1998" = list( output = "%02d/%02m/%Y", input = "%d/%m/%Y"),
       " 1/ 3/1998" = list( output = "%2m/%2d/%Y", input = "%m/%d/%Y"),
       " 3/ 1/1998" = list( output = "%2d/%2m/%Y", input = "%d/%m/%Y"),
       "1/3/98" = list( output = "%m/%d/%y", input = "%m/%d/%y"),
       "3/1/98" = list( output = "%d/%m/%y", input = "%d/%m/%y"),
       "1-3-1998" = list( output = "%m-%d-%Y", input = "%m-%d-%Y"),
       "3-1-1998" = list( output = "%d-%m-%Y", input = "%d-%m-%Y"),
       "01-03-1998" = list( output = "%02m-%02d-%Y", input = "%m-%d-%Y"),
       "03-01-1998" = list( output = "%02d-%02m-%Y", input = "%d-%m-%Y"),
       "Jan 3, 1998" = list( output = "%b %d, %Y", input = "%m %d, %Y"),
       "3 Jan 1998" = list( output = "%d %b %Y", input = "%d %m %Y"),
       "3Jan1998" = list( output = "%d%b%Y", input = "%d %m %Y"),
       "3Jan98" = list( output = "%d%b%C", input = "%d %m %y"),
       "Jan 03, 1998" = list( output = "%b %02d, %Y", input = "%m %d, %Y"),
       "03 Jan 1998" = list( output = "%02d %b %Y", input = "%d %m %Y"),
       "03Jan1998" = list( output = "%02d%b%Y", input = "%d %m %Y"),
       "January 3, 1998" = list( output = "%B %d, %Y", input = "%m %d, %Y"),
       "3 January 1998" = list( output = "%d %B %Y", input = "%d %m %Y"),
       "Sat Jan 3, 1998" = list( output = "%a %b %d, %Y", input = "%w %m %d, %Y"),
       "Sat 3 Jan 1998" = list( output = "%a %d %b %Y", input = "%w %d %m %Y"),
       "19980103" = list( output = "%Y%02m%02d", input = "%4Y%2m%2d" ),
       "01/03/1998 02:04 PM" = list( output = "%02m/%02d/%Y %02I:%02M %p", 
	 input = "%m/%d/%Y %H:%M %p"),
       "03/01/1998 14:04" = list( output = "%02d/%02m/%Y %02H:%02M", 
	 input = "%d/%m/%Y %H:%M"),
       "Jan 03, 1998 02:04 PM" = list( output = "%b %02d, %Y %02I:%02M %p", 
	 input = "%m %d, %Y %H:%M %p"),
       "03 Jan 1998 14:04" = list( output = "%02d %b %Y %02H:%02M", 
	 input = "%d %m %Y %H:%M"),
       "01/03/1998 02:04:32 PM (PST)" = 
       list(output = "%02m/%02d/%Y %02I:%02M:%02S %p (%z)", 
	 input = "%m/%d/%Y %H:%M:%S %p (%:)Z)"),
       "03/01/1998 14:04:32 (PST)" = 
       list(output = "%02d/%02m/%Y %02H:%02M:%02S (%z)", 
	 input = "%d/%m/%Y %H:%M:%S (%:)Z)"),
       "Jan 03, 1998 02:04:32 PM (PST)" = 
       list(output = "%b %02d, %Y %02I:%02M:%02S %p (%z)", 
	 input = "%m %d, %Y %H:%M:%S %p (%:)Z)"),
       "03 Jan 1998 14:04:32 (PST)" = 
       list(output = "%02d %b %Y %02H:%02M:%02S (%z)", 
	 input = "%d %m %Y %H:%M:%S (%:)Z)")
       )
       
  
