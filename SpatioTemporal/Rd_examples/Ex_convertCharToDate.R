##a vector of dates is returned as is
convertCharToDate( seq(as.Date("2012-01-01"),as.Date("2012-01-31"),by=5) )

##if given as chracter vectors Date is returned
convertCharToDate( c("2012-01-01","2012-01-05","2012-01-10","2012-01-12") )

##double is returned as is
convertCharToDate( rnorm(5) )

##other things result in NULL
convertCharToDate( c("a","b","c") )
convertCharToDate( c("2012-01-01", "2012-01-05", "a", "2012-01-12") )
convertCharToDate( c(1,2,3,"d") )
