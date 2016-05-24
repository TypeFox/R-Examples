require("TSmisc")
require("tfplot")

##################################################
##################################################

#### Data from Reserve Bank of Australia  ########

##################################################
##################################################

#  http://www.rba.gov.au/statistics/tables/index.html
#  http://www.rba.gov.au/statistics/tables/xls/d03hist.xls  #has money
#  http://www.rba.gov.au/statistics/tables/xls/g10hist.xls  #has GDP
#  http://www.rba.gov.au/statistics/tables/xls/f13hist.xls  #has interest rate (international too)
#  http://www.rba.gov.au/statistics/tables/xls/i01hist.xls  #has Int real GDP
#  http://www.rba.gov.au/statistics/tables/xls/i02hist.xls  #has Int CPI


####  Australian Money ####
  #  test file copied Nov. 29, 2010 from 
  #  http://www.rba.gov.au/statistics/tables/xls/d03hist.xls  
  
  testfile <- system.file("xlsExampleData/d03hist.xls", package = "TSmisc")

  con1 <- TSconnect("xls", dbname=testfile,
          map=list(ids  =list(i=11,     j="B:Q"), 
	           data =list(i=12:627, j="B:Q"), 
	           dates=list(i=12:627, j="A"),
                   names=list(i=4:7,    j="B:Q"), 
		   description = NULL,
		   tsrepresentation = function(data,dates){
		       ts(data,start=c(1959,7), frequency=12)}))

  con1u <- TSconnect("xls",
          dbname="http://www.rba.gov.au/statistics/tables/xls/d03hist.xls",
          map=list(ids  =list(i=11,     j="B:Q"), 
	           data =list(i=12:627, j="B:Q"), 
	           dates=list(i=12:627, j="A"),
                   names=list(i=4:7,    j="B:Q"), 
		   description = NULL,
		   tsrepresentation = function(data,dates){
		       ts(data,start=c(1959,7), frequency=12)}))

  z <- TSget("DMACN", con1)
  tfplot(z)

  z <- TSget(c("DMAM1N", "DMAM3N"), con1)
   
  con2 <- TSconnect("xls", dbname=testfile,
          map=list(ids  =list(i=11,     j="B:Q"), 
	           data =list(i=12:627, j="B:Q"), 
	           dates=list(i=12:627, j="A"),
                   names=list(i=4:7,    j="B:Q"), 
		   description = NULL,
		   tsrepresentation = function(data,dates){
	dt <- strptime(paste("01-",dates[1], sep=""), format="%d-%b-%Y")
	st <- c(1900+dt$year, dt$mon)
	ts(data,start=st, frequency=12)}))

  z <- TSget("DMACN", con2)

  require("zoo")
  
  con3 <- TSconnect("xls", dbname=testfile,
          map=list(ids  =list(i=11,     j="B:Q"), 
	           data =list(i=12:627, j="B:Q"), 
	           dates=list(i=12:627, j="A"),
                   names=list(i=4:7,    j="B:Q"), 
		   description = NULL,
		   tsrepresentation = function(data,dates){
		       zoo(data,order.by = as.Date(
			 paste("01-",dates, sep=""), format="%d-%b-%Y"))}))

  z <- TSget("DMACN", con3)

  TSrefperiod(z) 
  TSdescription(z) 

  unlink("Rplots.pdf")
