require("TSmisc")
require("tfplot")
require("zoo")

##################################################
##################################################

## test with only one series (which can cause a problem)
#### Data from Reserve Bank of Australia  ########

##################################################
##################################################


####  Australian Intervention copied Oct 28, 2013 ####
  # testfile <- "http://www.rba.gov.au/statistics/tables/xls/a05hist.xls"
  testfile <- system.file("xlsExampleData/a05hist.xls", package = "TSmisc")


# note that dates display in the spreadsheet like 02-Jan-1989, show in the
#   spreadsheet formula box like 1989-01-02, and import as "02/Jan/1989".
#   The last needs to be used.
#   If in doubt, look at con@dates
  con <- TSconnect("xls", dbname=testfile,
          map=list(ids  =list(i=11,      j="B"), 
	           data =list(i=12:6401, j="B"), 
	           dates=list(i=12:6401, j="A"),
                   names=list(i=4:5,     j="B"), 
		   description = NULL,
		   tsrepresentation = function(data,dates){
		     zoo(data,order.by = as.Date( dates, format="%d/%b/%Y"))}))

# The date format seems to change between the saved file an the URL
# For the URL use
#		     zoo(data,order.by = as.Date( dates, format="%d-%b-%Y"))}))


  z <- TSget("ARBANFXM", con)

  tfplot(z)
  tfplot(z, start=as.Date("2007-01-01"), end=as.Date("2009-06-01"))

 
  TSrefperiod(z) 
  TSdescription(z) 

  unlink("Rplots.pdf")
