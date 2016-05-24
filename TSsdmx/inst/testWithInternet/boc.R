
#################### Bank of Canada #################### 

require("TSsdmx")
require("tfplot")

if (FALSE) {

# This is not REST but URL calls return SDMX.

boc <- TSconnect("sdmx", dbname="BoC")
z <- TSget(c('CDOR', 'OIS'), boc)

# Finding series identifiers is difficult and
#    limited data is available (as of December 2010).

con <- TSconnect("sdmx", dbname="BoC") 

# identifiers can be extraced at ??

# z <- TSget("CDOR", con=con)
# above seems to end with funny <Series/>
# </Obs></Series><Series/></DataSet></CompactData></return>



z <- TSge(c("CDOR", "OIS", "SWAPPEDTOFLOAT"), con=con)

z <- TSgetBoC(c("CDOR", "OIS", "SWAPPEDTOFLOAT"))

tfplot(z, Title="From Bank of Canada")
TSdescription(z) 

#  not sure if these series exist, or some other problem
# options(TSconnection=con)
# 
# x <- TSget(c("TOTALSL","TOTALNS"), con, 
#        names=c("Total Consumer Credit Outstanding SA",
#                "Total Consumer Credit Outstanding NSA"))
# plot(x)
# tfplot(x)
# TSdescription(x) 

#http://credit.bank-banque-canada.ca/webservices?service=getSeriesSDMX&args=CDOR_-_-_FIRST_-_-_Last

TSgetBoC <- function(id, names=NULL){
   uri <- paste( "http://credit.bank-banque-canada.ca/webservices?service=getSeriesSDMX&args=",
   	    paste(id, "_-_-_", sep="", collapse=""),
   	    paste( "FIRST_-_-_Last",sep="", collapse=""), sep="")

   # should try to check <faultstring> 
   z <- getURLContent(uri)
   #NEED TO STRIP SOAP ENV
   #zz <- htmlParse(z)
   #zz <-   getNodeSet(htmlParse(z), "//series" )
   #zz <-   getNodeSet(zz, "//dataset" )
   #zz <-   getNodeSet(zz, "//CompactData" )
   #zz <-   getNodeSet(htmlParse(z, useInternalNodes = FALSE), "//return" )
   #zz <-   getNodeSet(htmlParse(z), "//return" )
   
   # should be able to parse for nmsp as in TSgetECB
   nmsp <-  c(ns="http://www.SDMX.org/resources/SDMXML/schemas/v2_0/message") 
   #        getNodeSet(htmlParse(z), "//series", nmsp )
   # length(getNodeSet(htmlParse(z), "//series", nmsp ) )
   # mode(getNodeSet(htmlParse(z), "//series", nmsp ) )

   #doc <- xmlParse(z)
   doc <- htmlParse(z)
   
   # DataSetParse assumes Xpath points to  Series nodes in doc 
   #   so getNodeSet(doc, Xpath, nmsp ) is a list with series as elements
   #   so getNodeSet(doc, "//series[@freq]", nmsp ) is a list with series as elements
   #r <- DataSetParse(doc,"//ns:Series[@FREQ]" ,nmsp)
   #r <- DataSetParse(doc,"//series[@freq]" ,nmsp)
   r <- DataSetParse(doc,"//series[@freq]" ,nmsp,
    obs=".//obs[@time_period]", timeperiod="time_period", value="obs_value")
   
   if(nseries(r) != length(id)) warning("some series not retrieved.")

   if(!is.null(names)) seriesNames(r) <- names

   r
   }

}
