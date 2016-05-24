#' Basic Price Information on natural Gas in Germany
#'
#' This function exposes an interface to the Netconnect Germany website
#' and allows the download of price data on natural gas.
#'
#' @param clss Class name for the object, default UrlData2
#'
#' @return UrlData2
#' @references 
#' \href{http://datenservice.net-connect-germany.de/Dokumente/NCG_XML_Interface_V1.6_de.pdf}{NCG}
#' @export
netConnectGermany <- function(clss="UrlData") urldata(
  clss=clss,
  resource="BasicPriceData",
  template="http://datenservice.net-connect-germany.de/XmlInterface/getXML.ashx?ReportId=BasicPriceData&Start=$(from)&End=$(to)",
  from=function(x=NULL) {if(is.null(x)) x <- Sys.time() - as.difftime(7, units="days"); strftime(as.Date(x), "%d-%m-%Y")},
  to=function(x=NULL) {if(is.null(x)) x <- Sys.time(); strftime(x, "%d-%m-%Y")},

  extract.fct=XML::xmlParse,
  transform.fct=function(dat) {
    namespaces <- c(
      sqlrowset="urn:schemas-microsoft-com:sql:SqlRowSet1",
      xsd="http://www.w3.org/2001/XMLSchema-instance",
      sqltypes="http://schemas.microsoft.com/sqlserver/2004/sqltypes"
      
    )
    nodes <- XML::getNodeSet(dat, "//sqlrowset:Price", namespaces=namespaces)
    dat <- sapply(nodes, function(n) c(
        XML::xmlValue(n[["Gasday"]]), 
        XML::xmlValue(n[["NCG_one_Day_Ahead_Settl_Price"]]), 
        XML::xmlValue(n[["GASPOOL_one_Day_Ahead_Settl_Price"]]),
        XML::xmlValue(n[["TTF"]]),
        XML::xmlValue(n[["Zeebruegge"]]),
        XML::xmlValue(n[["Unit"]])
    )) 
    dat <- data.frame(
        Gasday=strptime(dat[1,], "%Y-%m-%dT%H:%M:%S"), # 2010-10-05T06:00:00
        NCG_one_Day_Ahead_Settl_Price=as.numeric(dat[2,]), 
        GASPOOL_one_Day_Ahead_Settl_Price=as.numeric(dat[3,]), 
        TTF=as.numeric(dat[4,]), 
        Zeebruegge=as.numeric(dat[5,]), 
        Unit=dat[6,], 
        stringsAsFactors=FALSE
    )
    return(dat)
  }
)
