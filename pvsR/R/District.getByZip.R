##' Get district IDs according to the zip code
##' 
##' This function is a wrapper for the District.getByZip() method of the PVS API District class which grabs district IDs according to the ZIP code. The function sends a request with this method to the PVS API for all ZIP codes given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage District.getByZip(zip5, zip4=NULL)
##' @param zip5 a character string or list of character strings with the five-digit ZIP code
##' @param zip4 (optional) a character string or list of character strings with the expanded ZIP+4 code (default: All)
##' @return A data frame with a row for each district and columns with the following variables describing the district:\cr districtList.district*.districtId,\cr districtList.district*.name,\cr districtList.district*.officeId,\cr districtList.district*.stateId.
##' @references http://api.votesmart.org/docs/District.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get district IDs by ZIP Code
##' \dontrun{district <- District.getByZip(list(10001,10002),)}
##' \dontrun{district}
##' @export

District.getByZip <-
function (zip5, zip4=NULL) {
  
  if (length(zip4)==0) {
    # internal function
    District.getByZip.basic1 <- function (.zip5) {
      
      request <-  "District.getByZip?"
      inputs  <-  paste("&zip5=",.zip5,sep="")
      
      pvs.url <- paste("http://api.votesmart.org/",request,"key=",get('pvs.key',envir=.GlobalEnv),inputs,sep="") #generate url for request
      
      
      output.base <- xmlRoot(xmlTreeParse(pvs.url, useInternalNodes=TRUE))
      districts <-  removeChildren(output.base, kids=list(1,2))
      electionDistricts <- districts[["electionDistricts"]]
      districts <- removeChildren(districts, kids=list("electionDistricts"))
      
      output.districts <- data.frame(t(xmlSApply(districts, function(x) xmlSApply(x, xmlValue))), row.names=NULL)
      output.districts$Type <- "district"
      output.electionDistricts <- data.frame(t(xmlSApply(electionDistricts, function(x) xmlSApply(x, xmlValue))), row.names=NULL)
      output.electionDistricts$Type <- "electionDistrict"
      
      
      output <- rbind(output.districts, output.electionDistricts)      
      
      output$zip5 <- .zip5
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(zip5, FUN= function (s) {
      District.getByZip.basic1(.zip5=s)
    }
                          )
    
    
    output.list <- redlist(output.list)
    
    
    output <- dfList(output.list)
    
    
  } else {
    
    # internal function
    District.getByZip.basic2 <- function (.zip5, .zip4) {
      
      pvs.url <- paste("http://api.votesmart.org/",request,"key=",get('pvs.key',envir=.GlobalEnv),inputs,sep="") #generate url for request
      
      
      request <-  "District.getByZip?"
      inputs  <-  paste("&zip5=",.zip5, "&zip4=", .zip4, sep="")
      
      output.base <- xmlRoot(xmlTreeParse(pvs.url, useInternalNodes=TRUE))
      districts <-  removeChildren(output.base, kids=list(1,2))
      electionDistricts <- districts[["electionDistricts"]]
      districts <- removeChildren(districts, kids=list("electionDistricts"))
      
      output.districts <- data.frame(t(xmlSApply(districts, function(x) xmlSApply(x, xmlValue))), row.names=NULL)
      output.districts$Type <- "district"
      output.electionDistricts <- data.frame(t(xmlSApply(electionDistricts, function(x) xmlSApply(x, xmlValue))), row.names=NULL)
      output.electionDistricts$Type <- "electionDistrict"
      
      
      output <- rbind(output.districts, output.electionDistricts)
      
      output$zip5 <- .zip5
      output$zip4. <- .zip4
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(zip5, FUN= function (s) {
      lapply(zip4, FUN= function (c) {
        District.getByZip.basic2( .zip5=s, .zip4=c)
      }
             )
    }
                          )
    
    
    
    output.list <- redlist(output.list)
    
    output <- dfList(output.list)
    
    

    
    
    
  }
  
  
  
  
  
  output
  
  
  
}
