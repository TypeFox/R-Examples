##' Get details of a ballot measure
##' 
##' This function is a wrapper for the Measure.getMeasure() method of the PVS API Measure class which grabs details of a ballot measure. The function sends a request with this method to the PVS API for all measure IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Measure.getMeasure(measureId)
##' @param measureId a character string or list of character strings with the measure ID(s) (see references for details)
##' @return A data frame with a row for each ballot measure and columns with the following variables describing the ballot measure:\cr measure.measureId,\cr measure.measureCode,\cr measure.title,\cr measure.electionDate,\cr measure.electionType,\cr measure.source,\cr measure.url,\cr measure.summary,\cr measure.summaryUrl,\cr measure.measureText,\cr measure.textUrl,\cr measure.proUrl,\cr measure.conUrl,\cr measure.yes,\cr measure.no,\cr measure.outcome.
##' @references http://api.votesmart.org/docs/Measure.html\cr
##' Use Measure.getMeasuresByYearState() to get a list of measure IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get details on certain ballot measures
##' \dontrun{measure <- Measure.getMeasure(list(1632,1633))}
##' \dontrun{measure}

##' @export






Measure.getMeasure <-
function (measureId) {

  
# internal function
Measure.getMeasure.basic <- function (.measureId) {
  
request <-  "Measure.getMeasure?"
inputs  <-  paste("&measureId=",.measureId, sep="")
output  <-  pvsRequest(request,inputs)
output$measureId <-.measureId
output

}


  # Main function
  output.list <- lapply(measureId, FUN= function (b) {
    
      Measure.getMeasure.basic(.measureId=b)
           
        
    }
  )

# which list entry has the most columns, how many are these?
coln <- which.is.max(sapply(output.list, ncol));
max.cols <- max(sapply(output.list, ncol));

# give all list entries (dfs in list) the same number of columns and the same names
output.list2 <- lapply(output.list, function(x){
if (ncol(x) < max.cols) x <- data.frame(cbind(matrix(NA, ncol=max.cols-ncol(x), nrow = 1, ),x),row.names=NULL)
names(x) <- names(output.list[[coln]])
x
})

output <- do.call("rbind",output.list2)
output


}
