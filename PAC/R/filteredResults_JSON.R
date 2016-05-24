#' @import RJSONIO
NULL

#' Obtain the input to zoomable packed circles plot. The clusters are filtered by size.
#' 
#' @param toJSON_Input Variable name of the JSON input file.
#' @param CountFilter The minimum size of retained clusters in each sample. Set to 0 for no filtering.
#' @param filename The filename for zoomable packed circles plot. 
#' 
#' @return The JSON format input to D3 zoomable packed circles plot.
#' @export

filteredResults_JSON<-function(toJSON_Input, CountFilter = 1000, filename = "zoomablePackedCirclesInput.txt"){
  
  output1<-as.data.frame(toJSON_Input)
  output1[,3]<-as.numeric(as.character(output1[,3]))
  output1<-subset(output1, output1[,ncol(output1)] >= CountFilter)
  
  intermediate<-list(name="jsonOut",children=makeList(output1))
  jsonOut<-gsub("\\\n", "", toJSON(as.list(intermediate))) #jsonOut 
  
  write(jsonOut, file=filename,sep="")
}

