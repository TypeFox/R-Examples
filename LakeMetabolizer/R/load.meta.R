#'@title Loads a metadata file from the specified path
#'@description 
#'Parses a formatted metadata file. Useful for site-specific metadata 
#'that is not contained in the timeseries files.
#'
#'@usage
#'load.meta(fPath)
#'@param fPath The file path as a string
#'
#'@return A list with the metadata parsed from the file.
#'
#'@keywords IO
#'@keywords file
#'
#'@author Luke A Winslow
#'
#'@seealso 
#'\link{load.ts}
#'\link{load.all.data}
#'@export
load.meta = function(fPath){
  
  fid = file(fPath, 'rt')
  tmp = readLines(fid, 100, warn=FALSE)
  close(fid)
  
  tmp = strsplit(tmp, '[\t,]') #splits on tabs *and* commas
  
  header = tolower(tmp[[1]])
  
  id.col = which(header=='id')
  val.col = which(header=='value')
  
  if(length(id.col) == 0 || length(val.col) == 0){
    stop('Metadata file must contain "ID" and "Value" header elements.')
  }
  
  output = list()
  
  for(i in 2:length(tmp)){
    id = tolower(tmp[[i]][id.col])
    val = as.numeric(tmp[[i]][val.col])
    
    if(is.na(val)){
      val = tmp[[i]][val.col]
    }
    output[id] = val
  }
  
  return(output)  
  
}