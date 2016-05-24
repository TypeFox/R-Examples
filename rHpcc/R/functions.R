# This function is used to trim leading or trailing whitespaces
trim <- function (x) {
  gsub("^\\s+|\\s+$", "", x)
}

# Executes the ECL code on the cluster specified and returns the XML response
eclDirectCall <- function(hostName, port="8008", clusterName="thor", eclCode) {
  if(is.null(hostName) || !nzchar(trim(hostName))) {
    stop("HostName required")
  } else {
    urlString <- paste("http://",hostName,":",port,"/EclDirect/RunEcl",sep="");
    content <- postForm(urlString, eclText=enc2utf8(eclCode), cluster=enc2utf8(clusterName))
    content
  }
}

# Parses the XML returned from eclDirectCall() and allows you to download the result in either CSV or XML format
parseResults <- function(xmlResult, downloadPath, format){
  data
  if(missing(xmlResult)) {
    stop("Empty XML String.")
  } else {
    docRoot = xmlRoot(xmlTreeParse(xmlResult))
    nodes = getNodeSet(docRoot, "//Dataset")
    for(i in 1:length(nodes)) {
      datasetNode <- nodes[[i]]
      resultSetName = xmlGetAttr(datasetNode, "name")
      x <- array(1:length(datasetNode)*length(datasetNode[[1]]), dim=c(length(datasetNode),length(datasetNode[[1]])))
      for(j in 1:length(datasetNode)) {
        rowNode <- datasetNode[[j]]
        for(k in 1:length(rowNode)) {
          actualNode <- rowNode[[k]]
          x[j,k] <- xmlValue(actualNode)
        }
      }
      
      # Download the file if "downloadPath" argument is specified
      if(!missing(downloadPath)) {
        if(missing(format)) {
          fileName <- paste(resultSetName,".csv", sep = "")
          path <- paste(downloadPath,"/",fileName, sep = "")
          write.table(x,file=path,sep=",",row.names=F, col.names = F)
        } else {
          fileName <- paste(resultSetName,".xml", sep = "")
          path <- paste(downloadPath,"/",fileName, sep = "")
          write.table(xmlResult,file=path,sep=",",row.names=F, col.names = F)
        }   
      }      
      data <- x
      data = as.data.frame(data, stringsAsFactors=FALSE)
    }
    
    data
  } 
}

# Takes fixed-format file from the landing zone and distributes it across the nodes of the destination supercomputer.
# query <- sprayFixed(ip="127.0.0.1", filePath="/var/lib/HPCCSystems/mydropzone/sampleFile.txt", 
#             recordlength=255, clusterName="mythor", logicalFileName="IN::MyFile")
# eclDirectCall(hostName = "127.0.0.1", eclCode=query)
sprayFixed <- function(ip , filePath, recordlength, clusterName="mythor", logicalFileName) {
  if(missing(ip) || is.null(ip) || !nzchar(trim(ip))) {
    stop("ip required in function: sprayField()")
  } else if(missing(filePath) || is.null(filePath) || !nzchar(trim(filePath))) {
    stop("filePath required in function: sprayField()")
  } else if(missing(recordlength) || is.null(recordlength) || !nzchar(trim(recordlength))) {
    stop("recordlength required in function: sprayField()")
  } else if(missing(logicalFileName) || is.null(logicalFileName) || !nzchar(trim(logicalFileName))) {
    stop("logicalFileName required in function: sprayField()")
  } else {
    url <- paste("http://",ip,":8010/FileSpray",sep="")
    result <- paste("IMPORT STD; ","STD.File.SprayFixed('",ip,"','",filePath,"',",recordlength,",'",clusterName,"','",logicalFileName,"',-1,'",url,"');",sep="")
    
    result
  } 
}

# Takes a variable length file from the landing zone and distributes it across the nodes of the destination supercomputer
# query <- sprayVariable(ip="127.0.0.1", filePath="/var/lib/HPCCSystems/mydropzone/sampleFile.csv", clusterName="mythor", logicalFileName="IN::MyFile")
# eclDirectCall(hostName = "127.0.0.1", eclCode=query)
sprayVariable <- function(ip, filePath, clusterName="mythor", logicalFileName) {
  if(missing(ip) || is.null(ip) || !nzchar(trim(ip))) {
    stop("ip required in function: sprayField()")
  } else if(missing(filePath) || is.null(filePath) || !nzchar(trim(filePath))) {
    stop("filePath required in function: sprayField()")
  } else if(missing(logicalFileName) || is.null(logicalFileName) || !nzchar(trim(logicalFileName))) {
    stop("logicalFileName required in function: sprayField()")
  } else {
    url <- paste("http://",ip,":8010/FileSpray",sep="")
    result <- paste("IMPORT STD; STD.File.SprayVariable('",ip,"','",filePath,"',,,,,'",clusterName,"','",logicalFileName,"',-1,'",url,"');",sep="")
    result
  }
}
