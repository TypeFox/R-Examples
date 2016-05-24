##
# Copyright (c) 2010-2015 LabKey Corporation
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##

labkey.getFolders <- function(baseUrl, folderPath, includeEffectivePermissions=TRUE, includeSubfolders=FALSE, depth=50)
{
## Empty string/NULL checking

## Error if either baseUrl or folderPath
if(exists("baseUrl")==FALSE || exists("folderPath")==FALSE)
stop (paste("A value must be specified for both baseUrl and folderPath"))

## URL encoding of folder path  (if not already encoded)
if(folderPath==URLdecode(folderPath)) {folderPath <- URLencode(folderPath)}

## Formatting
baseUrl <- gsub("[\\]", "/", baseUrl)
folderPath <- gsub("[\\]", "/", folderPath)
if(substr(baseUrl, nchar(baseUrl), nchar(baseUrl))!="/"){baseUrl <- paste(baseUrl,"/",sep="")}
if(substr(folderPath, nchar(folderPath), nchar(folderPath))!="/"){folderPath <- paste(folderPath,"/",sep="")}
if(substr(folderPath, 1, 1)!="/"){folderPath <- paste("/",folderPath,sep="")}
if(includeSubfolders) {inclsf <- paste("1&depth=", depth, sep="")} else {inclsf <- "0"}
if(includeEffectivePermissions) {inclep <- "1"} else {inclep <- "0"}

## Construct url
myurl <- paste(baseUrl,"project",folderPath,"getContainers.view?","includeSubfolders=",inclsf,"&includeEffectivePermissions=",inclep, sep="")

## Set options
reader <- basicTextGatherer()
header <- basicTextGatherer()
myopts <- curlOptions(netrc=1, writefunction=reader$update, headerfunction=header$update, .opts=c(labkey.curlOptions()))

## Support user-settable options for debugging and setting proxies etc
if(exists(".lksession"))
{
	userOpt <- .lksession[["curlOptions"]] 
	if (!is.null(userOpt))
		{myopts<- curlOptions(.opts=c(myopts, userOpt))}
}

## Http get
handle <- getCurlHandle()
clist <- ifcookie()
if(clist$Cvalue==1) 
{
	mydata <- getURI(myurl, .opts=myopts, cookie=paste(clist$Cname,"=",clist$Ccont,sep=""))
}
else 
{
	myopts <- curlOptions(.opts=c(myopts, httpauth=1L))
	mydata <- getURI(myurl, .opts=myopts, curl=handle)
}

## Error checking, decode data and return data frame
h <- parseHeader(header$value())
status <- getCurlInfo(handle)$response.code
message <- h$statusMessage

if(status==500)
{decode <- fromJSON(mydata); message <- decode$exception; stop(paste("HTTP request was unsuccessful. Status code = ",status,", Error message = ",message,sep=""))}
if(status>=400)
{
  contTypes <- which(names(h)=='Content-Type')
  if(length(contTypes)>0 && (tolower(h[contTypes[1]])=="application/json;charset=utf-8" || tolower(h[contTypes[2]])=="application/json;charset=utf-8"))
  {
    decode <- fromJSON(mydata);
    message<-decode$exception;
    stop (paste("HTTP request was unsuccessful. Status code = ",status,", Error message = ",message,sep=""))
  } else
  {
    stop(paste("HTTP request was unsuccessful. Status code = ",status,", Error message = ",message,sep=""))
  }
}

decode <- fromJSON(mydata)
curfld <- decode
allpaths <- matrix(data=c(curfld$name, curfld$path, paste(curfld$effectivePermissions, collapse=",")), nrow=1, ncol=3, byrow=TRUE)
todo <- curfld$children[]
while (length(todo)>0)
{
	curfld<-todo[1][[1]]
	allpaths <- rbind(allpaths, c(curfld$name, curfld$path, paste(curfld$effectivePermissions, collapse=",")))
	todo<- c(todo, curfld$children[])
	todo<-todo[-1]
}

allpathsDF <- data.frame(allpaths, stringsAsFactors=FALSE)
colnames(allpathsDF) <- c("name", "folderPath", "effectivePermissions")

return(allpathsDF)

}


