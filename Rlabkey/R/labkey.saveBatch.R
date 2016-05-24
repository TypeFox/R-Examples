##
#  Copyright (c) 2010-2014 LabKey Corporation
# 
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
##

labkey.saveBatch <- function(baseUrl, folderPath, assayName, resultDataFrame, batchPropertyList=NULL, runPropertyList=NULL)
{  

## Error if any of baseUrl, folderPathare missing
if(exists("baseUrl")==FALSE || exists("folderPath")==FALSE || exists("assayName")==FALSE  || exists("resultDataFrame")==FALSE)
stop (paste("A value must be specified for each of baseUrl, folderPath, assayName, and resultDataFrame"))
## TODO.check for at least one of the instert blocks is not null


## Formatting
baseUrl <- gsub("[\\]", "/", baseUrl)
folderPath <- gsub("[\\]", "/", folderPath)
if(substr(baseUrl, nchar(baseUrl), nchar(baseUrl))!="/"){baseUrl <- paste(baseUrl,"/",sep="")}
if(substr(folderPath, nchar(folderPath), nchar(folderPath))!="/"){folderPath <- paste(folderPath,"/",sep="")}
if(substr(folderPath, 1, 1)!="/"){folderPath <- paste("/",folderPath,sep="")}

## URL encode folder path and assay name (if not already encoded) 
if(folderPath==URLdecode(folderPath)) {folderPath <- URLencode(folderPath)}
if(assayName==curlUnescape(assayName)) {assayNameParam <- curlEscape(assayName)}
else {assayNameParam <- assayName}

## Set options
reader <- basicTextGatherer()
header <- basicTextGatherer()
handle <- getCurlHandle()
headerFields <- c('Content-Type'="application/json;charset=utf-8")
clist <- ifcookie()
if(clist$Cvalue==1) {
    myopts<- curlOptions(cookie=paste(clist$Cname,"=",clist$Ccont, sep=""), writefunction=reader$update, headerfunction=header$update,
                        .opts=c(labkey.curlOptions()))
} else {
    myopts<- curlOptions(netrc=1, writefunction=reader$update, headerfunction=header$update,
                        .opts=c(labkey.curlOptions()))
}

## Support user-settable options for debugging and setting proxies etc
if(exists(".lksession"))
{
	userOpt <- .lksession[["curlOptions"]] 
	if (!is.null(userOpt))
		{myopts<- curlOptions(.opts=c(myopts, userOpt))}
}

## Translate assay name to an ID
myurl <- paste(baseUrl,"assay",folderPath,"assayList.view?name=", assayNameParam, sep="")

assayInfoJSON <- getURI(myurl, .opts=myopts, curl=handle)
h <- parseHeader(header$value())
status <- getCurlInfo(handle)$response.code
message <- h$statusMessage
if(status>=400)
  {stop(paste("Could not find assay by that name. Status code = ",status,", Error message = ",message,sep=""))}

assayDef <- NULL
assayInfo<- fromJSON(assayInfoJSON)
if (length(assayInfo) == 1 && length(assayInfo[[1]]) == 1)
{
	assayDef <- assayInfo[[1]][[1]]
	if (assayDef$name != assayName)
		{assayDef <- NULL}
	## TODO:  check assay domain def against dataframe	
}
if (is.null(assayDef))
	{stop(paste("Could not find an assay matching that name." ,sep=""))}
	
# build Assay object tree based on R lists
nrows <- nrow(resultDataFrame)
ncols <- ncol(resultDataFrame)
cnames <- colnames(resultDataFrame)
rowsVector <- vector(mode="list", length=nrows)
for(j in 1:nrows) {
	cvalues <- as.list(resultDataFrame[j,])
	names(cvalues) <- cnames
	rowsVector[[j]] <- cvalues
}

dataInputsArray <- vector(mode="list", length=0)

runsArray <- vector(mode="list", length=1)
runPropertyList <- c(runPropertyList, list("dataInputs"=dataInputsArray))

runsArray[[1]] <- c(runPropertyList, list("dataRows" = rowsVector))

batchPropertyList<- c(batchPropertyList, list("runs"= runsArray))
	
baseAssayList <- list(assayId=assayDef$id)
baseAssayList <- c(baseAssayList, list(batch=batchPropertyList))

## Now post form with batch object filled out
myurl <- paste(baseUrl,"assay",folderPath,"saveAssayBatch.view", sep="")
pbody<-toJSON(baseAssayList)
curlPerform(url=myurl, postFields=pbody, httpheader=headerFields, .opts=myopts, curl=handle)


## Error checking for incoming file
h <- parseHeader(header$value())
status <- getCurlInfo(handle)$response.code
message <- h$statusMessage
if(status==500) 
{decode <- fromJSON(reader$value()); message <- decode$exception; stop(paste("HTTP request was unsuccessful. Status code = ",status,", Error message = ",message,sep=""))}
if(status>=400)
{
    contTypes <- which(names(h)=='Content-Type')
	if(length(contTypes)>0 && (tolower(h[contTypes[1]])=="application/json;charset=utf-8" || tolower(h[contTypes[2]])=="application/json;charset=utf-8"))
    {
        decode <- fromJSON(reader$value());
        message<-decode$exception;
        stop (paste("HTTP request was unsuccessful. Status code = ",status,", Error message = ",message,sep=""))
    } else
	{
	    stop(paste("HTTP request was unsuccessful. Status code = ",status,", Error message = ",message,sep=""))
    }
}

newAssayInfo <- fromJSON(reader$value())

return(newAssayInfo)
}
                                                              
