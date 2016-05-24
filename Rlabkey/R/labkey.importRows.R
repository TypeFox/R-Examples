##
#  Copyright (c) 2013-2015 LabKey Corporation
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

labkey.importRows <- function(baseUrl, folderPath, schemaName, queryName, toImport)
{
## Default showAllRows=TRUE
showAllRows=TRUE

## Error if any of baseUrl, folderPath, schemName or toImport are missing
if(exists("baseUrl")==FALSE || exists("folderPath")==FALSE || exists("schemaName")==FALSE || exists("toImport")==FALSE)
stop (paste("A value must be specified for each of baseUrl, folderPath, schemaName and toImport."))

## Formatting
baseUrl <- gsub("[\\]", "/", baseUrl)
folderPath <- gsub("[\\]", "/", folderPath)
if(substr(baseUrl, nchar(baseUrl), nchar(baseUrl))!="/"){baseUrl <- paste(baseUrl,"/",sep="")}
if(substr(folderPath, nchar(folderPath), nchar(folderPath))!="/"){folderPath <- paste(folderPath,"/",sep="")}
if(substr(folderPath, 1, 1)!="/"){folderPath <- paste("/",folderPath,sep="")}

## URL encode folder path, JSON encode post body (if not already encoded)
toImport <- convertFactorsToStrings(toImport);
if(folderPath==URLdecode(folderPath)) {folderPath <- URLencode(folderPath)}
nrows <- nrow(toImport)
ncols <- ncol(toImport)
p1 <- toJSON(list(schemaName=schemaName, queryName=queryName, apiVersion=8.3))
cnames <- colnames(toImport)
p3 <- NULL
for(j in 1:nrows)
{
    cvalues <- as.list(toImport[j,])
	names(cvalues) <- cnames
	cvalues[is.na(cvalues)] = NULL
    p2 <- toJSON(cvalues)
    p3 <- c(p3, p2)
}
p3 <- paste(p3, collapse=",")
pbody <- paste(substr(p1,1,nchar(p1)-1),', \"rows\":[',p3,"] }",sep="")


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

## Support user-settable options for debuggin and setting proxies etc
if(exists(".lksession"))
{
	userOpt <- .lksession[["curlOptions"]]
	if (!is.null(userOpt))
		{myopts<- curlOptions(.opts=c(myopts, userOpt))}
}

## Post form
myurl <- paste(baseUrl,"query",folderPath,"importRows.api",sep="")
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
    } else{
        stop(paste("HTTP request was unsuccessful. Status code = ",status,", Error message = ",message,sep=""))
    }
}

newdata <- fromJSON(reader$value())

return(newdata)
}

