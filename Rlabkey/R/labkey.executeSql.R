##
#  Copyright (c) 2008-2015 Fred Hutchinson Cancer Research Center
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

labkey.executeSql <- function(baseUrl, folderPath, schemaName, sql, maxRows=NULL,
        rowOffset=NULL, colSort=NULL, showHidden=FALSE, colNameOpt='caption',
        containerFilter=NULL, parameters=NULL)
{

## Error if any of baseUrl, folderPath, schemaName or sql are missing
if(exists("baseUrl")==FALSE || exists("folderPath")==FALSE || exists("schemaName")==FALSE || exists("sql")==FALSE)
stop (paste("A value must be specified for each of baseUrl, folderPath, schemaName and sql."))

## URL encoding of schema and folder path (if not already encoded)
if(schemaName==curlUnescape(schemaName)) {schemaName <- curlEscape(schemaName)}
if(folderPath==URLdecode(folderPath)) {folderPath <- URLencode(folderPath)}
if(is.null(colSort)==FALSE) {if(colSort==curlUnescape(colSort)) colSort <- curlEscape(colSort)}
if(is.null(containerFilter)==FALSE) {if(containerFilter==curlUnescape(containerFilter)) containerFilter<- curlEscape(containerFilter)}

## Formatting
baseUrl <- gsub("[\\]", "/", baseUrl)
folderPath <- gsub("[\\]", "/", folderPath)
if(substr(baseUrl, nchar(baseUrl), nchar(baseUrl))!="/"){baseUrl <- paste(baseUrl,"/",sep="")}
if(substr(folderPath, nchar(folderPath), nchar(folderPath))!="/"){folderPath <- paste(folderPath,"/",sep="")}
if(substr(folderPath, 1, 1)!="/"){folderPath <- paste("/",folderPath,sep="")}

## Construct url
myurl <- paste(baseUrl,"query",folderPath,"executeSql.api?schemaName=",schemaName,"&apiVersion=8.3",sep="")
if(is.null(maxRows)==FALSE) {myurl <- paste(myurl,"&maxRows=",maxRows,sep="")}
if(is.null(maxRows)==TRUE) {myurl <- paste(myurl,"&showRows=all",sep="")}
if(is.null(rowOffset)==FALSE) {myurl <- paste(myurl,"&offset=",rowOffset,sep="")}
if(is.null(colSort)==FALSE) {myurl <- paste(myurl,"&query.sort=",colSort,sep="")}
if(is.null(parameters)==FALSE) {for(k in 1:length(parameters)) myurl <- paste(myurl,"&query.param.",parameters[k],sep="")}
if(is.null(containerFilter)==FALSE) {myurl <- paste(myurl,"&containerFilter=",containerFilter,sep="")}

## Set options
reader <- basicTextGatherer()
header <- basicTextGatherer()


handle <- getCurlHandle()
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
postForm(uri=myurl, "sql"=sql, .opts=myopts, curl=handle)



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
    }
    else {
        stop(paste("HTTP request was unsuccessful. Status code = ",status,", Error message = ",message,sep=""))
    }
}


newdata <- makeDF(rawdata=reader$value(), showHidden=showHidden, colNameOpt=colNameOpt)

return(newdata)
}

