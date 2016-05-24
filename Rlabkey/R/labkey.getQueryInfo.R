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

##labkey.getLookupDetails
labkey.getLookupDetails <- function(baseUrl, folderPath, schemaName, queryName, lookupKey)
{
if(exists("lookupKey")==FALSE )
{stop ("You must supply the key (name) value of a query field defined as a lookup type field.")}

lookupFields <- getQueryInfo(baseUrl=baseUrl, folderPath=folderPath, schemaName=schemaName, queryName=queryName,showDefaultView=FALSE, lookupKey=lookupKey)
return(lookupFields)
}

##Public getQueryDetails
labkey.getQueryDetails <- function(baseUrl, folderPath, schemaName, queryName)
{
queryDetails <- getQueryInfo(baseUrl=baseUrl, folderPath=folderPath, schemaName=schemaName, queryName=queryName,showDefaultView=FALSE)
return(queryDetails)
}

## Public getDefaultViewDetails
labkey.getDefaultViewDetails <- function(baseUrl, folderPath, schemaName, queryName)
{
viewDetails <- getQueryInfo(baseUrl=baseUrl, folderPath=folderPath, schemaName=schemaName, queryName=queryName,showDefaultView=TRUE)
return(viewDetails)
}

## internal reoutine that handles all of these
getQueryInfo <- function(baseUrl, folderPath, schemaName, queryName, showDefaultView=FALSE, lookupKey=NULL)
{
## Error if any of baseUrl, folderPath, schemName or queryName are missing
if(exists("baseUrl")==FALSE || exists("folderPath")==FALSE || exists("schemaName")==FALSE || exists("queryName")==FALSE )
{stop ("A value must be specified for each of baseUrl, folderPath, schemaName, and queryName.")}

if(is.null(lookupKey)==FALSE) {char <- nchar(lookupKey); if(char<1) {lookupKey<-NULL} }

## URL encoding (if not already encoded)
if(schemaName==curlUnescape(schemaName)) {schemaName <- curlEscape(schemaName)}
if(queryName==curlUnescape(queryName)) {queryName <- curlEscape(queryName)}
if(folderPath==URLdecode(folderPath)) {folderPath <- URLencode(folderPath)}
if(is.null(lookupKey)==FALSE) {if(lookupKey==curlUnescape(lookupKey)) lookupKey <- curlEscape(lookupKey)}

## Formatting
baseUrl <- gsub("[\\]", "/", baseUrl)
folderPath <- gsub("[\\]", "/", folderPath)
if(substr(baseUrl, nchar(baseUrl), nchar(baseUrl))!="/"){baseUrl <- paste(baseUrl,"/",sep="" )}
if(substr(folderPath, nchar(folderPath), nchar(folderPath))!="/"){folderPath <- paste(folderPath,"/",sep="")}
if(substr(folderPath, 1, 1)!="/"){folderPath <- paste("/",folderPath,sep="")}

## Construct url
myurl <- paste(baseUrl,"query",folderPath,"getQueryDetails.api?schemaName=", schemaName, "&queryName=", queryName, "&apiVersion=8.3", sep="")
if(is.null(lookupKey)==FALSE) {myurl <- paste(myurl,"&fk=",lookupKey,sep="")}

## Set options
reader <- basicTextGatherer()
header <- basicTextGatherer()
myopts <- curlOptions(netrc=1, writefunction=reader$update, headerfunction=header$update, .opts=c(labkey.curlOptions()))

## Support user-settable options for debuggin and setting proxies etc
if(exists(".lksession"))
{
	userOpt <- .lksession[["curlOptions"]] 
	if (!is.null(userOpt))
		{myopts<- curlOptions(.opts=c(myopts, userOpt))}
}

## Http get
handle <- getCurlHandle()
clist <- ifcookie()
if(clist$Cvalue==1) {mydata <- getURI(myurl, .opts=myopts, cookie=paste(clist$Cname,"=",clist$Ccont,sep=""))} else {mydata <- getURI(myurl, .opts=myopts, curl=handle)}

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

## If querying the default view, the metadata is in a differnt object in the json stream 
if (showDefaultView==TRUE) {qcs<-decode$defaultView$columns}
else {qcs <- decode$columns}

## parsed JSON stream has two types of problems related to nulls:  
## the value NULL as as named element of the parent node
## the absence of either a value or a name for some columns on some records
## etiher one can be detected by checking for is.null on a row-by row basis against the total set of column names

baseColumns <- c("name", "caption", "fieldKey", "type", "isNullable","isKeyField",
			"isAutoIncrement", "isVersionField","isHidden","isSelectable",
			"isUserEditable", "isReadOnly", "isMvEnabled","description")
lookupColumns <- c("keyColumn", "schemaName", "displayColumn", "queryName", "isPublic")

dmall <- matrix(nrow=0, ncol=20, byrow=TRUE)
if(length(qcs)>0)
{
	for (j in 1:length(qcs))
	{
		dmqrow<- matrix(data=decode$name[[1]], nrow=1, ncol=1, byrow=FALSE)
		for (nm in baseColumns) {
			if (is.null(qcs[[j]][[nm]])) {qcs[[j]][[nm]] <- NA}
			
			dmqrow<- matrix(data=cbind(dmqrow, qcs[[j]][[nm]]), nrow=1, byrow=FALSE)
		}			

		
		if (is.null(qcs[[j]]$lookup))
		{
			lookupinfo <- matrix(data=cbind(NA,NA,NA,NA,NA), ncol=5, byrow=FALSE)
		}			
		else
		{
			for (nm in lookupColumns) {
				if (is.null(qcs[[j]]$lookup[[nm]])) {qcs[[j]]$lookup[[nm]] <- NA}
			}			
			nm <- lookupColumns[1]
			lookupinfo<- as.matrix(qcs[[j]]$lookup[[nm]], nrow=1, byrow=FALSE)
			for (nm in lookupColumns[-1]) {
				lookupinfo<- matrix(data=cbind(lookupinfo, qcs[[j]]$lookup[[nm]]), nrow=1, byrow=FALSE)
			}			
		}		
		
		dmqrow<-cbind(dmqrow, lookupinfo)
		dmall <- rbind(dmall,dmqrow)		
	}
}
dfall <- as.data.frame(dmall, stringsAsFactors=FALSE)
colnames(dfall)<-c("queryName", "fieldName", "caption", "fieldKey", "type", "isNullable","isKeyField",
			"isAutoIncrement", "isVersionField","isHidden","isSelectable",
			"isUserEditable", "isReadOnly", "isMvEnabled", "description",
			"lookupKeyField","lookupSchemaName","lookupDisplayField", "lookupQueryName", "lookupIsPublic")

return(dfall)
}
