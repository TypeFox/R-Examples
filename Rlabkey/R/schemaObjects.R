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

	## an environment to hold metadata descriptions
	.lksite <- new.env(parent=emptyenv())
	.lksession <- new.env(parent=emptyenv())

	##########################################
	##  public functions
	##
	## a session holds current selected values of the site schema root and the user's currnt folder / container
	## These are put in two different environment spaces, the schema root is keyed by the URL because the schema is not going to be 
	## different if the URL is the same.  
	## The session holds pointers to the site and session variables and holds the baseUrl and current folder path
	##
	## Two option buckets are created at the session root so that they  can be read by network-level code without having to pass
	## them down through all the calls.

	getSession <-
	    function(baseUrl, folderPath="/home", curlOptions=NULL, lkOptions=NULL)
	{																					
		skey <-  gsub("[ :]*", "", as.character(date()))
		sitekey <- paste(baseUrl, "/", .getProjectFromPath(folderPath))
		.lksite[[sitekey]]<- list(NA)

		lks<- .LabkeySession(baseUrl=baseUrl, folderPath=folderPath, lksite=.lksite, skey=skey)

		.lksession[[skey]]<- list("baseUrl"=baseUrl, "folderPath"=folderPath, "validSchemas"=NA, "lkOptions" = NA)

		if (missing(lkOptions)) {lkOptions <- list(NA)}	
		.lksession$skey$lkOptions <- lkOptions

		# put curlOptions on a package-wide environment so we don't need to pass them everywhere
		.lksession[["curlOptions"]] <- curlOptions 

		.setupSchemas(lks, folderPath)	

		lks
	}
	
	## on extracting a schema (a named member of the session), fill in the schema cache
	##
##	`[[.LabkeySession`<- function(x, s, ...)
##	{
##		sc <- getSchema(x, s)
##		return (sc)	
##	}


	##   returns the set of fields accessible through a lookup col.  these are not cached
	##
	getLookups<-
		function(session,lookupField)
	{
		if(!inherits(session, "LabkeySession") )
			{stop("getLookups() requires a session object")}
						
		if (!inherits(lookupField, "LabkeyField")  || is.null(lookupField$lookupQueryName))
			{stop("getLookups() requires a LabkeyField object with a non-empty lookupQueryName property ")}	
					
		schemaName <- attr(lookupField, "schemaName")
		queryName <- attr(lookupField, "queryName")
		
		
		out <- labkey.getLookupDetails(attr(session,"baseUrl"), getFolderPath(session), schemaName, queryName, lookupField$fieldName )	
		fieldList <- list()
		for (r in row.names(out))
		{				
			rowList <- as.list(out[r,])
			fieldList[[as.character(rowList$fieldName)]]<- .LabkeyField(rowList, schemaName=schemaName, queryName=queryName)
		}

		return (.LabkeyQuery(x=fieldList, schemaName=schemaName, queryName=queryName))
	}

	#########################
	###  wrapper on labkey.selectRows that makes all these objects worthwhile
	getRows<-
		function(session, query, maxRows=NULL, colNameOpt='fieldname', ...)
	{
		if(!inherits(session, "LabkeySession") && inherits(query, "LabkeyQuery"))
			{stop("getRows() requires a session and a schema$query object")}
		if (length(attr(query, "schemaName"))==0 || length(attr(query, "queryName"))==0)
			{stop("Invalid query object")}

		## support session defaults
		skey <- attr(session, "skey")
		dflt <- .lksession$skey$lkOptions

		if (length(dflt)>0 && !is.na(dflt))
		{
			for (nm in names(dflt) )
			{
				if ((nm=="colNameOpt") && missing(colNameOpt)) {colNameOpt <- dflt$colNameOpt}
				if ((nm=="maxRows") && missing(maxRows)) {maxRows <- dflt$maxRows}
			}
		}
		lkdata <- labkey.selectRows(baseUrl=attr(session,"baseUrl"), folderPath=getFolderPath(session), schemaName=attr(query, "schemaName")
				, queryName=attr(query, "queryName"), maxRows=maxRows, colNameOpt=colNameOpt, ...)	
		return(lkdata)			

	}


	############################################
	##  list available schemas (given base Url and current folder path) 
	##
	lsSchemas <-
		function(session)
	{
		return(.getSchemasAsList(session))
	}


	############################################
	##  list available folders (given base Url and current folder path) 
	##
	lsFolders <-
		function(session)
	{
		folderPath <- getFolderPath(session)
		projectName <- .getProjectFromPath(folderPath)
		folders <- labkey.getFolders(attr(session,"baseUrl"), projectName, includeSubfolders=TRUE)
		return (sort(as.array(folders$folderPath)))	

	}



	############################################
	##  list available projects at base Url 
	##
	lsProjects <-
		function(baseUrl)
	{
		folders <- labkey.getFolders(baseUrl, "/", includeSubfolders=TRUE, depth=1)
		folders <- folders[(folders$folderPath != "/"),]
		return (sort(as.array(folders$folderPath)))	

	}

	##  getter for the folderPath 
	getFolderPath <-
	    function(session)
	{
		return(attr(session, "folderPath"))
	}

	############################################
	##  save a dataframe as an assay result
	##
	saveResults <-
		function(session, assayName, resultDataFrame, 
					batchPropertyList= list(name=paste("Batch ", as.character(date()))), 
					runPropertyList= list(name=paste("Assay Run ", as.character(date()))) )

	{
		assayInfo <- labkey.saveBatch(
				baseUrl=attr (session, "baseUrl"), 
				folderPath=getFolderPath(session), 
				assayName=assayName, 
				resultDataFrame=resultDataFrame, 
				batchPropertyList=batchPropertyList, 
				runPropertyList=runPropertyList)
				
		return(assayInfo)		
	}



	###########################################################################################
	##############################################
	##  Private functions

	##############################################
	##
	## .getProjectFromPath
	.getProjectFromPath <- function (folderPath)
	{
		## find the project name part
		path <- gsub("[\\]", "/", folderPath)
		pathParts <- strsplit(path,"/")[[1]]
		pathParts <- pathParts[pathParts!=""]
		if (length(pathParts)==0) 
			{projectName<- "/"} 
		else 
			{projectName<-pathParts[1]}
		return(projectName)	
	}

	##############################################
	##
	## get schema by name or index.  

	getSchema <-
	    function(session,schemaIndex)
	{	
		slist <- .getSchemasAsList(session)
		if (is.character(schemaIndex))
		{
			if (is.null(slist[[schemaIndex]]) )
				{ stop("Can't find schema by that name ") }
		}
		else
		{
			if(schemaIndex > length(slist))
				{stop("Can't find a schema by that number")}
		}

		sname <- slist[[schemaIndex]]

		## if new queriesList comes back empty, don't use it-- the cache will remain as it was
		newQueriesList <- .getQueriesList(session, schemaName=sname)
		
		lksite<- attr(session, "lksite")
		if (length(newQueriesList)>0) 
		{				
			lksite[[sname]]<- .LabkeySchema(x=newQueriesList, sname=sname)
		}

		return(lksite[[sname]])

	}


	##############################################
	##
	## load schemas from labkey server.  don't load queries until a specific schema is requested

	.setupSchemas <-
	    function(session, folderPath)
	{
		baseUrl<-attr(session,"baseUrl")
		out <- labkey.getSchemas(baseUrl, folderPath)
		## build a list of valid scheas for this folder context
		validSchemasDF <- cbind (out, date())
		validSchemasDF <- data.frame(validSchemasDF[order(validSchemasDF[,1]), ], stringsAsFactors=FALSE)
		colnames(validSchemasDF) <- c("schemaName", "timestamp")


		skey=attr(session, "skey") 
		.lksession$skey[["validSchemas"]] <- validSchemasDF	

		lksite<- attr(session, "lksite")
		for (sn in out$schemaName)
		{		
			if (length(lksite[[as.character(sn)]])==0) {
				lksite[[as.character(sn)]]<- structure(as.list(NA))
			}	
		}
	##    invisible(TRUE)                     # quiet success
	}


	.getQueriesList <- function(session, schemaName)
	{		
		out <- labkey.getQueries(attr(session,"baseUrl"), getFolderPath(session), schemaName)
		## if the existing cached schema passes the checks in .checkValid, then we return an empty queries list
		queriesList <- NULL	
		if (.checkValid(session, schemaName, out)==FALSE) {
			queries <- unique(out$queryName)	

			for (qy in queries)
			{	
				queryObjName <- as.character(qy)
				queriesList[[queryObjName]] <- .LabkeyQuery(x=.getQueryDetails(session, schemaName, qy), schemaName=schemaName
					, queryName=queryObjName)
			}
		}		
		return(queriesList)
	}

	.getQueryDetails <- function(session, schemaName, queryName)
	{	
		## get the default view info as this is what selectRows returns
		out <- labkey.getDefaultViewDetails(attr(session,"baseUrl"), getFolderPath(session), schemaName, queryName)	
		fieldList <- list()
		for (r in row.names(out))
		{				
			rowList <- as.list(out[r,-1])
			fieldList[[as.character(rowList$fieldName)]]<- .LabkeyField(rowList, schemaName, queryName)

		}
		return(fieldList)
	}


	.getSchemasAsList <- function(session)
	{
		out <- .LabkeySchemaList(NULL)
		skey=attr(session, "skey") 
		schemas <- .lksession$skey$validSchemas[,1]
		for (n in schemas) {out[[n]] <- n}
		return(out)
	}


	.checkValid <- function(session, schemaName, queriesDF)
	{	
		## check to see if the schema cache slots we are about to populate are identical down to the queryName level
		## a return of false on any query causes the whole schema to reload
		lksite<- attr(session, "lksite")
		sch <- lksite[[schemaName]]
		retval <- FALSE
		if(length(sch)==0) {}
		else if((length(sch)==1) && is.na(sch)) {}
		else {
			qnames <- names(sch)
			cachednames <- unique(queriesDF$queryName)
			if (length(qnames)!=length(cachednames)) {}
			else {
			 	if(all(qnames %in% cachednames))
			 		{retval<-TRUE}
			 }	
		}
		return (retval)
	}



	###########################################################################
	## Class declarations and print methods 
	##definition of a 'LabkeySession' 
	.LabkeySession <-
	    function(baseUrl, folderPath, lksite, skey)
	{
		 structure(list(NA), class="LabkeySession", baseUrl=baseUrl, folderPath=folderPath, lksite=lksite, skey=skey)
	}

	print.LabkeySession <-
	    function(x, ..., pad="")
	{
		cat("Base Url:  "	, attr(x,"baseUrl"), "\n")
		cat("Folder Path:  ", getFolderPath(x), "\n")
		cat("Available schemas: \n")			
		lsSchemas(x)

		cat("Available folders in current project:  \n\t")
		for (i in lsFolders(x)) {
			cat(as.character(i),"\n\t",sep="")
		}	
		cat("\nAvailable projects on server:  \n\t")
		for (i in lsProjects(attr(x,"baseUrl"))) {
			cat(as.character(i),"\n\t",sep="")
		}	
		cat("\n")
	}
	

	.LabkeySchemaList <-
	    function(x, ...)
	{
	    structure(as.list(x, ...), class="LabkeySchemaList")  
	}

	print.LabkeySchemaList <-
	    function(x,...)
	{  
		i = 0
		for (n in names(x)) {
			i=i+1
			cat(i,"\t", n, "\n",sep="")
		}
		cat("\n")
	}



	.LabkeySchema <-
	    function(x, sname, ...)
	{
	    structure(as.list(x, ...), class="LabkeySchema", schemaName=sname )  
	}

	print.LabkeySchema <-
	    function(x,...)
	{  
		cat("Schema: ", attr(x, "schemaName"), "\n")
		cat("Available queries: \n")
		i=1
		for (nm in names(x)) {
			cat(i,"\t",nm,"\n", sep="")
			i=i+1
		}
		cat("\n")
	}

	## wrap a list of LabKeyFields as a LabkeyQuery
	.LabkeyQuery <-
	    function(x, schemaName, queryName)
	{
	     structure(as.list(x), class="LabkeyQuery", schemaName=schemaName, queryName=queryName)
	}

	print.LabkeyQuery <-
	    function(x, ..., pad="")
	{		
		alldf <- as.data.frame(NULL, nrow=1, ncol=5)
		for (i in names(x)) 
		{    	
			keyinfo <- " "
			refinfo <- " "
			if (x[[i]]$isKeyField==TRUE  &&  !any(grepl("/", x[[i]]$fieldName, fixed=TRUE)))
			{
				keyinfo<- "PK" 
			}    		
			if (!is.na(x[[i]]$lookupQueryName)) 
			{
				refinfo <- paste(refinfo, "lookup to ",x[[i]]$lookupSchemaName,".", x[[i]]$lookupQueryName, sep="") 
			}

			outdf <- data.frame(cbind(x[[i]]$fieldName, x[[i]]$caption, x[[i]]$type, keyinfo, refinfo))
			alldf <- rbind(alldf, outdf)		
		}
		colnames(alldf) <- c("fieldName", "caption", "type", "key", "related query")
		print(alldf, right=FALSE)
	}

	.LabkeyField <-
	    function(x,schemaName, queryName)
	{
	    structure(as.list(x),
		      class="LabkeyField", schemaName=schemaName, queryName=queryName)
	}

	#print.LabkeyField <-
	#    function(x, ..., pad="")
	#{
	#	print(x)
	#}


## the following two functions add an entry for the Rlabkey Users guide to the Vignettes menu
RlabkeyUsersGuide <- function(view=TRUE)
{
     f <- system.file("doc", "usersguide.pdf", package = "Rlabkey")
    if (view) {
        if (.Platform$OS.type == "windows") 
            shell.exec(f)
        else system(paste(Sys.getenv("R_PDFVIEWER"), f, "&"))
    }
    return(f)
}

.onLoad<- function(libname, pkgname)
{	
	try(utils::winMenuAddItem("Vignettes", "Rlabkey", "RlabkeyUsersGuide()"), TRUE)
}

.onUnload<- function(libpath)
{	
	numItems= as.integer(0);
	
	try (numItems<- length(utils::winMenuItems("Vignettes")), TRUE)
	if (numItems > 0)
	{
		try(utils::winMenuDelItem("Vignettes", "Rlabkey"), TRUE)
	}	
}
