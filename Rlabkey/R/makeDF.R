##
#  Copyright (c) 2008-2014 Fred Hutchinson Cancer Research Center
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

 makeDF <- function(rawdata, colSelect=NULL, showHidden, colNameOpt)
{

    decode <- fromJSON(rawdata)

	## Check for invalid colSelect name (with labkey 8.3 this returns lsid column only)
	if(is.null(colSelect)==FALSE){
	if(length(decode$columnModel)==1 & decode$columnModel[[1]]$header!=colSelect) 
		{stop(paste('The column names in the query "',decode$queryName,'" do not match one or more of the names specified in the colSelect variable. Be sure you are using the column name and not the column label. See the documentation for more details.',sep=''))}}


  	## Get column names in proper order, associated header index, hidden tag, and data type
  	cnames <- NULL
  	hindex <- NULL
  	hide <- NULL
  	for(j in 1:length(decode$columnModel))
  	{   	
		## three different ways to refer to columns exist today
		## selectRows and executeSQL by default return the field caption (also called the  "label")
		## When sepcifying colSelect, colFilter, or ExecuteSql, you must use field names.
		## when running R "views"  at the server, the field_name is modified to use underscores and lower cased.  
		## This also makes it a legal name in R, which can be useful.

  		if (colNameOpt=="caption")
  			{cname <- decode$columnModel[[j]]$header}
  		else if (colNameOpt == "fieldname")
  			{cname <- decode$columnModel[[j]]$dataIndex}
  		else if (colNameOpt == "rname" )
  			{cname <- .getRNameFromName(decode$columnModel[[j]]$dataIndex, existing=cnames) }
  		else
  			{ stop("Invalid colNameOpt option.  Valid values are caption, fieldname, and rname.") }
  			
		cnames <- c(cnames, cname)
  	        hindex <- c(hindex, decode$columnModel[[j]]$dataIndex)
  	        hide <- c(hide, decode$columnModel[[j]]$hidden)}
  	refdf <- data.frame(cnames,hindex,hide)
  	
	## Check for no rows returned, put data in data frame 
  	if(length(decode$rows)<1)
		{tohide <- length(which(refdf$hide==TRUE))
		totalcol <- length(refdf$cnames)
		if(showHidden==FALSE)
			{emptydf <- as.data.frame(rep(list(num=double(0)), each=(totalcol-tohide)))
			colnames(emptydf) <- refdf$cnames[refdf$hide==FALSE]
			warning("Empty data frame was returned. Query may be too restrictive.", call.=FALSE)
			return(emptydf)}else
		{emptydf <- as.data.frame(rep(list(num=double(0)), each=(totalcol)))
		colnames(emptydf) <- refdf$cnames
		warning("Empty data frame was returned. Query may be too restrictive.", call.=FALSE)}
		return(emptydf)}
		
	if(length(decode$rows)>0) {
		hold.dat <- NULL
    		hold.dat <- matrix(sapply(sapply(decode$rows,filterrow),rbind), nrow=length(decode$rows), byrow=TRUE)
    		hold.dat <- as.data.frame(hold.dat,stringsAsFactors=FALSE)
		tmprow <- filterrow(decode$rows[[1]])
		names(hold.dat) <- names(tmprow)   			
	}

	## Order data
	oindex <- NULL
	## number of cols selected may be more or less than described in metadata
  	for(k in 1:length(cnames)){oindex <- rbind(oindex, which(names(hold.dat)==refdf$hindex[k]))}

  	refdf$oindex <- oindex
	newdat <- as.data.frame(hold.dat[,refdf$oindex],stringsAsFactors=FALSE)

  	refdf$type <- NULL
  	for(p in 1:dim(refdf)[1]) {   
  	    ind <- which(refdf$hindex==decode$metaData$fields[[p]]$name)
  	    refdf$type[ind] <- decode$metaData$fields[[p]]$type
  	}

  	## Delete hidden column(s) unless showHidden=TRUE
	if (showHidden==TRUE || is.null(decode$metaData$id)) {}  else { 
		hide.ind <- which(refdf$hide==TRUE)
		if(length(hide.ind)>0 ) {
			if(length(hide.ind) == length(newdat)) {
				stop("No visible columns selected.  Use the showHidden=TRUE to see these columns.") 
			} else {
				newdat <- newdat[,-hide.ind]
				refdf <- refdf[-hide.ind,]
				cnames <- cnames[-hide.ind] 
			}
		} else {}		
	}

	## Set mode for multiple columns of data (this also removes list factor)
	if(is.null(dim(newdat))==FALSE)
  	{
  	    for(j in 1:ncol(newdat))
  	    {
  	        mod <- refdf$type[j]
  	        try(
                if(mod=="date") { newdat[,j] <- .parseDate(newdat[,j])} else
                if(mod=="string"){	suppressWarnings(mode(newdat[,j]) <- "character")} else
                if(mod=="int"){ suppressWarnings(mode(newdat[,j]) <- "numeric")} else
                if(mod=="boolean"){suppressWarnings(mode(newdat[,j]) <- "logical")} else
                if(mod=="float"){suppressWarnings(mode(newdat[,j]) <- "numeric")} else
                {print("MetaData field type not recognized.")}
            , silent=TRUE)
        }
	    newdat <- as.data.frame(newdat, stringsAsFactors=FALSE); colnames(newdat)<-cnames
    }

	## Set mode for single column of data
	if(is.null(dim(newdat))==TRUE & length(newdat)>1)
	{
	    mod <- refdf$type
	    try(
            if(mod=="date"){ newdat <- .parseDate(newdat)}else
            if(mod=="string"){suppressWarnings(mode(newdat) <- "character")} else
            if(mod=="int"){ suppressWarnings(mode(newdat) <- "numeric")} else
            if(mod=="boolean"){suppressWarnings(mode(newdat) <- "logical")} else
            if(mod=="float"){suppressWarnings(mode(newdat) <- "numeric")} else
            {print("MetaData field type not recognized.")}
        , silent=TRUE)
	    newdat <- as.data.frame(newdat, stringsAsFactors=FALSE); colnames(newdat)<-cnames[1]
    }

return(newdat)
}

## need to get rid of hidden hrefs within the row, R doesn't use them and their presence causes problems
## also consolidate null handling here
filterrow<-function(row)
{
	filtered <- NULL
	for (x in 1:length(row)) {		
		valname <- names(row[x])
		if ((nchar(valname)>11) && (substr(valname,1,11) == as.character("_labkeyurl_"))) {
			next
		}
		if (is.null(row[x][[valname]])) { row[x][[valname]]<-NA }
		filtered <- c(filtered, row[x])
	}
return(filtered)

}

.getRNameFromName <- function(lkname, existing=NULL)
{
	rname <- tolower(chartr(" /", "__", lkname))
	
	if (length(existing)>0) 
	{ 
		for (i in 1:99)
		{
			if(length(existing[rname == existing]) ==0 )
				{break;}
			else 
				{rname<- c(rname + as.character(i))}
        } 
  	}    	
  	return (rname)
}

.parseDate <- function(s)
{
    s <- as.character(s);
    d <- tryCatch(as.Date(s),error = function(e) NA);
    if (any(is.na(d)))
    {
        d[is.na(d)] <- as.Date(s[is.na(d)], "%d %b %Y %H:%M:%S");
    }
    return(d);
}

# converts the factors of a data frame to characters
convertFactorsToStrings <- function(df)
{
    factors <- sapply(df, is.factor)
    df[factors] <- lapply(df[factors], as.character)
    return(df);
}