##
## setup classes that will be used by the package
##

require(RJSONIO)
require(RCurl)

## - - - - - - - - - - - - - - - - - - - -

setClass( "FactualConnection" ,
	representation(
		apiKey="character" ,
		baseURL="character"
	)
) ;

## sic - instead of a setMethod( "initialize" ... ) we will use 
## a wrapper function to set the defaults.  The wrapper lets 
## end-users create the object as a simple function call (aka,
## without having to deal with OO syntax) and by not using
## an initializer, we keep all of the default values in one
## place.

createFactualConnection <- function( apiKey = NULL , apiVersion = NULL , baseURL = NULL ){

	## fail if apiKey is NA, NULL, or blank
	if( is.null( apiKey ) || "" == apiKey ){
		stop( "Please provide your factual API key as the 'apiKey' parameter" )
	}

	## confirm that at least one of apiVersion and baseURL is NA
	if( ! is.null( apiVersion ) && ! is.null( baseURL ) ){
		stop( "Please define only one of 'apiVersion' or 'baseURL'" )
	}

	result <- new( "FactualConnection" )

	result@apiKey <- apiKey

	if( ! is.null( baseURL ) ){
		result@baseURL <- baseURL
	}else if( ! is.null( apiVersion ) ){
		result@baseURL <- paste( "http://api.factual.com/" , apiVersion , sep="" )
	}else{
		result@baseURL <- "http://api.factual.com/v2"
	}

	return( result )
}


## - - - - - - - - - - - - - - - - - - - -

setClass( "FactualSchemaResult" ,
	representation(
		url="character" ,
		tableID="character" , ## called "key" in output from Factual
		name="character" ,

		status="character" ,
		message="character" ,

		table.meta="list",

		fields="data.frame" ,

		processingTime="numeric" ,
		runTime="numeric" ,
		fetchedAt="POSIXct"
	) ,

	prototype(
		fetchedAt=Sys.time()
	)

) ;


processFactualSchemaJSON <- function( factual.output ){

	result <- new( "FactualSchemaResult" )

	result@status <- factual.output$status

	if( ! is.null( factual.output$message ) ){
		result@message <- factual.output$message
	}

	if( "error" == result@status ){
		warning( "Factual returned an error; there will be no data.  Please see @message." )
	}else{

		result@tableID <- factual.output$schema$key
		result@name <- factual.output$schema$name

		result@table.meta <- list(
			isPrivate = factual.output$schema$isPrivate ,
			geoEnabled = factual.output$schema$geoEnabled ,
			totalRowCount = factual.output$schema$totalRowCount ,
			rating = factual.output$schema$rating ,
			createdByUserId = factual.output$schema$createdByUserId ,
			sourceTables = factual.output$schema$sourceTables ,
			updatedAt <- as.POSIXct( factual.output$schema$updatedAt ) ,
			createdAt <- as.POSIXct( factual.output$schema$createdAt ) ,
			description = factual.output$schema$description ,
			isDownloadable = factual.output$schema$isDownloadable ,
			source = factual.output$schema$source ,
			creator = factual.output$schema$creator ,
			fieldRefs = factual.output$schema$fieldRefs
		)

		result@fields <- as.data.frame( do.call( rbind , factual.output$schema$fields ) )
	}

	return( result )

}


factualGetSchema <- function( connection , tableID , verbose=FALSE ){

	##	URL format:
	##	http://api.factual.com/v2/tables/TABLE_ID/schema?api_key=API_KEY
	factual.url <- paste( connection@baseURL , "tables" , tableID , "schema" , sep="/" )
	factual.url <- paste( factual.url , '?api_key=' , connection@apiKey , sep="" )

	if( verbose ){
		cat( "URL: \"" , factual.url , "\"" , "\n" , sep="" )
	}

	time.start <- Sys.time()
	factual.output <- fromJSON( factual.url )
	time.end <- Sys.time()

	result <- processFactualSchemaJSON( factual.output )
	result@url <- factual.url
	result@runTime <- as.numeric( time.end - time.start )

	return( result )

}

## - - - - - - - - - - - - - - - - - - - -

setClass( "FactualReadResult" ,
	representation(
		url="character" ,

		tableID="character" ,
		status="character" ,
		message="character" ,

		resultRows="numeric" , ## called "rows" in Factual output
		tableRows="numeric" , ## called "total_rows" in Factual output

		results="data.frame" ,
		
		processingTime="numeric" ,
		runTime="numeric" ,
		fetchedAt="POSIXct"
	) ,

	prototype(
		fetchedAt=Sys.time()
	)
) ;


processFactualReadJSON <- function( factual.output , reformatNames=FALSE ){

	result <- new( "FactualReadResult" )

	if( ! is.null( factual.output$message ) ){
		result@message <- factual.output$message
	}


	result@status <- factual.output$status

	if( "error" == result@status ){
		warning( "Factual returned an error; there will be no data.  Please see @message." )
	}else{
		result@resultRows <- factual.output$response$rows
		result@tableRows <- factual.output$response$total_rows
		result@results <- as.data.frame( do.call( rbind , factual.output$response$data ) )
		names( result@results ) <- unlist( factual.output$response$fields )

		## keep only alphanumeric characters
		if( reformatNames ){
			names( result@results ) <- gsub( "[^a-zA-Z0-9]" , "" , names( result@results ) )
		}

	}

	return( result )

}


factualRead <- function( connection , tableID , sort=NULL , limit=NULL , offset=NULL , filters=NULL , subject_key=NULL , reformatNames=FALSE , verbose=FALSE ){

	##
	## URL format:
	##	http://api.factual.com/v2/tables/TABLE_ID/read?api_key=API_KEY
	## 
	## plus, any of the following params:
	## 	sort=[{SORT_JSON}]
	## 	limit=LIMIT_NUM
	## 	offset=OFFSET_NUM
	## 	filters={FILTERS_JSON}
	## 	subject_key=SUBJECT_KEY_STRING
	##

	factual.url <- paste( connection@baseURL , "tables" , tableID , "read" , sep="/" )
	factual.url <- paste( factual.url , '?api_key=' , connection@apiKey , sep="" )


	if( ! is.null( sort ) ){
		stopifnot( is.character( sort ) )
		factual.url <- paste( factual.url , '&sort=' , curlEscape( sort ) , sep="" )
	}

	if( ! is.null( limit ) ){
		stopifnot( is.numeric( limit ) )
		factual.url <- paste( factual.url , '&limit=' , limit , sep="" )
	}

	if( ! is.null( offset ) ){
		stopifnot( is.numeric( offset ) )
		factual.url <- paste( factual.url , '&offset=' , offset , sep="" )
	}

	if( ! is.null( filters ) ){
		stopifnot( is.character( filters ) )
		factual.url <- paste( factual.url , '&filters=' , curlEscape( filters ) , sep="" )
	}

	if( ! is.null( subject_key ) ){
		stopifnot( is.character( subject_key ) )
		factual.url <- paste( factual.url , '&subject_key=' , curlEscape( subject_key ) , sep="" )
	}

	if( verbose ){
		cat( "URL: \"" , factual.url , "\"" , "\n" , sep="" )
	}

	time.start <- Sys.time()
	factual.output <- fromJSON( factual.url )
	time.end <- Sys.time()

	result <- processFactualReadJSON( factual.output , reformatNames )
	result@url <- factual.url
	result@tableID <- tableID
	result@runTime <- as.numeric( time.end - time.start )

	return( result )

}


## - - - - - - - - - - - - - - - - - - - -
