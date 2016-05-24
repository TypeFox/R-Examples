
readWKT = function( text, id = NULL, p4s = NULL) {
    
    #wkt = str_replace(text,"\n","")
    wkt = gsub("\n", "", text)

    if (length(wkt) != 1) stop("WKT must have length 1")
    
	pat = 
	"POINT|LINESTRING|LINEARRING|POLYGON|MULTIPOINT|MULTILINESTRING|MULTIPOLYGON"
    # m =  str_extract_all(wkt, pat)
	m = strsplit(wkt, pat)
    # ngeoms =  length( m[[1]] )
    ngeoms =  length( m[[1]] ) - 1 # "" before first match is always the first
    if(is.null(ngeoms)) ngeoms = 0
    
    if(is.null(id)) {
        if (ngeoms == 0) {
            id = c()
        } else {
            id = 1:ngeoms
        }
    } 
    
    # if the number of ids doesn't take into account sub geometries in geometry collection then create subids
    if( length(id) == 1 & ngeoms != 1)
        id = paste(id,1:ngeoms,sep=".")
    
    if( length(id) != ngeoms )
        stop("number of WKT geometries does not match number of ids")        

    p4s = checkP4S(p4s)
    
    id = as.character(id)
        
    tryCatch(res <- .Call("rgeos_readWKT", .RGEOS_HANDLE, wkt, p4s, 
                            id, PACKAGE="rgeos"), 
             error = function(e) { stop( paste( "Unable to parse: ",wkt,"\n",
                                                "GEOS reported: \"", e$message,"\"",sep=""),call.=FALSE) } )
    
	#if (length(unique(row.names(res))) != ngeoms)
	#if (!is.null(res)) {
	#	if ( sum(sapply(row.names(res),length)) != ngeoms )
    #    	warning(paste("Number of geometries does not match between object and text, check WKT validity.",wkt),call.=FALSE)
    #}

    return( res )
}


writeWKT = function( spgeom, byid = FALSE) {

    stopifnot(is.logical(byid))
    byid = as.logical(byid)
    
    res <- .Call("rgeos_writeWKT", .RGEOS_HANDLE, spgeom, byid)
    
    return(res) 
}
