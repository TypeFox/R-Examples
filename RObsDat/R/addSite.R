addSite <- function(Code, Name, x, y, Elevation=rep(0, length(Code)), LatLongDatum, LocalProjection=NULL, isLocal=NULL, VerticalDatum=NULL,PositionAccuracy=rep(0, length(Code)), State=NULL, County=NULL, Comment=NULL){
	#mandatory fields: Code, Name, Lat, Long, LatLongDatum

	#check from reference tables: SpatialReferences -> LatLongDatum, LocalProjection
	stopifnot(length(LatLongDatum) == length(Code))
	stopifnot(length(Name) == length(Code))
	stopifnot(length(x) == length(Code))
	stopifnot(length(y) == length(Code))
	stopifnot(length(Elevation) == length(Code))
	stopifnot(length(LatLongDatum) == length(Code))
	stopifnot(length(PositionAccuracy) == length(Code))


	if(!is.null(isLocal)){
		stopifnot(length(isLocal) == length(Code))
	}
	if(!is.null(State)){
		stopifnot(length(State) == length(Code))
	}
	if(!is.null(County)){
		stopifnot(length(County) == length(Code))
	}
	if(!is.null(Comment)){
		stopifnot(length(Comment) == length(Code))
	}

	#Check for existing entries
	theUnexisting <- rep(TRUE, length(Code))
	if(NROW(existing <- getMetadata("Site",Name=Name, exact=TRUE))>0){
		warning(paste("Skiping existing Site(s):", paste(existing[,3], collapse="; "), "\n"))
		theUnexisting <- !(Name %in% existing[,3])

	}

	if(!any(theUnexisting)){
		return()
	}


	SpatialReferenceID <- getID("SpatialReference", LatLongDatum)


	if(!is.null(LocalProjection)){
		stopifnot(length(LocalProjection) == length(Code))
		SpatialReferenceID2 <- getID("SpatialReference", LocalProjection)
	} else {
		#ToDo
		SpatialReferenceID2 <- rep(1, length(Code))
	}

	#check from reference tables: VerticalDatum
	if(!is.null(VerticalDatum)){
		stopifnot(length(VerticalDatum) == length(Code))
		VertDatID <- getID("VerticalDatum",VerticalDatum)
	} else {
		VertDatID <- rep(getID("VerticalDatum", "Unknown"), length(Code))
	}

	#transform coordinates
	#depending on value of isLocal
	todo("ToDo automatic coordinate conversion")

	IaddSite(getOption("odm.handler"),Code = Code[theUnexisting], Name =Name[theUnexisting], Latitude=x[theUnexisting], Longitude=y[theUnexisting], Elevation=Elevation[theUnexisting], LatLongDatum=SpatialReferenceID[theUnexisting], LocalProjection=SpatialReferenceID2[theUnexisting], VerticalDatum=VertDatID[theUnexisting],PosAccuracy=PositionAccuracy[theUnexisting], State=State[theUnexisting], County=County[theUnexisting], Comments=Comment[theUnexisting])


}
