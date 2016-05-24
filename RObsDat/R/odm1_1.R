setClass("odm1_1",
	representation= representation(con = "DBIConnection", checked = "logical"),
	prototype = list(con=NULL, checked=FALSE)
	)
setClass("odm1_1Ver", contains="odm1_1"	)


setGeneric("IdbState", function(object) {standardGeneric("IdbState")})
setGeneric("IgetSite", function(object, ID=NULL, Code=NULL, Name=NULL, x=NULL, y=NULL, Elevation=NULL, LatLongDatum=NULL, exact=FALSE ) { standardGeneric("IgetSite")}) 
setGeneric("IgetUnits", function(object, ID=NULL, Name=NULL, Type=NULL, Abbreviation=NULL, exact=FALSE ) { standardGeneric("IgetUnits")}) 
setGeneric("IaddUnits", function(object, Name, Type, Abbreviation ) { standardGeneric("IaddUnits")}) 
setGeneric("IaddSpatialReference", function(object, ID, SRSName, SRSID, IsGeographic, Notes ) { standardGeneric("IaddSpatialReference")}) 
setGeneric("IgetVariable", function(object, ID=NULL, Code=NULL, Name=NULL, Speciation=NULL, Unit=NULL, Medium=NULL,exact=FALSE, ...  ) { standardGeneric("IgetVariable")}) 
setGeneric("IgetQualifier", function(object, ID=NULL, Code=NULL, Description=NULL, ...  ) { standardGeneric("IgetQualifier")}) 
setGeneric("IgetLabMethod", function(object, ID=NULL, LabName=NULL, Organization=NULL, MethodName=NULL, Description=NULL, Link=NULL, ...  ) { standardGeneric("IgetLabMethod")}) 
setGeneric("IgetMethod", function(object, ID=NULL, Description=NULL, ...  ) { standardGeneric("IgetMethod")}) 
setGeneric("IgetOffsetType", function(object, ID=NULL, Description=NULL, Units=NULL, ...  ) { standardGeneric("IgetOffsetType")}) 
setGeneric("IaddOffsetType", function(object, Description=NULL, Units=NULL, ...  ) { standardGeneric("IaddOffsetType")}) 
setGeneric("IgetSample", function(object, ID=NULL, LabSampleCode=NULL, ...  ) { standardGeneric("IgetSample")}) 
setGeneric("IgetSource", function(object, ID=NULL, Organization=NULL, Description=NULL, Citation=NULL, ...  ) { standardGeneric("IgetSource")}) 
setGeneric("IgetISOMetadata", function(object, ID=NULL, Title=NULL, Abstract=NULL, TopicCategory=NULL, ...  ) { standardGeneric("IgetISOMetadata")}) 
setGeneric("IgetQualityControlLevel", function(object, ID=NULL, Code=NULL, Definition=NULL, Explanation=NULL, ...  ) { standardGeneric("IgetQualityControlLevel")}) 
setGeneric("IaddQualityControlLevel", function(object, ID, Code, Definition, Explanation, ...  ) { standardGeneric("IaddQualityControlLevel")}) 
setGeneric("IaddLabMethod", function(object, ID, LabName, Organization, MethodName, Description, Link, ...  ) { standardGeneric("IaddLabMethod")}) 

for(i in CVtables()){
	code <- paste('setGeneric("Iget',i,'", function(object, Term=NULL,  Definition=NULL, exact=FALSE, ...  ) { standardGeneric("Iget',i,'")})', sep="")
	eval(parse(text=code))
}
for(i in CVtables()){
	#ToDo: does not work for VariableName
	code <- paste('setMethod("Iget',i,'", signature(object = "odm1_1"), function(object,  Term=NULL,  Definition=NULL, exact=FALSE, ...){ return(IgetCV(object, "',i,'", Term, Definition, exact=exact)) })', sep="")
	eval(parse(text=code))
}
for(i in CVtables()){
	code <- paste('setMethod("Iget',i,'", signature(object = "NULL"), h.m)', sep="")
	eval(parse(text=code))
}

 
setGeneric("IgetDataValues", function(object,ID=NULL, from=NULL, to=NULL, tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, ...  ) { standardGeneric("IgetDataValues")}) 
setGeneric("IgetSpatialReference", function(object,ID=NULL, SRSID=NULL, SRSName=NULL, IsGeographic=NULL, Notes=NULL, exact=FALSE) { standardGeneric("IgetSpatialReference")}) 
setGeneric("IgetOldDataValues", function(object,ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, VersionID, exact,...  ) { standardGeneric("IgetOldDataValues")}) 
setGeneric("IgetDeletedDataValues", function(object,ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, VersionID, exact,...  ) { standardGeneric("IgetDeletedDataValues")}) 
setGeneric("IaddDataValues", function(object, localDateTime, values, TZ, SiteID, VariableID, Offset=NULL, OffsetTypeID=NULL, CensorCode, QualifierID=NULL, MethodID, SourceID, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID, exact, ...  ) { standardGeneric("IaddDataValues")}) 
setGeneric("IaddSite", function(object, Code, Name, Latitude, Longitude, LatLongDatum, Elevation=NULL, VerticalDatum=NULL, LocalX=NULL,LocalY=NULL, LocalProjection=NULL, PosAccuracy=NULL, State=NULL, County=NULL, Comments=NULL) {standardGeneric("IaddSite")})
setGeneric("IaddVariable", function(object, Code, Name, Speciation, Unit, SampleMedium,ValueType, IsRegular, TimeSupport, TimeUnits, DataType, GeneralCategory, NoDataValue) {standardGeneric("IaddVariable")})
setGeneric("IaddSource", function(object, Organization, SourceDescription, SourceLink, ContactName, Phone, Email, Address, City, State, ZipCode, Citation, Metadata) {standardGeneric("IaddSource")})
setGeneric("IaddISOMetadata", function(object, TopicCategory, Title, Abstract, ProfileVersion, MetadataLink) {standardGeneric("IaddISOMetadata")})
setGeneric("IupdateDataValues", function(object, ValueID, localDateTime, value, TZ, SiteID, VariableID, Offset=NULL, OffsetTypeID=NULL, CensorCode, QualifierID=NULL, MethodID, SourceID, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID, ...  ) { standardGeneric("IupdateDataValues")}) 
setGeneric("IdeleteDataValues", function(object, ValueID) { standardGeneric("IdeleteDataValues")}) 
setGeneric("IarchiveDataValues", function(object, ValueID,reason) { standardGeneric("IarchiveDataValues")}) 
setGeneric("IaddDataVersion", function(object, reason) { standardGeneric("IaddDataVersion")}) 
setGeneric("IgetCurrentDataVersion", function(object) { standardGeneric("IgetCurrentDataVersion")}) 
setGeneric("IgetDataVersions", function(object) { standardGeneric("IgetDataVersions")}) 
setGeneric("IaddCV", function(object, table, term, definition){standardGeneric("IaddCV")})
setGeneric("IgetCV", function(object, table, term, definition, exact=FALSE){standardGeneric("IgetCV")})
setGeneric("IaddSynonym", function(object, phrase, table, id){standardGeneric("IaddSynonym")})
setGeneric("IgetSynonymID", function(object, phrase, table){standardGeneric("IgetSynonymID")})
setGeneric("IgetSynonyms", function(object, phrase, table){standardGeneric("IgetSynonyms")})

setGeneric("IgetNo", function(object, table) {standardGeneric("IgetNo")})

#copy set generic from above
# mark and replace with
# '<,'>s/setGeneric/setMethod/g
# :'<,'>s/,.*/, signature(object = NULL), h.m)   

setMethod("IgetCV", signature(object = "NULL"), h.m)
setMethod("IgetSynonymID", signature(object = "NULL"), h.m)
setMethod("IgetSynonyms", signature(object = "NULL"), h.m)
setMethod("IdbState", signature(object = "NULL"), h.m)
setMethod("IaddCV", signature(object = "NULL"), h.m)
setMethod("IgetSite", signature(object = "NULL"), h.m)
setMethod("IgetUnits", signature(object = "NULL"), h.m)
setMethod("IaddUnits", signature(object = "NULL"), h.m)
setMethod("IaddSpatialReference", signature(object = "NULL"), h.m)
setMethod("IgetVariable", signature(object = "NULL"), h.m)
setMethod("IgetQualifier", signature(object = "NULL"), h.m)
setMethod("IgetLabMethod", signature(object = "NULL"), h.m)
setMethod("IgetMethod", signature(object = "NULL"), h.m)
setMethod("IgetOffsetType", signature(object = "NULL"), h.m)
setMethod("IaddOffsetType", signature(object = "NULL"), h.m)
setMethod("IgetSample", signature(object = "NULL"), h.m)
setMethod("IgetSource", signature(object = "NULL"), h.m)
setMethod("IgetISOMetadata", signature(object = "NULL"), h.m)
setMethod("IgetQualityControlLevel", signature(object = "NULL"), h.m)
setMethod("IaddQualityControlLevel", signature(object = "NULL"), h.m)
setMethod("IgetCensorCode", signature(object = "NULL"), h.m)
setMethod("IgetDataValues", signature(object = "NULL"), h.m)
setMethod("IgetSpatialReference", signature(object = "NULL"), h.m)
setMethod("IgetOldDataValues", signature(object = "NULL"), h.m)
setMethod("IgetDeletedDataValues", signature(object = "NULL"), h.m)
setMethod("IaddDataValues", signature(object = "NULL"), h.m)
setMethod("IaddSite", signature(object = "NULL"), h.m)
setMethod("IaddVariable", signature(object = "NULL"), h.m)
setMethod("IupdateDataValues", signature(object = "NULL"), h.m)
setMethod("IdeleteDataValues", signature(object = "NULL"), h.m)
setMethod("IarchiveDataValues", signature(object = "NULL"), h.m)
setMethod("IaddDataVersion", signature(object = "NULL"), h.m)
setMethod("IgetCurrentDataVersion", signature(object = "NULL"), h.m)
setMethod("IgetDataVersions", signature(object = "NULL"), h.m)
setMethod("IgetNo", signature(object = "NULL"), h.m)
setMethod("IaddSource", signature(object = "NULL"), h.m)
setMethod("IaddISOMetadata", signature(object = "NULL"), h.m)

check.version <- function(object, version){
		if(getOption("odm.handler")@checked==TRUE) return()

		if(mdbExistsTable(object@con, "ODMVersion")){
			query = "SELECT * FROM ODMVersion"
		        if(getOption("verbose.queries", default=FALSE)) print(query)
			res <- dbGetQuery(object@con, query)
			if(NROW(res)==0){
				stop("No Database Version stored in table ODMVersion. Please specify the version (e.g. 'dbGetQuery(options('odm.handler')[[1]]@con, 'INSERT INTO ODMVersion (VersionNumber) values (\"1.1Ver\")') ).\n")
			}
			if(res[,1]!=version){
				stop("Invalid Database version. Expected ",version,", obtained ",res[,1])
			}

		} else {
			warning("Creating database structure")
			run.sql.script(object@con, system.file("odm1_1_raw.sql", package="RObsDat"))
			if(version=="1.1Ver"){
			     run.sql.script(object@con, system.file("odm1_1_addVersion.sql", package="RObsDat"))
			     todo("Fix problem with constraints for DataValuesRepository")
			}
			warning("Updating controlled vocabularies")
			updateCV()
		}
		if(!mdbExistsTable(object@con, "Synonyms")){
			warning("Creating synonym table")
			query = "CREATE TABLE Synonyms (phrase TEXT, RecID NUMERIC, tab TEXT)"
		        if(getOption("verbose.queries", default=FALSE)) print(query)
			dbGetQuery(object@con, query)
		}
		han <- getOption("odm.handler")
		han@checked <- TRUE
		options(odm.handler=han)
}
 
setMethod("IdbState", signature(object = "odm1_1"), 
	function(object){
		check.version(object, "1.1")
	}
)
setMethod("IdbState", signature(object = "odm1_1Ver"), 
	function(object){
		check.version(object, "1.1Ver")
	}
)
setMethod("IgetISOMetadata", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, Title=NULL, Abstract=NULL, TopicCategory=NULL, exact=FALSE){
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "MetadataID", exact=TRUE)
			w.o <- expand.where(w.o, Title, "Title", exact=exact)
			w.o <- expand.where(w.o, Abstract, "Abstract", exact=exact)
			w.o <- expand.where(w.o, TopicCategory, "TopicCategory", exact=TRUE)

			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)

			query <- paste("SELECT * FROM ISOMetadata", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1] <- c("ID")}
			return(res)
		}
	)


setMethod("IgetUnits", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, Name=NULL, Type=NULL, Abbreviation=NULL, exact=FALSE){
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "UnitsID", exact=TRUE)
			w.o <- expand.where(w.o, Name, "UnitsName", exact=exact)
			w.o <- expand.where(w.o, Type, "UnitsType", exact=exact)
			w.o <- expand.where(w.o, Abbreviation, "UnitsAbbreviation", exact=exact)

			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)

			query <- paste("SELECT * FROM Units", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(!is.null(res) & NCOL(res)>0) {names(res)[1:4] <- c("ID", "Name", "Type", "Abbreviation")}
			return(res)
		}
	)


setMethod("IgetSite", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, Code=NULL, Name=NULL, x=NULL, y=NULL, Elevation=NULL, LatLongDatum=NULL, exact=FALSE){
			#the where.object is used to assemble the where clause
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "SiteID", exact=TRUE)
			w.o <- expand.where(w.o, Code, "SiteCode", exact=exact)
			w.o <- expand.where(w.o, Name, "SiteName", exact=exact)

			if(!is.null(x)){
				stop("ToDo: Implement getSite for x")
			}
			if(!is.null(y)){
				stop("ToDo: Implement getSite for y")
			}
			if(!is.null(LatLongDatum)){
				stop("ToDo: Implement getSite for LatLongDatum")
			}
			if(!is.null(Elevation)){
				stop("ToDo: Implement getSite for Elevation")
			}

			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)

			query <- paste("SELECT * FROM Sites ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1:3] <- c("ID", "Code", "Name")}
			return(res)
		}
	)


setMethod("IgetVariable", 
		signature(object = "odm1_1"),
		function(object, ID=NULL, Code=NULL, Name=NULL, Speciation=NULL, Unit=NULL, Medium=NULL, exact=FALSE, ... ){
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "VariableID", exact=TRUE)
			w.o <- expand.where(w.o, Code, "VariableCode", exact=exact)
			w.o <- expand.where(w.o, Name, "VariableName", exact=exact)
			if(!is.null(Speciation)){
			    theSpeciation <- getID("Speciation", Speciation)
			    w.o <- expand.where(w.o, theSpeciation, "Speciation", exact=TRUE)
			}
			if(!is.null(Unit)){
			    UnitID <- getID("Units", Unit)
			    w.o <- expand.where(w.o, UnitID, "VariableUnitsID", exact=TRUE)
			}
			if(!is.null(Medium)){
			    theMedium <- getID("Medium", Medium)
			    w.o <- expand.where(w.o, theMedium, "Medium", exact=TRUE)
			}
			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)

			query <- paste("SELECT * FROM Variables ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1:3] <- c("ID", "Code", "Name")}
			return(res)
		}
	)

setMethod("IgetLabMethod", 
		signature(object = "odm1_1"),
		function(object, ID=NULL, LabName=NULL, Organization=NULL, MethodName=NULL, Description=NULL, Link=NULL, ...){
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "LabMethodID", exact=TRUE)
			w.o <- expand.where(w.o, LabName, "LabName", exact=FALSE)
			w.o <- expand.where(w.o, Organization, "LabOrganization", exact=FALSE)
			w.o <- expand.where(w.o, MethodName, "LabMethodName", exact=FALSE)
			w.o <- expand.where(w.o, Link, "LabMethodLink", exact=FALSE)
			w.o <- expand.where(w.o, Description, "LabMethodDescription", exact=FALSE)
			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)
			query <- paste("SELECT * FROM Qualifiers ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1] <- c("ID")}
			return(res)
		}
	)

setMethod("IgetQualifier", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, Code=NULL, Description=NULL, ...){
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "QualifierID", exact=TRUE)
			w.o <- expand.where(w.o, Code, "QualifierCode", exact=TRUE)
			w.o <- expand.where(w.o, Description, "QualifierDescription", exact=FALSE)
			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)
			query <- paste("SELECT * FROM Qualifiers ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1:3] <- c("ID", "Code", "Description")}
			return(res)
		}
	)

setMethod("IgetSpatialReference",
	       signature(object = "odm1_1"),
       	       function(object,ID=NULL, SRSID=NULL, SRSName=NULL, IsGeographic=NULL, Notes=NULL,exact=FALSE) { 

			#the where.object is used to assemble the where clause
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "SpatialReferenceID", exact=TRUE, isnumeric=TRUE)
			w.o <- expand.where(w.o, SRSID, "SRSID", exact=TRUE, isnumeric=TRUE)
			w.o <- expand.where(w.o, SRSName, "SRSName", exact=exact)
			w.o <- expand.where(w.o, IsGeographic, "IsGeographic", exact=TRUE)
			w.o <- expand.where(w.o, Notes, "Notes", exact=exact)

			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)
			
			query <- paste("SELECT * FROM SpatialReferences", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1] <- c("ID")}
			return(res)
	       }) 
setMethod("IgetQualityControlLevel", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, Code=NULL, Definition=NULL, Explanation=NULL, ...){
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "QualityControlLevelID", exact=TRUE)
			w.o <- expand.where(w.o, Code, "QualityControlLevelCode", exact=FALSE)
			w.o <- expand.where(w.o, Definition, "Definition", exact=FALSE)
			w.o <- expand.where(w.o, Explanation, "Explanation", exact=FALSE)
			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)
			query <- paste("SELECT * FROM QualityControlLevels ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1:2] <- c("ID", "Code")}
			return(res)
		}
	)
setMethod("IgetSample", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, LabSampleCode=NULL, ...){
		  	w.o <- list(where.clause = "", the.and  = "")

			w.o <- expand.where(w.o, ID, "SampleID", exact=TRUE)
			w.o <- expand.where(w.o, LabSampleCode, "LabSampleCode", exact=FALSE)
			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)
			query <- paste("SELECT * FROM Samples ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1] <- c("ID")}
			return(res)
		}
	)
setMethod("IaddOffsetType", signature=(object = "odm1_1"),
	  function(object, Description, Units) {
			for(rownum in seq(along=Description)){
				#no Foreign Key
				theDescription <- sv(Description, rownum)
				theUnit <- svk(Units, "Unit", rownum, object)

				insert.query <- paste("INSERT INTO OffsetTypes (OffsetDescription, OffsetUnitsID) VALUES (\"", paste(
						theDescription ,
						theUnit,
						sep="\", \""), "\")", sep="")
			 	run.query(object, insert.query )
			}
	  }
)
setMethod("IgetOffsetType", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, Description=NULL, Units=NULL, ...){
		  	w.o <- list(where.clause = "", the.and  = "")

			w.o <- expand.where(w.o, ID, "OffsetTypeID", exact=TRUE)
			w.o <- expand.where(w.o, Description, "OffsetDescription", exact=FALSE)
			w.o <- expand.where(w.o, Units, "OffsetUnitsID", exact=FALSE)

			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)
			query <- paste("SELECT * FROM OffsetTypes ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1:3] <- c("ID", "UnitsID", "Description")}
			return(res)
		}
	)
setMethod("IgetMethod", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, Description=NULL, ...){
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "MethodID", exact=TRUE)
			w.o <- expand.where(w.o, Description, "MethodDescription", exact=FALSE)
			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)
			query <- paste("SELECT * FROM Methods ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1:3] <- c("ID", "Code", "Description")}
			return(res)
		}
	)
setMethod("IgetSource", 
		signature(object = "odm1_1"),
		function(object,  ID=NULL, Organization=NULL, Description=NULL, Citation=NULL, ...){
		  	w.o <- list(where.clause = "", the.and  = "")
			w.o <- expand.where(w.o, ID, "SourceID", exact=TRUE)
			w.o <- expand.where(w.o, Organization, "Organization", exact=FALSE)
			w.o <- expand.where(w.o, Citation, "Citation", exact=FALSE)
			w.o <- expand.where(w.o, Description, "SourceDescription", exact=FALSE)
			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)

			query <- paste("SELECT * FROM Sources ", where.clause)
			res <- run.query(object, query )
			#ToDo: Define standard return values and order
			if(NCOL(res)>0) {names(res)[1:3] <- c("ID", "Organization", "Description")}
			return(res)
		}
	)
setMethod("IgetDeletedDataValues", 
		signature(object = "odm1_1Ver"),
		function(object,ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, VersionID, exact=FALSE, ...  ){

			where.clause <- assembleDataWhereClause(ID=ID, from=from, to=to, tz=tz, SiteID=SiteID, VariableID=VariableID, Offset=Offset, OffsetTypeID=OffsetTypeID, CensorCode=CensorCode, QualifierID=QualifierID, MethodID=MethodID, SourceID=SourceID, SampleID=SampleID, DerivedFromID=DerivedFromID, QualityControlLevelID=QualityControlLevelID, exact=exact)
			if(where.clause==""){
				where.clause <- "WHERE"
			} else {
				where.clause <- paste(where.clause, "AND")
			}

			the.query <- paste("SELECT * FROM DataValuesRepository NATURAL JOIN (SELECT ValueID, max(VersionID) AS VersionID FROM DataValuesRepository R ", where.clause, " NOT EXISTS (SELECT * FROM DataValues N WHERE N.ValueID = R.ValueID ) GROUP BY ValueID) AS RecSelect ")
			to.ret <- run.query(object, the.query )
			return(to.ret)
		}
	)

setMethod("IgetOldDataValues", 
		signature(object = "odm1_1Ver"),
		function(object,ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, VersionID, exact=FALSE, ...  ){

			where.clause <- assembleDataWhereClause(ID=ID, from=from, to=to, tz=tz, SiteID=SiteID, VariableID=VariableID, Offset=Offset, OffsetTypeID=OffsetTypeID, CensorCode=CensorCode, QualifierID=QualifierID, MethodID=MethodID, SourceID=SourceID, SampleID=SampleID, DerivedFromID=DerivedFromID, QualityControlLevelID=QualityControlLevelID, exact=exact)
			if(where.clause==""){
				where.clause <- "WHERE"
			} else {
				where.clause <- paste(where.clause, "AND")
			}

			the.query <- paste("SELECT * FROM DataValuesRepository NATURAL JOIN (SELECT ValueID, min(VersionID) AS VersionID FROM DataValuesRepository ", where.clause, " VersionID >= ",VersionID," GROUP BY ValueID) AS RecSelect ")
			to.ret <- run.query(object, the.query )
			return(to.ret)
		}
	)
setMethod("IgetDeletedDataValues", 
		signature(object = "odm1_1"),
		function(object,ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, VersionID, ...  ){

			#Do nothing because no version
		}
)

setMethod("IgetOldDataValues", 
		signature(object = "odm1_1"),
		function(object,ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, VersionID, ...  ){

			#Do nothing because no version
		}
)
setMethod("IgetDataValues", 
		signature(object = "odm1_1"),
		function(object,ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, exact=FALSE, ...  ){

			where.clause <- assembleDataWhereClause(ID=ID, from=from, to=to, tz=tz, SiteID=SiteID, VariableID=VariableID, Offset=Offset, OffsetTypeID=OffsetTypeID, CensorCode=CensorCode, QualifierID=QualifierID, MethodID=MethodID, SourceID=SourceID, SampleID=SampleID, DerivedFromID=DerivedFromID, QualityControlLevelID=QualityControlLevelID, exact=exact)

			query <- paste("SELECT * FROM DataValues ", where.clause)
			to.ret <- run.query(object, query )
			#ToDo: Define standard return values and order
			return(to.ret)
		}
	)
setMethod("IdeleteDataValues", 
		signature(object = "odm1_1"),
		function(object, ValueID){
			the.query <- paste("DELETE FROM DataValues 
					 WHERE ValueID = ", paste(ValueID, collapse=" OR ValueID = "), sep="")
				run.query(object, the.query)

		}
)
setMethod("IarchiveDataValues", 
		signature(object = "odm1_1"),
		function(object, ValueID,reason){
			#Do nothing as this does not have a version management system
		}
)
setMethod("IgetDataVersions", 
		signature(object = "odm1_1"),
		function(object){
			stop("Versions not implemented")
		}
)
setMethod("IgetCurrentDataVersion", 
		signature(object = "odm1_1"),
		function(object){
			stop("Versions not implemented")
		}
)
setMethod("IaddDataVersion", 
		signature(object = "odm1_1"),
		function(object, reason){
			stop("Versions not implemented")
		}
)

setMethod("IgetDataVersions", 
		signature(object = "odm1_1Ver"),
		function(object){
			#Do nothing as this does not have a version management system
			the.query <- "SELECT * FROM Versions"
			to.ret <- run.query(object, the.query )
			return(to.ret)
		}
)

setMethod("IgetCurrentDataVersion", 
		signature(object = "odm1_1Ver"),
		function(object){
			#Do nothing as this does not have a version management system
			the.query <- "SELECT MAX( VersionID ) FROM Versions"
			to.ret <- run.query(object, the.query )
			return(as.numeric(to.ret))
		}
)

setMethod("IaddDataVersion", 
		signature(object = "odm1_1Ver"),
		function(object, reason){
			# Eintrag in VersionsTabelle ValidUntilDate zum alten Datensatz
			the.query <- paste("UPDATE Versions SET ValidUntil = ", sqlstatements(object, "now") ," WHERE VersionID =", IgetCurrentDataVersion(object))
			run.query(object, the.query )
			# Neuer Eintrag in VersionsTabelle mit neuer Begruendung
			the.query <- paste("INSERT INTO Versions (VersionComment) Values (\"",reason,"\")", sep="")
			run.query(object, the.query )

			# Achtung: beim Abholen von der Versionshistory 
			#	muss die Begruendung for die Aenderung jeweils 
			#	bei der naechsten Versionsnummer nachgeschaut werden.
			#the.query <- paste("SELECT ", sqlstatements(object, "last_id"))
			the.query <- paste("SELECT MAX(VersionID) FROM Versions")
			to.ret <- run.query(object, the.query )
			return(as.numeric(to.ret))
		}
)
setMethod("IarchiveDataValues", 
		signature(object = "odm1_1Ver"),
		function(object, ValueID,reason){
			stopifnot(!is.null(reason))
			stopifnot(length(ValueID)>0)
			# daten kopieren mit aktueller Datenversion als ValidUntilID
			currentVersion <-  IgetCurrentDataVersion(object)
			if(is.na(currentVersion)){
				the.query <- paste("INSERT INTO Versions (VersionComment) Values ('Initial Version')", sep="")
				run.query(object, the.query )
				currentVersion <-  IgetCurrentDataVersion(object)
			}
			the.query = paste("INSERT INTO DataValuesRepository SELECT ValueID, VersionID, DataValue, ValueAccuracy, LocalDateTime, UTCOffset, DateTimeUTC, SiteID, VariableID, OffsetValue, OffsetTypeID, CensorCode, QualifierID, MethodID, SourceID, SampleID, DerivedFromID, QualityControlLevelID FROM DataValues, (SELECT ",currentVersion," as VersionID) as Version WHERE ValueID = ", paste(ValueID, collapse=" OR ValueID = "), sep="")
			run.query(object, the.query )
			# Version updaten
			IaddDataVersion(object, reason)

		}
)
setMethod("IupdateDataValues", 
		signature(object = "odm1_1"),
		function(object, ValueID, localDateTime, value, TZ, SiteID, VariableID, Offset=NULL, OffsetTypeID=NULL, CensorCode, QualifierID=NULL, MethodID, SourceID, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID, valueAccuracy=NULL,...  ){
				insert.query <- paste("UPDATE DataValues SET DataValue = \"",value,
					"\", ValueAccuracy = ", valueAccuracy,
				       	", LocalDateTime = \"", localDateTime,
				       	"\", UTCOffset = \"", tz2offset(TZ),
					"\", DateTimeUTC = \"", strftime(localDateTime, tz=TZ),
					"\", SiteID = ", SiteID, 
					", VariableID = ", VariableID,
					", OffsetValue = ", Offset,
					", OffsetTypeID = ", OffsetTypeID,
				        ", CensorCode = \"", CensorCode,
					"\", QualifierID = ", QualifierID,
				        ", MethodID = ", MethodID, 
					", SourceID = ", SourceID,
					", SampleID = ", SampleID, 
					", DerivedFromID = ", DerivedFromID,
					", QualityControlLevelID = ", QualityControlLevelID,
					" WHERE ValueID = ", ValueID, sep="")
				run.query(object, insert.query)

		}
)
setMethod("IaddVariable", signature=(object = "odm1_1"),
	  function(object, Code, Name, Speciation, Unit, SampleMedium,ValueType, IsRegular, TimeSupport, TimeUnits, DataType, GeneralCategory, NoDataValue) {
			for(rownum in seq(along=Name)){
				#no Foreign Key
		  		theCode <- sv(Code, rownum)
		  		theNoDataValue <- sv(NoDataValue, rownum)
				theIsRegular <- sv(IsRegular, rownum)
				theTimeSupport <- sv(TimeSupport, rownum)
				#Handle CV-Fields: VariableName, Speciation, SampleMedium,
			        #	ValueType, DataType, GeneralCategory
				#Other Fields with foreign key: VariableUnits, TimeUnits
				theName <- svk(Name, "Name", rownum, object)
				theSpeciation <- svk(Speciation, "Speciation", rownum, object)
				theUnit <- svk(Unit, "Unit", rownum, object)
				theSampleMedium <- svk(SampleMedium, "SampleMedium", rownum, object)
				theValueType <- svk(ValueType, "ValueType", rownum, object)
				theTimeUnits <- svk(TimeUnits, "TimeUnits", rownum, object)
				theDataType <- svk(DataType, "DataType", rownum, object)
				theGeneralCategory <- svk(GeneralCategory, "GeneralCategory", rownum, object)


				insert.query <- paste("INSERT INTO Variables (VariableCode, VariableName, Speciation, VariableUnitsID, SampleMedium,ValueType, IsRegular, TimeSupport, TimeUnitsID, DataType, GeneralCategory, NoDataValue) VALUES (\"", paste(
						theCode,
						theName ,
						theSpeciation ,
						theUnit ,
						theSampleMedium ,
						theValueType ,
						theIsRegular ,
						theTimeSupport ,
						theTimeUnits ,
						theDataType ,
						theGeneralCategory ,
						sep="\", \""), "\",", theNoDataValue  ," )", sep="")
			 	run.query(object, insert.query )
			}
	  }
)
setMethod("IaddQualityControlLevel",
		signature(object = "odm1_1"),
		function(object,  ID, Code, Definition, Explanation) {
		for(rownum in seq(along=ID)){
			insert.query <- paste("INSERT INTO QualityControlLevels (QualityControlLevelID, QualityControlLevelCode, Definition, Explanation) VALUES (\"", paste(
					ID[rownum],
					Code[rownum],
					Definition[rownum],
					Explanation[rownum],
					sep="\", \""), "\")", sep="")
			run.query(object, insert.query )

		}
	}
)
setMethod("IaddISOMetadata",
		signature(object = "odm1_1"),
		function(object, TopicCategory, Title, Abstract, ProfileVersion, MetadataLink) {
		for(rownum in seq(along=Title)){
			theTopicCategory <- svk(TopicCategory, "TopicCatgeoryCV", rownum,object)
			insert.query <- paste("INSERT INTO ISOMetadata (TopicCategory, Title, Abstract, ProfileVersion, MetadataLink) VALUES (\"", paste(
					theTopicCategory,
					Title[rownum],
					Abstract[rownum],
					ProfileVersion[rownum],
					MetadataLink[rownum],
					sep="\", \""), "\")", sep="")
			run.query(object, insert.query )

		}
	}
)
setMethod("IaddSource",
		signature(object = "odm1_1"),
		function(object, Organization, SourceDescription, SourceLink, ContactName, Phone, Email, Address, City, State, ZipCode, Citation, Metadata) {
		for(rownum in seq(along=SourceDescription)){
			theMetadata <- svk(Metadata, "ISOMetadata", rownum,object)
			insert.query <- paste("INSERT INTO Sources (Organization, SourceDescription, SourceLink, ContactName, Phone, Email, Address, City, State, ZipCode, Citation, MetadataID) VALUES (\"", paste(
					Organization[rownum],
					SourceDescription[rownum],
					SourceLink[rownum],
					ContactName[rownum],
					Phone[rownum],
					Email[rownum],
					Address[rownum],
					City[rownum],
					State[rownum],
					ZipCode[rownum],
					Citation[rownum],
					theMetadata,
					sep="\", \""), "\")", sep="")
			run.query(object, insert.query )

		}
	}
)

setMethod("IaddSite", 
		signature(object = "odm1_1"),
		function(object, Code, Name, Latitude, Longitude, LatLongDatum, Elevation=NULL, VerticalDatum=NULL, LocalX=0,LocalY=0, LocalProjection=NULL, PosAccuracy=0, State=NULL, County=NULL, Comments=NULL){
			for(rownum in seq(along=Name)){

				#no Foreign Key
				#replace with
				# :'<,'>s/\i*$/the& <- sv(&, rownum)/  
				theElevation <- sv(Elevation, rownum)
				theLocalX <- sv(LocalX, rownum)
				theLocalY <- sv(LocalY, rownum)
				thePosAccuracy <- sv(PosAccuracy, rownum)
				theState <- sv(State, rownum)
				theCounty <- sv(County, rownum)
				theComment <- sv(Comments, rownum)
				#with Foreign Key 
				# :'<,'>s/\i*$/the& <- svk(&, "&", rownum, object)/   
				theVerticalDatum <- svk(VerticalDatum, "VerticalDatumCV",rownum,object)
				theLocalProjection <- svk(LocalProjection, "SpatialReference", rownum,object)

				insert.query <- paste("INSERT INTO Sites (SiteCode, SiteName, Latitude, Longitude, LatLongDatumID, Elevation_m, VerticalDatum, LocalX, LocalY, LocalProjectionID, PosAccuracy_m, State, County, Comments) VALUES (", paste(

						paste('"',Code[rownum], '"',sep=""),
						paste('"',Name[rownum], '"',sep=""),
						Latitude[rownum],
						Longitude[rownum],
						LatLongDatum[rownum],
						theElevation,
						paste('"',theVerticalDatum, '"',sep=""),
						theLocalX,
						theLocalY,
						theLocalProjection,
						thePosAccuracy,
						paste('"',theState, '"',sep=""),
						paste('"',theCounty, '"',sep=""),
						paste('"',theComment, '"',sep=""),
						sep=", "), ")", sep="")
			 	run.query(object, insert.query )
			}

		}
)
setMethod("IaddDataValues", 
		signature(object = "odm1_1"),
		function(object, localDateTime, values, TZ, SiteID, VariableID, Offset=NULL, OffsetTypeID=NULL, CensorCode, QualifierID=NULL, MethodID, SourceID, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID, valueAccuracy=NULL,...  ){
			for(rownum in seq(along=values)){
				
				#with Foreign Key 
				theOffsetTypeID <- svk(OffsetTypeID, "OffsetType", rownum,object)
				theQualifierID <- svk(QualifierID, "Qualifier", rownum,object)
				theSampleID <- svk(SampleID, "Sample", rownum,object)
				theSiteID <- svk(SiteID, "SiteID", rownum, object)
				theVariableID <- svk(VariableID, "VariableID", rownum, object)
				theMethodID <- svk(MethodID, "MethodID", rownum, object)
				theCensorCode <- svk(CensorCode, "CensorCodeCV", rownum, object)
				theSourceID <- svk(SourceID, "SourceID", rownum, object)
				theQualityControlLevelID <- svk(QualityControlLevelID, "QualityControlLevelID", rownum, object)

				#without Foreign Key
				theOffset <- sv(Offset, rownum)
				thevalue <- sv(values, rownum)
				if(thevalue=="NULL") next # skip NA values from importing
				thevalueAccuracy <- sv(valueAccuracy, rownum)
				thelocalDateTime <- sv(localDateTime, rownum)
				theTZ <- tz2offset(sv(TZ, rownum))
				theDerivedFromID <- sv(DerivedFromID, rownum)

				insert.query <- paste("INSERT INTO DataValues (DataValue, ValueAccuracy, LocalDateTime, UTCOffset, DateTimeUTC, SiteID, VariableID, OffsetValue, OffsetTypeID, CensorCode, QualifierID, MethodID, SourceID, SampleID, DerivedFromID, QualityControlLevelID) VALUES (", paste(thevalue,
						thevalueAccuracy,
						paste("'", strftime(thelocalDateTime, "%Y-%m-%d %H:%M:%S", tz=as.character(theTZ) ), "'", sep=""),
						theTZ , 
						paste("'", strftime(thelocalDateTime, tz="GMT", ), "'", sep=""), 
						theSiteID, theVariableID,
						theOffset,
						theOffsetTypeID,
						paste("'", theCensorCode, "'", sep=""),
						theQualifierID,
						theMethodID,
						theSourceID,
						theSampleID,
						theDerivedFromID,
						theQualityControlLevelID ,sep=", "), ")", sep="")
			 	run.query(object, insert.query )
			}
		}
	)
setMethod("IgetNo",
	  signature(object= "odm1_1"),
	  function(object, table){
		#lookup table for database information
		tab.def <- matrix(
				 #db table, db attribute
				c("SpatialReferences", "SpatialReferences", "SRSName", "SRSID, IsGeographic, Notes, SpatialReferenceID","0, 0, '', 'No'", "SpatialReferenceID",
	  			  "OffsetType", "OffsetTypes", "OffsetDescription", "OffsetUnitsID", "1", "OffsetTypeID",
				  "Qualifier", "Qualifiers", "QualifierCode", "QualifierDescription", "'No Qualifier - to enable optional fields with foreign keys'", "QualifierID" ,
				  "LabMethod", "LabMethods", "LabMethodName", "LabMethodDescription, LabName, LabOrganization, LabMethodLink", "'No LabMethod - to enable optional fields with foreign keys', 'NoLab', 'NoOrganization', 'NoLink'", "LabMethodID" ,

				  #does not work because samples depends also on LabMethods
				  "Sample", "Samples", "LabSampleCode", "", "", "SampleID",
				  "Method", "Methods", "MethodDescription", "MethodLink", "''", "MethodID",
				  "QualityControlLevel", "QualityControlLevel", "QualityControlLevelCode", "Definition, Explanation, QualityControlLevelID", "'No Code', 'Default entry if no code is available', 999", "QualityControlLevelID",
				  "ISOMetadata", "ISOMetadata", "Title", "TopicCategory, Abstract", "'Unknown', 'Entry with unknown ISOMetadata'", "MetadataID"
				),
		  	byrow=TRUE, ncol=6)
		CVtab <- paste(CVtables(), "CV", sep="")
		lc <- length(CVtab)
		CVentries <- cbind(CVtab, CVtab, rep("Definition", lc), rep("Term",lc), rep("'no'",lc), rep("Term", lc))
		tab.def <- rbind(tab.def, CVentries)
		tab.def <- as.data.frame(tab.def, stringsAsFactors=FALSE)
		names(tab.def) <- c("entry", "tab", "col", "other.fields", "default.values", "primary.key")
		tab.def.bak <- tab.def
		tab.def <- tab.def[tab.def$entry %in% c(table, paste(table,"CV",sep="")),]
		if(NROW(tab.def)!=1){
			cat("Existing tables:\n")
			print(tab.def.bak$tab)
		       	stop(paste("IgetNo: No definitions for table", table))
		}

		#Foreign keys
		#generate?
		fk.list <- list("Sample" = c("LabMethod", "SampleType"))

		query <- paste("SELECT * FROM ", tab.def$tab ," WHERE ", tab.def$col,"='No",tab.def$tab,"'", sep="")
		res <- run.query(object, query)
		if(NROW(res)==0){
			fks <- fk.list[[table]]
			if(!is.null(fks)){
				for(thefk in fks){
					novalue <- IgetNo(object, thefk)
					if(thefk %in% CVtables()){
						fieldname <- thefk
					} else {
						fieldname <- paste(thefk, "ID", sep="")
					}
					if(tab.def$other.fields!=""){
					       	thesep = ", "
					} else {
						thesep = ""
					}
					tab.def$other.fields <- paste(tab.def$other.fields, thesep, fieldname, sep="")
					tab.def$default.values <- paste(tab.def$default.values, thesep, "'", novalue, "'", sep="")


				}
			}
			ins.query <- paste("INSERT INTO ", tab.def$tab ," (",tab.def$col, ", ", tab.def$other.fields,") Values ('No",tab.def$tab,"',", tab.def$default.values,")", sep="" )
			res <- run.query(object, ins.query)
			res <- run.query(object, query)

		}

		#necessary because postres does not support capital letters
		coln <- which(tolower(tab.def$primary.key) %in% tolower(names(res)))

		toret <-res[[coln]] 
		stopifnot(!is.null(toret))
		return(toret)
	  }
)

setMethod("IgetCV",
	  signature(object= "odm1_1"),
	  function(object, table, term, definition, exact=FALSE){
		  #check for valid tables
		  stopifnot(table %in% CVtables())
		        definition <- substr(definition, 1, 255)
		  	where.object <- list(where.clause = "", the.and  = "")

			where.object <- expand.where(where.object, term, "Term", exact=exact)
			where.object <- expand.where(where.object, definition, "Definition", exact=exact)

			where.clause <- where.object$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)

			query <- paste("SELECT * FROM ",table,"CV", where.clause, sep="")
			res <- run.query(object, query )
			return(res)
	  }
)
setMethod("IaddCV",
	  signature(object= "odm1_1"),
	  function(object, table, term, definition){
		  #check for valid tables
		  stopifnot(table %in% CVtables())
		        definition <- substr(definition, 1, 255)
			query <- paste("INSERT INTO ",table,"CV (Term, Definition) Values (\"", term, "\",\"",definition,"\")"  , sep="")
			res <- run.query(object, query )
			return(res)
	  }
)
setMethod("IaddUnits",
	  signature(object= "odm1_1"),
	  function(object, Name, Type, Abbreviation){
			query <- paste("INSERT INTO Units (UnitsName, UnitsType, UnitsAbbreviation) Values (\"",Name, "\",\"", Type ,"\",\"",Abbreviation,"\")"  , sep="")
			res <- run.query(object, query )
			return(res)
	  }
)

setMethod("IaddSpatialReference",
	  signature(object= "odm1_1"),
	  function(object, ID, SRSName, SRSID, IsGeographic, Notes){
		  #check for valid tables
			query <- paste("INSERT INTO SpatialReferences (SpatialReferenceID, SRSID, SRSName, IsGeographic, Notes) Values (\"", ID, "\",\"", SRSID, "\",\"",SRSName, "\",\"", IsGeographic ,"\",\"",Notes,"\")"  , sep="")
			res <- run.query(object, query )
			return(res)
	  }
)

setMethod("IaddSynonym",
	       signature(object= "odm1_1"),
       	       function(object, phrase, table, id){
		       query <- paste('INSERT INTO Synonyms (phrase, RecID, tab) Values ("',
				       phrase, '","', id, '","', table,'")', sep='')
			res <- run.query(object, query)
	       })
setMethod("IgetSynonymID",
	       signature(object= "odm1_1"),
	       function(object, phrase, table){
		       query <- paste('SELECT RecID FROM Synonyms WHERE phrase = "',
				       phrase, '" AND tab="', table, '"', sep='')
			res <- run.query(object, query)
			if(NCOL(res)==0){
				return(NULL)
			} else {
				return(res[,1])
			}

	       })
setMethod("IgetSynonyms",
	       signature(object= "odm1_1"),
	       function(object){
		       query <- paste('SELECT * FROM Synonyms ')
			res <- run.query(object, query)
			if(NCOL(res)==0){
				return(NULL)
			} else {
				return(res)
			}

	       })
