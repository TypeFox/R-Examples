assembleDataWhereClause <- function(ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), SiteID=NULL, VariableID=NULL, Offset=NULL, OffsetTypeID=NULL, CensorCode=NULL, QualifierID=NULL, MethodID=NULL, SourceID=NULL, SampleID=NULL, DerivedFromID=NULL, QualityControlLevelID=NULL, exact=FALSE){

			tz <- match.arg(tz)
			#the where.object is used to assemble the where clause
		  	w.o <- list(where.clause = "", the.and  = "")

			#w.o <- expand.where(w.o, term, "Term", exact=exact)

			#assembling standard terms
			w.o <- expand.where(w.o, ID, "ValueID", exact=TRUE)
			w.o <- expand.where(w.o, SiteID, "SiteID", exact=TRUE)
			w.o <- expand.where(w.o, QualifierID, "QualifierID", exact=TRUE)
			w.o <- expand.where(w.o, VariableID, "VariableID", exact=TRUE)
			w.o <- expand.where(w.o, CensorCode, "CensorCode", exact=TRUE)
			w.o <- expand.where(w.o, MethodID, "MethodID", exact=TRUE)
			w.o <- expand.where(w.o, SourceID, "SourceID", exact=TRUE)
			w.o <- expand.where(w.o, SampleID, "SampleID", exact=TRUE)
			w.o <- expand.where(w.o, DerivedFromID, "DerivedFromID", exact=TRUE)
			w.o <- expand.where(w.o, QualityControlLevelID, "QualityControlLevelID", exact=TRUE)

			
			if(tz=="local"){
				datefield <- "LocalDateTime"
			} else {
				datefield <- "DateTimeUTC"
			}

			if(!is.null(from)){
				w.o$where.clause <- paste(w.o$where.clause, w.o$the.and, datefield, " >= '", strftime(from, tz="GMT"), "'", sep="")
				w.o$the.and <- " AND "
			}
			if(!is.null(to)){
				w.o$where.clause <- paste(w.o$where.clause, w.o$the.and, datefield, " <= '", strftime(to,tz="GMT"), "'", sep="")
				w.o$the.and <- " AND "
			}

			if(!is.null(Offset)  &  !is.null(OffsetTypeID)){
				if(!is.na(Offset) &  !is.na(OffsetTypeID)){
					w.o$where.clause <- paste(w.o$where.clause, w.o$the.and, "OffsetValue like '%", Offset, "%' AND OffsetTypeID = ", OffsetTypeID, sep="")
					w.o$the.and <- " AND "
				}
			}

			where.clause <- w.o$where.clause
			if(where.clause!="") where.clause <- paste(" WHERE ", where.clause)

			return(where.clause)
}
