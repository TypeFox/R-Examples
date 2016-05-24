`convertScaleScore` <-
function(tmp.data.for.equate,
	tmp.year.for.equate,
	equate.list,
	conversion.type,
	equating.method,
	state) {

	SCALE_SCORE <- SCALE_SCORE_EQUATED <- GRADE <- CONTENT_AREA <- YEAR <- NULL

	### Define variables

	tmp.unique.years <- unique(tmp.data.for.equate$YEAR)

	if (conversion.type=="NEW_TO_OLD") {
		tmp.years.for.equate <- tmp.unique.years[tmp.unique.years >= tmp.year.for.equate]
	} else {
		tmp.years.for.equate <- tmp.unique.years[tmp.unique.years < tmp.year.for.equate]
	}
	if (paste("SCALE_SCORE_EQUATED", equating.method, conversion.type, sep="_") %in% names(tmp.data.for.equate)) tmp.data.for.equate[,paste("SCALE_SCORE_EQUATED", equating.method, conversion.type, sep="_"):=NULL, with=FALSE]


	### Create scale.score.concordance lookup

	for (i.iter in names(equate.list)) {
		i <- unlist(strsplit(i.iter, "[.]"))[1]
		for (j.iter in names(equate.list[[i.iter]])) {
			j <- unlist(strsplit(j.iter, "_"))[2]
			tmp.data.for.equate[YEAR %in% tmp.years.for.equate & CONTENT_AREA==i & GRADE==j,
				paste("SCALE_SCORE_EQUATED", toupper(equating.method), conversion.type, sep="_"):=equate.list[[paste(i, tmp.year.for.equate, sep=".")]][[paste("GRADE", j, sep="_")]][[toupper(equating.method)]][[conversion.type]][['interpolated_function']](SCALE_SCORE), with=FALSE]
		}
	}
	tmp.data.for.equate[!YEAR %in% tmp.years.for.equate, paste("SCALE_SCORE_EQUATED", toupper(equating.method), conversion.type, sep="_") := SCALE_SCORE, with=FALSE]

	return(tmp.data.for.equate)
} ### END convertScaleScore function
