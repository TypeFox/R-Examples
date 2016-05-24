`getYearsContentAreasGrades` <- 
function(state,
	years,
	content_areas=NULL,
	content_areas_domains=NULL,
	earliest_year_reported=NULL) {

	CONTENT_AREA <- NULL

	tmp.list <- list()
	if (is.null(content_areas)) tmp.content_areas <- content_areas_domains else tmp.content_areas <- content_areas
	for (i in intersect(names(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]]), tmp.content_areas)) {
		if (!is.null(content_areas_domains) && !is.null(SGP::SGPstateData[[state]][["Student_Report_Information"]][['Content_Areas_Domains']])) {
			tmp.dt <- data.table(GRADE=as.character(unique(unlist(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][
					grep(i, SGP::SGPstateData[[state]][["Student_Report_Information"]][['Content_Areas_Domains']])]))))
		} else {
			tmp.dt <- data.table(GRADE=as.character(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[i]]))
		}

		if (!is.null(earliest_year_reported[[i]])) {
			tmp.years.diff <- as.numeric(tail(unlist(strsplit(tail(sort(years), 1), "_")), 1))-as.numeric(tail(unlist(strsplit(earliest_year_reported[[i]], "_")), 1))
			tmp.dt <- CJ(tmp.dt$GRADE, intersect(yearIncrement(tail(sort(years), 1), c(-seq(tmp.years.diff), 0)), years))
		} else {
			tmp.dt <- CJ(tmp.dt$GRADE, years)
		}

		setnames(tmp.dt, c("GRADE", "YEAR")) 
		tmp.list[[i]] <- data.table(CONTENT_AREA=i, tmp.dt)
	}
	tmp.dt <- rbindlist(tmp.list, fill=TRUE)
	setkeyv(tmp.dt, c("CONTENT_AREA", "GRADE", "YEAR"))
	return(tmp.dt[!is.na(CONTENT_AREA)])
} ## END getYearsContentAreasGrades
