`getPanelData` <- 
function(sgp.data,
	sgp.type,
	sgp.iter,
	sgp.csem=NULL,
	sgp.scale.score.equated=NULL,
	sgp.targets=NULL,
	SGPt=NULL) {

	YEAR <- CONTENT_AREA <- VALID_CASE <- V3 <- V5 <- ID <- GRADE <- SCALE_SCORE <- YEAR_WITHIN <- tmp.timevar <- FIRST_OBSERVATION <- LAST_OBSERVATION <- ACHIEVEMENT_LEVEL <- DATE <- NULL

	if (is(sgp.data, "DBIObject")) {
		con <- dbConnect(SQLite(), dbname = "Data/tmp_data/TMP_SGP_Data.sqlite")
		var.names <- dbListFields(con, "sgp_data")
		sqlite.tf <- TRUE
	} else {
		var.names <- names(sgp.data)
		sqlite.tf <- FALSE
	}

	if ("STATE" %in% var.names) state <- "STATE" else state <- NULL


	###
	### sgp.percentiles
	###
	
	if (sgp.type=="sgp.percentiles") {

		if (!is.null(tmp.exclude.lookup <- sgp.iter$sgp.exclude.sequences)) {
			if (is.data.table(tmp.exclude.lookup)) {
				tmp.exclude.lookup <- setkey(data.table(tmp.exclude.lookup[, list(VALID_CASE, CONTENT_AREA, YEAR, GRADE)]))
			} else stop("Element 'sgp.exclude.sequences' of sgp.config must be a data table with variables 'VALID_CASE', 'CONTENT_AREA', 'YEAR', and 'GRADE'.")
			if (sqlite.tf) {
				tmp.exclude.ids <- NULL
				for (y in 1:nrow(tmp.exclude.lookup)) {
					tmp.exclude.ids <- c(tmp.exclude.ids, unlist(dbGetQuery(con, paste(
						"select ID from sgp_data where CONTENT_AREA in ('", paste(tmp.exclude.lookup[y]$CONTENT_AREA, collapse="', '"), "')",
												 " AND GRADE in ('", paste(tmp.exclude.lookup[y]$GRADE, collapse="', '"), "')",
												 " AND YEAR in ('", paste(unique(tmp.exclude.lookup[y]$YEAR), collapse="', '"), "')", sep=""))))
				}
				tmp.exclude.ids <- unique(tmp.exclude.ids)
			} else tmp.exclude.ids <- unique(sgp.data[tmp.exclude.lookup][['ID']])
		} else tmp.exclude.ids <- as.character(NA)
	
		if ("YEAR_WITHIN" %in% var.names) {
			tmp.lookup <- data.table(V1="VALID_CASE", tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])),
				tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])), sgp.iter[["sgp.grade.sequences"]],
				tail(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.grade.sequences"]])), FIRST_OBSERVATION=as.integer(NA), LAST_OBSERVATION=as.integer(NA))
			tmp.lookup[grep("FIRST", V5, ignore.case=TRUE), FIRST_OBSERVATION:=1L]; tmp.lookup[grep("LAST", V5, ignore.case=TRUE), LAST_OBSERVATION:=1L]; tmp.lookup[,V5:=NULL]
			setnames(tmp.lookup, paste("V", 1:4, sep=""), c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))
			
			tmp.lookup.list <- list()
			for (i in unique(sgp.iter[["sgp.panel.years.within"]])) {
				if (sqlite.tf) {
					setkeyv(tmp.lookup, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					suppressWarnings(tmp.lookup.list[[i]] <- data.table(dbGetQuery(con, 
						paste("select * from sgp_data where CONTENT_AREA in ('", paste(tmp.lookup[get(i)==1]$CONTENT_AREA, collapse="', '"), "')",
						" AND GRADE in ('", paste(tmp.lookup[get(i)==1]$GRADE, collapse="', '"), "')", " AND ", i, " in (1)", 
						" AND YEAR in ('", paste(tmp.lookup[get(i)==1]$YEAR, collapse="', '"), "')", sep="")), 
						key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))[tmp.lookup[get(i)==1], nomatch=0][!ID %in% tmp.exclude.ids][,
							'tmp.timevar':=paste(YEAR, CONTENT_AREA, i, sep="."), with=FALSE][, 
							c("ID", "GRADE", "SCALE_SCORE", "YEAR_WITHIN", "tmp.timevar", sgp.csem, sgp.scale.score.equated), with=FALSE])
				} else {
					setkeyv(sgp.data, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					setkeyv(tmp.lookup, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					suppressWarnings(tmp.lookup.list[[i]] <- data.table(sgp.data[tmp.lookup[get(i)==1], nomatch=0][!ID %in% tmp.exclude.ids][,
							'tmp.timevar':=paste(YEAR, CONTENT_AREA, i, sep="."), with=FALSE][,
						c("ID", "GRADE", "SCALE_SCORE", "YEAR_WITHIN", "tmp.timevar", sgp.csem, sgp.scale.score.equated), with=FALSE], key="ID")) ### Could be NULL and result in a warning
				}
			}
			if (sqlite.tf) dbDisconnect(con)
			
			if (tail(sgp.iter[['sgp.panel.years']], 1)==head(tail(sgp.iter[['sgp.panel.years']], 2), 1)) {
				setkey(tmp.lookup.list[[1]], ID); setkey(tmp.lookup.list[[2]], ID)
				tmp.ids <- intersect(tmp.lookup.list[[1]][['ID']], tmp.lookup.list[[2]][['ID']])
				tmp.ids <- tmp.ids[tmp.lookup.list[[1]][tmp.ids][['YEAR_WITHIN']] < tmp.lookup.list[[2]][tmp.ids][['YEAR_WITHIN']]]
				tmp.lookup.list <- lapply(tmp.lookup.list, function(x) x[tmp.ids])
			}
			return(ddcast(rbindlist(tmp.lookup.list),
					ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", "YEAR_WITHIN", sgp.csem, sgp.scale.score.equated), sep="."))
		} else {
			tmp.lookup <- SJ("VALID_CASE", tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])),
				tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])), sgp.iter[["sgp.grade.sequences"]])
			# ensure lookup table is ordered by years.  NULL out key after sorted so that it doesn't corrupt the join in ddcast.
			setkey(tmp.lookup, V3)
			setkey(tmp.lookup, NULL)

			if (sqlite.tf) {
				tmp.data <- data.table(dbGetQuery(con, 
					paste("select * from sgp_data where CONTENT_AREA in ('", paste(tmp.lookup$V2, collapse="', '"), "')", 
					" AND GRADE in ('", paste(tmp.lookup$V4, collapse="', '"), "')",
					" AND YEAR in ('", paste(tmp.lookup$V3, collapse="', '"), "')", sep="")), key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))
				dbDisconnect(con)
				if (is.null(SGPt)) {
					return(ddcast(tmp.data[tmp.lookup, nomatch=0][!ID %in% tmp.exclude.ids][,'tmp.timevar':=paste(YEAR, CONTENT_AREA, sep="."), with=FALSE],
							ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", sgp.csem, sgp.scale.score.equated), sep="."))
				} else {
					return(ddcast(tmp.data[tmp.lookup, nomatch=0][!ID %in% tmp.exclude.ids][,'tmp.timevar':=paste(YEAR, CONTENT_AREA, sep="."), with=FALSE][,
							c("TIME", "TIME_LAG"):=list(as.numeric(DATE), as.numeric(DATE-c(NA, DATE[-.N]))), by=ID],
							ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", sgp.csem, sgp.scale.score.equated, "TIME", "TIME_LAG"), sep="."))
				}
			} else {
				if (is.null(SGPt)) {
					return(ddcast(sgp.data[tmp.lookup, nomatch=0][!ID %in% tmp.exclude.ids][,'tmp.timevar':=paste(YEAR, CONTENT_AREA, sep="."), with=FALSE], 
							ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", sgp.csem, sgp.scale.score.equated), sep="."))
				} else {
					return(ddcast(sgp.data[tmp.lookup, nomatch=0][!ID %in% tmp.exclude.ids][,
							'tmp.timevar':=paste(YEAR, CONTENT_AREA, sep="."), with=FALSE][,
							c("TIME", "TIME_LAG"):=list(as.numeric(DATE), as.numeric(DATE-c(NA, DATE[-.N]))), by=ID],
							ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", sgp.csem, sgp.scale.score.equated, "TIME", "TIME_LAG"), sep="."))
				}
			}
		}
	} ### END if (sgp.type=="sgp.percentiles")


	###
	### sgp.projections & sgp.projections.baseline
	###
	
	if (sgp.type %in% c("sgp.projections", "sgp.projections.baseline")) {

		if (sgp.type=="sgp.projections") {
			sgp.projection.content.areas.label <- "sgp.projection.content.areas"
			sgp.projection.grade.sequences.label <- "sgp.projection.grade.sequences"
			sgp.projection.panel.years.label <- "sgp.projection.panel.years"
		} else {
			sgp.projection.content.areas.label <- "sgp.projection.baseline.content.areas"
			sgp.projection.grade.sequences.label <- "sgp.projection.baseline.grade.sequences"
			sgp.projection.panel.years.label <- "sgp.projection.baseline.panel.years"
		}

		if ("YEAR_WITHIN" %in% var.names) {
			tmp.lookup <- data.table(V1="VALID_CASE", tail(sgp.iter[[sgp.projection.content.areas.label]], length(sgp.iter[[sgp.projection.grade.sequences.label]])),
				sapply(head(sgp.iter[["sgp.panel.years"]], length(sgp.iter[[sgp.projection.grade.sequences.label]])), yearIncrement, tail(sgp.iter$sgp.panel.years.lags, 1)),
				sgp.iter[[sgp.projection.grade.sequences.label]], head(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[[sgp.projection.grade.sequences.label]])), 
				FIRST_OBSERVATION=as.integer(NA), LAST_OBSERVATION=as.integer(NA))
			tmp.lookup[grep("FIRST", V5, ignore.case=TRUE), FIRST_OBSERVATION:=1L]; tmp.lookup[grep("LAST", V5, ignore.case=TRUE), LAST_OBSERVATION:=1L]; tmp.lookup[,V5:=NULL]
			setnames(tmp.lookup, paste("V", 1:4, sep=""), c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))
			
			tmp.lookup.list <- list()
			for (i in unique(sgp.iter[["sgp.panel.years.within"]])) {
				if (sqlite.tf) {
					setkeyv(tmp.lookup, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					suppressWarnings(tmp.lookup.list[[i]] <- data.table(dbGetQuery(con, 
						paste("select * from sgp_data where CONTENT_AREA in ('", paste(tmp.lookup[get(i)==1]$CONTENT_AREA, collapse="', '"), "')", 
						" AND GRADE in ('", paste(tmp.lookup[get(i)==1]$GRADE, collapse="', '"), "')", " AND ", i, " in (1)",
						" AND YEAR in ('", paste(tmp.lookup[get(i)==1]$YEAR, collapse="', '"), "')", sep="")), key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))[tmp.lookup[get(i)==1], nomatch=0][,
						'tmp.timevar':=paste(YEAR, CONTENT_AREA, i, sep="."), with=FALSE][, 
							c("ID", "GRADE", "SCALE_SCORE", "YEAR_WITHIN", "tmp.timevar", sgp.scale.score.equated), with=FALSE])
				} else {
					setkeyv(sgp.data, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					setkeyv(tmp.lookup, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					suppressWarnings(tmp.lookup.list[[i]] <- data.table(sgp.data[tmp.lookup[get(i)==1], nomatch=0][,'tmp.timevar':=paste(YEAR, CONTENT_AREA, i, sep="."), with=FALSE][,
						c("ID", "GRADE", "SCALE_SCORE", "YEAR_WITHIN", "tmp.timevar", sgp.scale.score.equated), with=FALSE], key="ID")) ### Could be NULL and result in a warning
				}
			}
			if (sqlite.tf) dbDisconnect(con)

			if (tail(sgp.iter[['sgp.panel.years']], 1)==head(tail(sgp.iter[['sgp.panel.years']], 2), 1)) {
				tmp.ids <- intersect(tmp.lookup.list[[1]][['ID']], tmp.lookup.list[[2]][['ID']])
				tmp.ids <- tmp.ids[tmp.lookup.list[[1]][tmp.ids][['YEAR_WITHIN']] < tmp.lookup.list[[2]][tmp.ids][['YEAR_WITHIN']]]
				tmp.lookup.list <- lapply(tmp.lookup.list, function(x) x[tmp.ids])
			}
			if (is.null(sgp.targets)) {
				tmp.data <- ddcast(rbindlist(tmp.lookup.list),
						ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", "YEAR_WITHIN", state, sgp.scale.score.equated, SGPt), sep=".")
				setnames(tmp.data, tail(sort(grep("YEAR_WITHIN", names(tmp.data), value=TRUE)), 1), "YEAR_WITHIN")
				if (length(setdiff(grep("YEAR_WITHIN", names(tmp.data), value=TRUE), "YEAR_WITHIN")) > 0) {
					tmp.data[,setdiff(grep("YEAR_WITHIN", names(tmp.data), value=TRUE), "YEAR_WITHIN"):=NULL, with=FALSE]
				}
				if ("STATE" %in% var.names && dim(tmp.data)[1]!=0) {
					setnames(tmp.data, tail(sort(grep("STATE", names(tmp.data), value=TRUE)), 1), "STATE")
					if (length(setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE")) > 0) {
						tmp.data[,setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE"):=NULL, with=FALSE]
					}
				}
				return(tmp.data)
			} else {
				tmp.data <- ddcast(rbindlist(tmp.lookup.list),
					ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", "YEAR_WITHIN", state, sgp.scale.score.equated, SGPt), sep=".")[
					sgp.targets[CONTENT_AREA==tail(sgp.iter[[sgp.projection.content.areas.label]], 1) & YEAR==tail(sgp.iter[["sgp.panel.years"]], 1)], nomatch=0][,
					!c("CONTENT_AREA", "YEAR"), with=FALSE]
				setnames(tmp.data, tail(sort(grep("YEAR_WITHIN", names(tmp.data), value=TRUE)), 1), "YEAR_WITHIN")
				if (length(setdiff(grep("YEAR_WITHIN", names(tmp.data), value=TRUE), "YEAR_WITHIN")) > 0) {
					tmp.data[,setdiff(grep("YEAR_WITHIN", names(tmp.data), value=TRUE), "YEAR_WITHIN"):=NULL, with=FALSE]
				}
				if ("STATE" %in% var.names && dim(tmp.data)[1]!=0) {
					setnames(tmp.data, tail(sort(grep("STATE", names(tmp.data), value=TRUE)), 1), "STATE")
					if (length(setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE")) > 0) {
						tmp.data[,setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE"):=NULL, with=FALSE]
					}
				}
				return(tmp.data)
			}
		} ### END if ("YEAR_WITHIN" %in% var.names)

		tmp.lookup <- SJ("VALID_CASE", tail(sgp.iter[[sgp.projection.content.areas.label]], length(sgp.iter[[sgp.projection.grade.sequences.label]])),
			tail(sgp.iter[[sgp.projection.panel.years.label]], length(sgp.iter[[sgp.projection.grade.sequences.label]])), sgp.iter[[sgp.projection.grade.sequences.label]])
		# ensure lookup table is ordered by years.  NULL out key after sorted so that it doesn't corrupt the join in ddcast.
		setkey(tmp.lookup, V3)
		setkey(tmp.lookup, NULL)

		if (is.null(sgp.targets)) {
			if (sqlite.tf) {
				tmp.data <- ddcast(
					data.table(dbGetQuery(con, 
						paste("select * from sgp_data where CONTENT_AREA in ('", paste(tmp.lookup$V2, collapse="', '"), "')", 
						" AND GRADE in ('", paste(tmp.lookup$V4, collapse="', '"), "')", " AND YEAR in ('", paste(tmp.lookup$V3, collapse="', '"), "')", sep="")), 
						key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))[tmp.lookup, nomatch=0][,
						'tmp.timevar':=paste(YEAR, CONTENT_AREA, sep="."), with=FALSE],
					ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", state, sgp.scale.score.equated, SGPt), sep=".")
				dbDisconnect(con)
			} else {
				tmp.data <- ddcast(sgp.data[tmp.lookup, nomatch=0][,'tmp.timevar' := paste(YEAR, CONTENT_AREA, sep="."), with=FALSE],
					ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", state, sgp.scale.score.equated, SGPt), sep=".")
			}

			if ("STATE" %in% var.names && dim(tmp.data)[1]!=0) {
				setnames(tmp.data, tail(sort(grep("STATE", names(tmp.data), value=TRUE)), 1), "STATE")
				if (length(setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE")) > 0) tmp.data[,setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE"):=NULL, with=FALSE]
			}
			return(tmp.data)
		} else {
			if (sqlite.tf) {
				tmp.data <- ddcast(
					data.table(dbGetQuery(con, paste("select * from sgp_data where CONTENT_AREA in ('", paste(tmp.lookup$V2, collapse="', '"), "')", 
						" AND GRADE in ('", paste(tmp.lookup$V4, collapse="', '"), "')", " AND YEAR in ('", paste(tmp.lookup$V3, collapse="', '"), "')", sep=""))
						, key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))[tmp.lookup, nomatch=0][,'tmp.timevar':=paste(YEAR, CONTENT_AREA, sep="."), with=FALSE],
					ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", state, sgp.scale.score.equated, SGPt), sep=".")[
						sgp.targets[CONTENT_AREA==tail(sgp.iter[[sgp.projection.content.areas.label]], 1) & 
						YEAR==tail(sgp.iter[[sgp.projection.panel.years.label]], 1)], nomatch=0][,!c("CONTENT_AREA", "YEAR"), with=FALSE]
				dbDisconnect(con)
			} else {
				tmp.data <- ddcast(sgp.data[tmp.lookup, nomatch=0][, 'tmp.timevar' := paste(YEAR, CONTENT_AREA, sep="."), with=FALSE],
					ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", state, sgp.scale.score.equated, SGPt), sep=".")[
						sgp.targets[CONTENT_AREA==tail(sgp.iter[[sgp.projection.content.areas.label]], 1) & 
						YEAR==tail(sgp.iter[[sgp.projection.panel.years.label]], 1)], nomatch=0][,!c("CONTENT_AREA", "YEAR"), with=FALSE]
			}

			if ("STATE" %in% var.names && dim(tmp.data)[1]!=0) {
				setnames(tmp.data, tail(sort(grep("STATE", names(tmp.data), value=TRUE)), 1), "STATE")
				if (length(setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE")) > 0) tmp.data[,setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE"):=NULL, with=FALSE]
			}
			return(tmp.data)
		}
	} ### END if (sgp.type=="sgp.projections")


	###
	### sgp.projections.lagged
	###
	
	if (sgp.type %in% c("sgp.projections.lagged", "sgp.projections.lagged.baseline")) {

		if (sgp.type=="sgp.projections.lagged") {
			sgp.projection.content.areas.label <- "sgp.projection.content.areas"
			sgp.projection.grade.sequences.label <- "sgp.projection.grade.sequences"
			sgp.projection.panel.years.label <- "sgp.projection.panel.years"
		} else {
			sgp.projection.content.areas.label <- "sgp.projection.baseline.content.areas"
			sgp.projection.grade.sequences.label <- "sgp.projection.baseline.grade.sequences"
			sgp.projection.panel.years.label <- "sgp.projection.baseline.panel.years"
		}

		if ("YEAR_WITHIN" %in% var.names) {
			if (sqlite.tf) {
				tmp.ids <- data.table(dbGetQuery(con, paste("select ID from sgp_data where CONTENT_AREA in ('", tail(sgp.iter[["sgp.content.areas"]], 1), "')", 
					" AND GRADE in ('", tail(sgp.iter[["sgp.grade.sequences"]], 1), "')",
					" AND YEAR in ('", tail(sgp.iter[["sgp.panel.years"]], 1), "')", sep="")))
				tmp.data <- data.table(data.table(dbGetQuery(con, 
					paste("select * from sgp_data where GRADE in ('", paste(sgp.iter[["sgp.grade.sequences"]], collapse="', '"), "')",
						" AND YEAR in ('", paste(sgp.iter[["sgp.panel.years"]], collapse="', '"), "')", sep="")), 
					key="ID")[tmp.ids], key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))
			} else {
				setkeyv(sgp.data, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", tail(sgp.iter[["sgp.panel.years.within"]], 1)))
				tmp.ids <- sgp.data[SJ("VALID_CASE", tail(sgp.iter[["sgp.content.areas"]], 1), tail(sgp.iter[["sgp.panel.years"]], 1), 
					tail(sgp.iter[["sgp.grade.sequences"]], 1), 1)][,"ID", with=FALSE]
				tmp.data <- data.table(sgp.data, key="ID")[tmp.ids]
			}
			tmp.lookup <- data.table(V1="VALID_CASE", tail(sgp.iter[[sgp.projection.content.areas.label]], length(sgp.iter[[sgp.projection.grade.sequences.label]])),
				head(sgp.iter[["sgp.panel.years"]], length(sgp.iter[[sgp.projection.grade.sequences.label]])), sgp.iter[[sgp.projection.grade.sequences.label]],
				head(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[[sgp.projection.grade.sequences.label]])), FIRST_OBSERVATION=as.integer(NA), LAST_OBSERVATION=as.integer(NA))
			tmp.lookup[grep("FIRST", V5, ignore.case=TRUE), FIRST_OBSERVATION:=1L]; tmp.lookup[grep("LAST", V5, ignore.case=TRUE), LAST_OBSERVATION:=1L]; tmp.lookup[,V5:=NULL]
			setnames(tmp.lookup, paste("V", 1:4, sep=""), c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))
			tmp.lookup.list <- list()
			for (i in unique(sgp.iter[["sgp.panel.years.within"]])) {
				if (sqlite.tf) {
					setkeyv(tmp.lookup, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					suppressWarnings(tmp.lookup.list[[i]] <- data.table(dbGetQuery(con, 
						paste("select * from sgp_data where CONTENT_AREA in ('", paste(tmp.lookup[get(i)==1]$CONTENT_AREA, collapse="', '"), "')", 
							" AND GRADE in ('", paste(tmp.lookup[get(i)==1]$GRADE, collapse="', '"), "')", " AND ", i, " in (1)",
							" AND YEAR in ('", paste(tmp.lookup[get(i)==1]$YEAR, collapse="', '"), "')", sep="")), 
						key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))[tmp.lookup[get(i)==1], nomatch=0][,
							'tmp.timevar':=paste(YEAR, CONTENT_AREA, i, sep="."), with=FALSE][, 
								c("ID", "GRADE", "SCALE_SCORE", "ACHIEVEMENT_LEVEL", "YEAR_WITHIN", "tmp.timevar", sgp.scale.score.equated), with=FALSE])
				} else {
					setkeyv(sgp.data, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					setkeyv(tmp.lookup, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", i))
					suppressWarnings(tmp.lookup.list[[i]] <- data.table(sgp.data[tmp.lookup[get(i)==1], nomatch=0][,'tmp.timevar':=paste(YEAR, CONTENT_AREA, i, sep="."), with=FALSE][,
						c("ID", "GRADE", "SCALE_SCORE", "ACHIEVEMENT_LEVEL", "YEAR_WITHIN", "tmp.timevar", sgp.scale.score.equated), with=FALSE], key="ID")) ### Could be NULL and result in a warning
				}
			}
			if (sqlite.tf) dbDisconnect(con)
			
			achievement.level.prior.vname <- paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep=".")	
			if (is.null(sgp.targets)) {
				tmp.data <- ddcast(rbindlist(tmp.lookup.list),
					ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", "ACHIEVEMENT_LEVEL", "YEAR_WITHIN", state, sgp.scale.score.equated, SGPt), sep=".")

				setnames(tmp.data, names(tmp.data)[grep(achievement.level.prior.vname, names(tmp.data))], achievement.level.prior.vname)
				setnames(tmp.data, tail(sort(grep("YEAR_WITHIN", names(tmp.data), value=TRUE)), 1), "YEAR_WITHIN")

				if (length(setdiff(grep("YEAR_WITHIN", names(tmp.data), value=TRUE), "YEAR_WITHIN")) > 0) {
					tmp.data[,setdiff(grep("YEAR_WITHIN", names(tmp.data), value=TRUE), "YEAR_WITHIN"):=NULL, with=FALSE]
				}
				
				if ("STATE" %in% var.names && dim(tmp.data)[1]!=0) {
					setnames(tmp.data, tail(sort(grep("STATE", names(tmp.data), value=TRUE)), 1), "STATE")
					if (length(setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE")) > 0) {
						tmp.data[,setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE"):=NULL, with=FALSE]
					}
				}
				return(tmp.data)
			} else {
				tmp.data <- ddcast(rbindlist(tmp.lookup.list),
					ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", "ACHIEVEMENT_LEVEL", "YEAR_WITHIN", state, sgp.scale.score.equated, SGPt), sep=".")[
						sgp.targets[CONTENT_AREA==tail(sgp.iter[["sgp.content.areas"]], 1) & 
						YEAR==tail(sgp.iter[["sgp.panel.years"]], 1)], nomatch=0][,!c("CONTENT_AREA", "YEAR"), with=FALSE]

				setnames(tmp.data, names(tmp.data)[grep(achievement.level.prior.vname, names(tmp.data))], achievement.level.prior.vname)
				setnames(tmp.data, tail(sort(grep("YEAR_WITHIN", names(tmp.data), value=TRUE)), 1), "YEAR_WITHIN")

				if (length(setdiff(grep("YEAR_WITHIN", names(tmp.data), value=TRUE), "YEAR_WITHIN")) > 0) {
					tmp.data[,setdiff(grep("YEAR_WITHIN", names(tmp.data), value=TRUE), "YEAR_WITHIN"):=NULL, with=FALSE]
				}
				
				if ("STATE" %in% var.names && dim(tmp.data)[1]!=0) {
					setnames(tmp.data, tail(sort(grep("STATE", names(tmp.data), value=TRUE)), 1), "STATE")
					if (length(setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE")) > 0) {
						tmp.data[,setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE"):=NULL, with=FALSE]
					}
				}
				return(tmp.data)
			}
		} else {
			if (is.null(sgp.targets)) {
				if (sqlite.tf) {
					tmp.ids <- data.table(dbGetQuery(con, paste("select ID from sgp_data where CONTENT_AREA in ('", tail(sgp.iter[["sgp.content.areas"]], 1), "')", 
						" AND GRADE in ('", tail(sgp.iter[["sgp.grade.sequences"]], 1), "')",
						" AND YEAR in ('", tail(sgp.iter[["sgp.panel.years"]], 1), "')", sep="")))
					tmp.data <- data.table(data.table(dbGetQuery(con, 
						paste("select * from sgp_data where GRADE in ('", paste(sgp.iter[["sgp.grade.sequences"]], collapse="', '"), "')",
							" AND YEAR in ('", paste(sgp.iter[["sgp.panel.years"]], collapse="', '"), "')", sep="")), 
						key="ID")[tmp.ids], key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))[SJ("VALID_CASE", sgp.iter[[sgp.projection.content.areas.label]],
						tail(head(sgp.iter[["sgp.panel.years"]], -1), length(sgp.iter[[sgp.projection.grade.sequences.label]])),
						sgp.iter[[sgp.projection.grade.sequences.label]]), nomatch=0][,'tmp.timevar':=paste(YEAR, CONTENT_AREA, sep="."), with=FALSE]
					tmp.data <- ddcast(tmp.data,
							ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", "ACHIEVEMENT_LEVEL", state, sgp.scale.score.equated, SGPt), sep=".")
					dbDisconnect(con)
				} else {
					tmp.data <- ddcast(
						data.table(
							data.table(sgp.data, key="ID")[
								sgp.data[SJ("VALID_CASE", 
								tail(sgp.iter[["sgp.content.areas"]], 1), 
								tail(sgp.iter[["sgp.panel.years"]], 1), 
								tail(sgp.iter[["sgp.grade.sequences"]], 1))][,"ID", with=FALSE]], key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))[
							SJ("VALID_CASE", sgp.iter[[sgp.projection.content.areas.label]],
								tail(head(sgp.iter[["sgp.panel.years"]], -1), length(sgp.iter[[sgp.projection.grade.sequences.label]])),
								sgp.iter[[sgp.projection.grade.sequences.label]]), nomatch=0][,
								'tmp.timevar' := paste(YEAR, CONTENT_AREA, sep="."), with=FALSE],
							ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", "ACHIEVEMENT_LEVEL", state, sgp.scale.score.equated, SGPt), sep=".")
				}

				if ("STATE" %in% var.names && dim(tmp.data)[1]!=0) {
					setnames(tmp.data, tail(sort(grep("STATE", names(tmp.data), value=TRUE)), 1), "STATE")
					if (length(setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE")) > 0) {
						tmp.data[,setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE"):=NULL, with=FALSE]
					}
				}
				return(tmp.data)
			} else {
				if (sqlite.tf) {
					tmp.ids <- data.table(dbGetQuery(con, paste("select ID from sgp_data where CONTENT_AREA in ('", tail(sgp.iter[["sgp.content.areas"]], 1), "')", 
						" AND GRADE in ('", tail(sgp.iter[["sgp.grade.sequences"]], 1), "')",
						" AND YEAR in ('", tail(sgp.iter[["sgp.panel.years"]], 1), "')", sep="")))
					tmp.data <- data.table(data.table(dbGetQuery(con, 
						paste("select * from sgp_data where GRADE in ('", paste(sgp.iter[["sgp.grade.sequences"]], collapse="', '"), "')",
							" AND YEAR in ('", paste(sgp.iter[["sgp.panel.years"]], collapse="', '"), "')", sep="")), 
						key="ID")[tmp.ids], key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))[SJ("VALID_CASE", sgp.iter[[sgp.projection.content.areas.label]],
						tail(head(sgp.iter[["sgp.panel.years"]], -1), length(sgp.iter[[sgp.projection.grade.sequences.label]])),
						sgp.iter[[sgp.projection.grade.sequences.label]]), nomatch=0][,'tmp.timevar':=paste(YEAR, CONTENT_AREA, sep="."), with=FALSE]

						tmp.data <- ddcast(tmp.data,
							ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", state, sgp.scale.score.equated, SGPt), sep=".")[
								sgp.targets[CONTENT_AREA==tail(sgp.iter[["sgp.content.areas"]], 1) & 
								YEAR==tail(sgp.iter[["sgp.panel.years"]], 1)], nomatch=0][, !c("CONTENT_AREA", "YEAR"), with=FALSE]
					dbDisconnect(con)
				} else {
					tmp.data <- ddcast(
						data.table(
							data.table(sgp.data, key="ID")[
								sgp.data[SJ("VALID_CASE", 
								tail(sgp.iter[["sgp.content.areas"]], 1), 
								tail(sgp.iter[["sgp.panel.years"]], 1), 
								tail(sgp.iter[["sgp.grade.sequences"]], 1))][,"ID", with=FALSE]], 
						key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))[
						SJ("VALID_CASE", sgp.iter[[sgp.projection.content.areas.label]], 
							tail(head(sgp.iter[["sgp.panel.years"]], -1), length(sgp.iter[[sgp.projection.grade.sequences.label]])),
							sgp.iter[[sgp.projection.grade.sequences.label]]), nomatch=0][,
							'tmp.timevar' := paste(YEAR, CONTENT_AREA, sep="."), with=FALSE],
						ID ~ tmp.timevar, value.var=c("GRADE", "SCALE_SCORE", state, sgp.scale.score.equated, SGPt), sep=".")[
							sgp.targets[CONTENT_AREA==tail(sgp.iter[["sgp.content.areas"]], 1) & 
							YEAR==tail(sgp.iter[["sgp.panel.years"]], 1)], nomatch=0][, !c("CONTENT_AREA", "YEAR"), with=FALSE]
				}

				if ("STATE" %in% var.names && dim(tmp.data)[1]!=0) {
					setnames(tmp.data, tail(sort(grep("STATE", names(tmp.data), value=TRUE)), 1), "STATE")
					if (length(setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE")) > 0) {
						tmp.data[,setdiff(grep("STATE", names(tmp.data), value=TRUE), "STATE"):=NULL, with=FALSE]
					}
				}
				return(tmp.data)
			}
		}
	} ### END if (sgp.type=="sgp.projections.lagged")
} ## END getPanelData
