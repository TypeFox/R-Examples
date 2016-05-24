`studentGrowthPlot` <-
function(Scale_Scores,                        ## Vector of Scale Scores
	Plotting_Scale_Scores,                ## Score used for plotting, if missing, then Scale_Scores are used for plotting,
                                              ## if supplied Scale_Scores used for text
	Achievement_Levels,                   ## NOTE: Achievement_Levels must/should be supplied as factors with appropriate level codings
	SGP,                                  ## Vector of SGPs
	SGP_Levels,                           ## Vector of SGP Levels
	Grades,                               ## Vector of Grade levels for student
	Content_Areas,                        ## Vector of Content Areas for student
	Cuts,                                 ## Vector of NY1, NY2, and NY3 cutscores
	Plotting_Cuts,                        ## Vector of NY1, NY2, and NY3 cutscores used for plotting (transformed if non-vertical/equated scale)
	SGP_Targets,                          ## Vector of CUKU, CUKU_Current, MUSU, MUSU_Current (multi) year targets
	SGP_Scale_Score_Targets,              ## Vector of CUKU, CUKU_Current, MUSU, MUSU_Current scale score targets
	Plotting_SGP_Scale_Score_Targets,     ## Vector of CUKU, CUKU_Current, MUSU, MUSU_Current scale score targets for plotting (transformed if non-vertical/equated scale)
	Cutscores,                            ## data.frame of long formatted achievement level cutscores
	Years,                                ## Vector of years corresponding to Scale_Scores, Content_Areas, ... arguments supplied
	Report_Parameters) {                  ## list containing Current_Year, Content_Area, Content_Area_Title, State, Denote_Content_Area, SGP_Targets, Configuration, Language, Assessment_Transition,
                                              ## Fan


	############################################
	### Create relevant variables
	############################################

	YEAR <- CONTENT_AREA <- GRADE <- CUTSCORES <- CUTLEVEL <- level_1_1_curve <- level_2_1_curve <- NULL ## To prevent R CMD check warnings

	if (is.null(Report_Parameters$Assessment_Transition)) {
		achievement.level.labels <- list(SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["Achievement_Level_Labels"]])
		number.achievement.level.regions <- sapply(achievement.level.labels, length)
		level.to.get.cuku <- list(which.max(SGP::SGPstateData[[Report_Parameters$State]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient")-1)
		level.to.get.musu <- list(which.max(SGP::SGPstateData[[Report_Parameters$State]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient"))
	} else {
		achievement.level.labels <- SGP::SGPstateData[[Report_Parameters$State]][["Assessment_Program_Information"]][["Assessment_Transition"]][
			grep("Achievement_Level_Labels", names(SGP::SGPstateData[[Report_Parameters$State]][["Assessment_Program_Information"]][["Assessment_Transition"]]))]
		number.achievement.level.regions <- sapply(achievement.level.labels, length)
		achievement.levels.proficiency <- SGP::SGPstateData[[Report_Parameters$State]][["Assessment_Program_Information"]][["Assessment_Transition"]][
			grep("Achievement_Levels", names(SGP::SGPstateData[[Report_Parameters$State]][["Assessment_Program_Information"]][["Assessment_Transition"]]))]
		level.to.get.cuku <- sapply(achievement.levels.proficiency, function(x) which.max(x[['Proficient']]=="Proficient")-1)
		level.to.get.musu <- sapply(achievement.levels.proficiency, function(x) which.max(x[['Proficient']]=="Proficient"))
	}
	number.growth.levels <- length(SGP::SGPstateData[[Report_Parameters$State]][["Growth"]][["Levels"]])
	growth.level.labels <- SGP::SGPstateData[[Report_Parameters$State]][["Growth"]][["Levels"]]
	growth.level.cutscores <- SGP::SGPstateData[[Report_Parameters$State]][["Growth"]][["Cutscores"]][["Cuts"]]
	growth.level.cutscores.text <- SGP::SGPstateData[[Report_Parameters$State]][["Growth"]][["Cutscores"]][["Labels"]]
	content.area.label <- SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["Content_Areas_Labels"]][[Report_Parameters$Content_Area_Title]]

	if (!is.null(SGP::SGPstateData[[Report_Parameters$State]][["SGP_Configuration"]][["content_area.projection.sequence"]][[Report_Parameters$Content_Area]])) {
		grades.content_areas.reported.in.state <- data.frame(
				GRADE=SGP::SGPstateData[[Report_Parameters$State]][["SGP_Configuration"]][["grade.projection.sequence"]][[Report_Parameters$Content_Area]],
				YEAR_LAG=c(1, SGP::SGPstateData[[Report_Parameters$State]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[Report_Parameters$Content_Area]]),
				CONTENT_AREA=SGP::SGPstateData[[Report_Parameters$State]][["SGP_Configuration"]][["content_area.projection.sequence"]][[Report_Parameters$Content_Area]],
				stringsAsFactors=FALSE
			)
	} else {
		grades.content_areas.reported.in.state <- data.frame(
				GRADE=SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["Grades_Reported"]][[Report_Parameters$Content_Area]],
				YEAR_LAG=c(1, diff(as.numeric(SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["Grades_Reported"]][[Report_Parameters$Content_Area]]))),
				CONTENT_AREA=Report_Parameters$Content_Area,
				stringsAsFactors=FALSE
			)
	}

	grades.content_areas.reported.in.state$GRADE_NUMERIC <- (as.numeric(grades.content_areas.reported.in.state$GRADE[2])-1)+c(0, cumsum(tail(grades.content_areas.reported.in.state$YEAR_LAG, -1)))

	test.abbreviation <- SGP::SGPstateData[[Report_Parameters$State]][["Assessment_Program_Information"]][["Assessment_Abbreviation"]]

	if (identical(SGP::SGPstateData[[Report_Parameters$State]][['Assessment_Program_Information']][['Test_Season']], "Fall")) {
		test.season <- SGP::SGPstateData[[Report_Parameters$State]][['Assessment_Program_Information']][['Test_Season']]
	} else {
		test.season <- NULL
	}

	achievement.level.region.colors <- lapply(number.achievement.level.regions, function(x) paste("grey", round(seq(62, 91, length=x)), sep=""))

	border.color <- "grey25"
	if (is.null(SGP::SGPstateData[[Report_Parameters$State]][["SGP_Configuration"]][["arrow.legend.color"]])) {
		arrow.legend.color <- rev(diverge_hcl(number.growth.levels, h = c(180, 40), c = 255, l = c(20, 100)))
	} else {
		arrow.legend.color <- SGP::SGPstateData[[Report_Parameters$State]][["SGP_Configuration"]][["arrow.legend.color"]]
	}
	if (is.null(SGP::SGPstateData[[Report_Parameters$State]][["SGP_Configuration"]][["sgp.target.types"]])) {
		sgp.target.types <- c("Scale_Score_Targets_CUKU", "Scale_Score_Targets_MUSU", "Scale_Score_Targets_Current_CUKU", "Scale_Score_Targets_Current_MUSU")
	} else {
		sgp.target.types <- SGP::SGPstateData[[Report_Parameters$State]][["SGP_Configuration"]][["sgp.target.types"]]
	}
	missing.data.symbol <- "--"
	if (!is.null(my.year.span <- SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["sgPlot.year.span"]])) {
		studentGrowthPlot.year.span <- my.year.span
	} else studentGrowthPlot.year.span <- 5
	if (is.null(Report_Parameters[['Denote_Content_Area']]) || Report_Parameters[['Denote_Content_Area']]==FALSE) {
		legend.fill.color <- "white"
	} else {
		legend.fill.color <- rgb(0,0,1,0.25)
	}

	if (is.null(Report_Parameters[['Configuration']][['Connect_Points']])) {
		connect.points <- "Arrows"
	} else {
		connect.points <- Report_Parameters[['Configuration']][['Connect_Points']]
	}

	if (!is.null(Report_Parameters[['SGP_Targets']])) {
		if (all(c("sgp.projections", "sgp.projections.lagged") %in% Report_Parameters[['SGP_Targets']]) | all(c("sgp.projections.baseline", "sgp.projections.lagged.baseline") %in% Report_Parameters[['SGP_Targets']])) tmp.target.types <- names(unlist(SGP_Targets)[!is.na(unlist(SGP_Targets))])
		if (identical("sgp.projections", Report_Parameters[['SGP_Targets']]) | identical("sgp.projections.baseline", Report_Parameters[['SGP_Targets']])) tmp.target.types <- grep("Current", names(unlist(SGP_Targets)[!is.na(unlist(SGP_Targets))]), value=TRUE)
		if (identical("sgp.projections.lagged", Report_Parameters[['SGP_Targets']]) | identical("sgp.projections.lagged.baseline", Report_Parameters[['SGP_Targets']])) tmp.target.types <- grep("Current", names(unlist(SGP_Targets)[!is.na(unlist(SGP_Targets))]), value=TRUE, invert=TRUE)
	}

	if (!is.null(SGP::SGPstateData[[Report_Parameters$State]][['SGP_Configuration']][['sgPlot.show.content_area.progression']]) |
		!is.null(Report_Parameters[['Configuration']][['Assessment_Transition']])) {
			sgPlot.show.content_area.progression <- SGP::SGPstateData[[Report_Parameters$State]][['SGP_Configuration']][['sgPlot.show.content_area.progression']]
	} else {
		if (is.null(Report_Parameters[['Configuration']][['Show_Content_Area_Progression']])) {
			if (length(unique(Content_Areas[!is.na(Content_Areas)])) > 1 || !all(Content_Areas[!is.na(Content_Areas)]==Report_Parameters$Content_Area)) {
				sgPlot.show.content_area.progression <- TRUE
			} else {
				sgPlot.show.content_area.progression <- FALSE
			}
		} else {
			sgPlot.show.content_area.progression <- Report_Parameters[['Configuration']][['Show_Content_Area_Progression']]
		}
	}

	if (is.null(SGP::SGPstateData[[Report_Parameters$State]][['SGP_Configuration']][['Show_Fan_Cut_Scores']])) {
		show.fan.cutscores <- FALSE
	} else {
		show.fan.cutscores <- TRUE
	}

	if (is.null(Report_Parameters[['Configuration']][['Font_Size']])) {
		title.ca.size <- 1.8
		legend.size <- 0.5
		bottom.right.vp.size <- 1.2
		bottom.left.vp.size <- 0.7
	} else {
		title.ca.size <- Report_Parameters[['Configuration']][['Font_Size']][['title.ca.size']]
		legend.size <- Report_Parameters[['Configuration']][['Font_Size']][['legend.size']]
		bottom.right.vp.size <- Report_Parameters[['Configuration']][['Font_Size']][['bottom.right.vp.size']]
		bottom.left.vp.size <- Report_Parameters[['Configuration']][['Font_Size']][['bottom.left.vp.size']]
	}

	if (is.null(Report_Parameters[['Configuration']][['Language']])) {
		achievement.label <- "Achievement"
		achievement_level.label <- "Achievement Level"
		achievement_target.label <- "Achievement Target"
		growth.label <- "Growth"
		growth_percentile.label <- "Growth Percentile"
		growth_level.label <- "Growth Level"
		growth_target.label <- "Growth Target"
		level.label <- "Level"
		percentiles.label <- "Percentiles"
		scale_score.label <- "Scale Score"
		grade.label <- "Grade"
		CU.label <- "Catch Up"
		KU.label <- "Keep Up"
		MU.label <- "Move Up"
		SU.label <- "Stay Up"
		target.label <- "Target"
	}

	if (identical(toupper(Report_Parameters[['Configuration']][['Language']]), "SPANISH")) {
		achievement.label <- "Resultado"
		achievement_level.label <- "Nivel de Capacitaci\uF3n" # Hex code for accented o is \uF3 - http://www.ascii.cl/htmlcodes.htm
		achievement_target.label <- "Meta de Capacitaci\uF3n"
		growth.label <- "Crecimiento"
		growth_percentile.label <- "Porcentaje de Crecimiento"
		growth_level.label <- "Nivel de Crecimiento"
		growth_target.label <- "Meta de Crecimiento"
		level.label <- "Niveles"
		percentiles.label <- "El Percentil"
		scale_score.label <- "Escala"
		grade.label <- "Grado"
		CU.label <- "Alcanzar"
		KU.label <- "Guadar"
		MU.label <- "Avanzar"
		SU.label <- "Mantener"
		target.label <- "Meta"
		SGP_Levels <- SGP::SGPstateData[[Report_Parameters$State]][["Growth"]][["Levels"]][match(SGP_Levels, SGP::SGPstateData[[paste(head(unlist(strsplit(Report_Parameters$State, "_")), -1), collapse="_")]][["Growth"]][["Levels"]])]
	}

	if (length(content_area.label.pieces <- strsplit(content.area.label, " ")[[1]])==1) split.content_area.tf <- FALSE else split.content_area.tf <- TRUE

	if (is.null(Report_Parameters[['Fan']])) {
		show.fan <- TRUE
	} else {
		show.fan <- eval(parse(text=Report_Parameters[['Fan']]))
	}

	Cutscores <- copy(Cutscores)


	##############################
	### Utility functions
	##############################

	ach.level.labels <- function(perlevel){
		tmp <- unlist(sapply(achievement.level.labels, names), use.names=FALSE)[match(perlevel, unlist(achievement.level.labels))]
		tmp[is.na(tmp) & !is.na(perlevel)] <- perlevel[is.na(tmp) & !is.na(perlevel)]
		tmp[is.na(tmp)] <- missing.data.symbol
		return(tmp)
	}

	sgp.level.labels <- function(sgp_level){
		sgp_level[is.na(sgp_level)] <- missing.data.symbol
		return(sgp_level)
	}

	arrow.color <- function(sgp){
		arrow.legend.color[findInterval(sgp, growth.level.cutscores)+1]
	}

	get.my.cutscore.year <- function(state, content_area, year) {
		if (!is.null(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Transformed_Achievement_Level_Cutscores"]][[content_area]])) {
			return(NA)
		} else {
			year <- tail(sort(c(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Earliest_Year_Reported"]][[content_area]], year)), 1)
			tmp.cutscore.years <-
				sapply(strsplit(names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]])[grep(content_area, names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]]))], "[.]"), '[', 2)
			if (year %in% tmp.cutscore.years) {
				return(year)
			} else {
				if (year==sort(c(year, tmp.cutscore.years))[1]) {
					return(NA)
				} else {
					return(sort(tmp.cutscore.years)[which(year==sort(c(year, tmp.cutscore.years)))-1])
				}
			}
		}
	}

	interpolate.grades <- function(grades, content_areas, data.year.span) {
		last.number <- function (x) {
			if (sum(!is.na(x)) > 0) return(max(which(!is.na(x)))) else return (0)
		}

		first.number <- function (x) {
			if (sum(!is.na(x)) > 0 ) return(min(which(!is.na(x)))) else return (0)
		}

		convert.grades <- function(grades, content_areas, to="GRADE_NUMERIC", lookup=grades.content_areas.reported.in.state) {
			if (to=="GRADE_NUMERIC") {
				return(as.numeric(lookup$GRADE_NUMERIC[match(paste(grades, content_areas), paste(lookup$GRADE, lookup$CONTENT_AREA))]))
			}
			if (to=="GRADE") {
				return(as.character(lookup$GRADE[match(grades, lookup$GRADE_NUMERIC)]))
			}
		}

		extend.grades <- function (x) {
			tmp.grades.numeric <- x
			tmp.content_areas <- grades.content_areas.reported.in.state$CONTENT_AREA[match(x, grades.content_areas.reported.in.state$GRADE_NUMERIC)]
			tmp.grades <- convert.grades(x, to="GRADE")
			tmp.head <- min(which(x[1] <= grades.content_areas.reported.in.state$GRADE_NUMERIC), na.rm=TRUE)
			tmp.tail <- min(which(tail(x, 1) <= grades.content_areas.reported.in.state$GRADE_NUMERIC), na.rm=TRUE)
			if (tmp.head==1) {
				tmp.grades <- c("GRADE_LOWER", tmp.grades); tmp.content_areas <- c("PLACEHOLDER", tmp.content_areas); tmp.grades.numeric <- c(NA, tmp.grades.numeric)
			} else {
				tmp.grades <- c(grades.content_areas.reported.in.state$GRADE[tmp.head-1], tmp.grades)
				tmp.content_areas <- c(grades.content_areas.reported.in.state$CONTENT_AREA[tmp.head-1], tmp.content_areas)
				tmp.grades.numeric <- c(grades.content_areas.reported.in.state$GRADE_NUMERIC[tmp.head-1], tmp.grades.numeric)
			}
			if (tmp.tail==dim(grades.content_areas.reported.in.state)[1]) {
				tmp.grades <- c(tmp.grades, "GRADE_UPPER");
				tmp.content_areas <- c(tmp.content_areas, "PLACEHOLDER");
				tmp.grades.numeric <- c(tmp.grades.numeric, NA)
			} else {
				tmp.grades <- c(tmp.grades, grades.content_areas.reported.in.state$GRADE[tmp.tail+1])
				tmp.content_areas <- c(tmp.content_areas, grades.content_areas.reported.in.state$CONTENT_AREA[tmp.tail+1])
				tmp.grades.numeric <- c(tmp.grades.numeric, grades.content_areas.reported.in.state$GRADE_NUMERIC[tmp.tail+1])
			}
			tmp.grades[is.na(tmp.grades)] <- tmp.grades.numeric[is.na(tmp.grades)]
			return(data.frame(GRADE=tmp.grades, GRADE_NUMERIC=tmp.grades.numeric, CONTENT_AREA=tmp.content_areas, stringsAsFactors=FALSE))
		}

		first.scale.score <- first.number(head(grades, data.year.span-1))
		last.scale.score <- last.number(grades)
		grades <- convert.grades(grades, content_areas)

		if (first.scale.score == 0) {
			year_span <- 0
			return (list(
				interp.df = data.frame(GRADE=head(unique(Cutscores$GRADE), 7), CONTENT_AREA=Report_Parameters$Content_Area, stringsAsFactors=FALSE),
				year_span=year_span,
				years=yearIncrement(Report_Parameters$Current_Year, -5:1)))
		} else {
			if (last.scale.score < data.year.span) {
				grades[(last.scale.score+1):data.year.span] <- (grades[last.scale.score]-1):(grades[last.scale.score] - (data.year.span - last.scale.score))
				grades[grades < min(grades.content_areas.reported.in.state$GRADE_NUMERIC)] <- min(grades.content_areas.reported.in.state$GRADE_NUMERIC)
			}

			if (first.scale.score > 1) {
				grades[1:(first.scale.score-1)] <- (grades[first.scale.score] + (first.scale.score - 1)):(grades[first.scale.score]+1)
				grades[grades > max(grades.content_areas.reported.in.state$GRADE_NUMERIC)] <- max(grades.content_areas.reported.in.state$GRADE_NUMERIC)
				if (any(is.na(grades))) {
					grades[which(is.na(grades))] <- approx(grades, xout=which(is.na(grades)))$y
					grades <- as.integer(grades)
				}
				if (!grades[1] %in% grades.content_areas.reported.in.state$GRADE_NUMERIC) {
					grades[1] <- grades.content_areas.reported.in.state$GRADE_NUMERIC[which.min(grades[1] > grades.content_areas.reported.in.state$GRADE_NUMERIC)-1]
				}
				if (any(!grades %in% grades.content_areas.reported.in.state$GRADE_NUMERIC)) {
					for (tmp.missing.grades in which(!grades %in% grades.content_areas.reported.in.state$GRADE_NUMERIC)) {
						grades[tmp.missing.grades] <-
							grades.content_areas.reported.in.state$GRADE_NUMERIC[which.min(grades[tmp.missing.grades] > grades.content_areas.reported.in.state$GRADE_NUMERIC)-1]
					}
				}
			}

			if (any(is.na(grades))) {
				tmp.na <- which(is.na(grades))
				for (i in tmp.na) {
					grades[i] <- grades.content_areas.reported.in.state$GRADE_NUMERIC[match(grades[i-1], grades.content_areas.reported.in.state$GRADE_NUMERIC)]
				}
				if (length(intersect(tmp.na, which(!is.na(suppressWarnings(as.numeric(grades)))))) > 0) {
					tmp.indices <- intersect(tmp.na, which(!is.na(suppressWarnings(as.numeric(grades)))))
					grades[tmp.indices] <- NA
					grades[tmp.indices] <- round(approx(suppressWarnings(as.numeric(grades)), xout=tmp.indices)$y)
				}
			}

			if (grades[1] == max(grades.content_areas.reported.in.state$GRADE_NUMERIC)) {
				year_span <- data.year.span
				temp.grades.content_areas <- extend.grades(rev(grades))
				return (list(
					interp.df = temp.grades.content_areas,
					year_span=year_span,
					increment_for_projection_current=0,
					years=yearIncrement(Report_Parameters$Current_Year, seq(1-max(which(grades[1]==temp.grades.content_areas$GRADE_NUMERIC)), length=dim(temp.grades.content_areas)[1]))))
			} else {
				year.increment.for.projection.current <- grades.content_areas.reported.in.state$YEAR_LAG[which(grades[1]==grades.content_areas.reported.in.state$GRADE_NUMERIC)+1]
				year_span <- max(min(last.scale.score, data.year.span-1), min(grades[1]-min(grades.content_areas.reported.in.state$GRADE_NUMERIC)+1, data.year.span-1))-
					(year.increment.for.projection.current-1)
				temp.grades <- c(rev(head(grades, year_span)),
					head(seq(grades.content_areas.reported.in.state$GRADE_NUMERIC[match(grades[1], grades.content_areas.reported.in.state$GRADE_NUMERIC)]+1, length=data.year.span),
						data.year.span-year_span))
				temp.grades.content_areas <- extend.grades(temp.grades)
				return (list(
					interp.df = temp.grades.content_areas,
					year_span=year_span,
					increment_for_projection_current=year.increment.for.projection.current,
					years=yearIncrement(Report_Parameters$Current_Year, seq(1-max(which(grades[1]==temp.grades.content_areas$GRADE_NUMERIC)), length=dim(temp.grades.content_areas)[1]))))
			}
				}
	} ### END interpolate.grades function

	year.function <- function(year, add.sub, vec.length, output.type="numeric", season=NULL) {
		if (is.null(season)) {
			if (length(grep("_", year) > 0)) {
				tmp <- as.numeric(unlist(strsplit(as.character(year), "_")))+add.sub
				if (output.type=="numeric") {
					return(seq(from=tmp[2], length=vec.length))
				} else {
					return(paste(seq(from=tmp[1], length=vec.length), "-", seq(from=tmp[2], length=vec.length), sep=""))
				}
			} else {
				return(seq(from=as.numeric(year)+add.sub, length=vec.length))
			}
		} else {
			if (length(grep("_", year) > 0)) {
		 		tmp <- as.numeric(unlist(strsplit(as.character(year), "_")))+add.sub
				if (output.type=="numeric") {
					return(seq(from=tmp[2], length=vec.length))
				} else {
					return(paste(season, seq(from=tmp[1], length=vec.length)))
				}
			} else {
				return(paste(season, seq(from=as.numeric(year)+add.sub, length=vec.length)))
			}
		}
	}

	stextGrob <- function (label, r=0.1, x = x, y = y, 
		just = "centre", hjust = NULL, vjust = NULL, rot = 0, check.overlap = FALSE, 
		default.units = "native", name = NULL, gp = gpar(), vp = NULL){
		# http://stackoverflow.com/questions/7734535/control-font-thickness-without-changing-font-size

		let <- textGrob("a", gp=gp, vp=vp)
		wlet <- grobWidth(let)
		hlet <- grobHeight(let)

		tg <- textGrob(label=label, x=x, y=y, gp=gpar(col="white"), just = just, hjust = hjust, vjust = vjust, rot = rot,
				check.overlap = check.overlap, default.units = default.units)

		tgl <- c(lapply(seq(0, 2*pi, length=36), function(theta){
		  textGrob(label=label,x=x+cos(theta)*r*wlet, y=y+sin(theta)*r*hlet, gp=gpar(col="black"),
				just = just, hjust = hjust, vjust = vjust, rot = rot, check.overlap = check.overlap, default.units = default.units)
		  }), list(tg))

		g <- gTree(children=do.call(gList, tgl), vp=vp, name=name, gp=gp)
	}

	grid.stext <- function(...){
		g <- stextGrob(...)
		grid.draw(g)
		invisible(g)
	}

	###
	### END Utility functions
	###

	grade.values <- interpolate.grades(Grades, Content_Areas, studentGrowthPlot.year.span)

	if (!is.null(Report_Parameters[['Assessment_Transition']][['Assessment_Transition_Type']])) {
		if (identical(toupper(Report_Parameters[['Assessment_Transition']][['Assessment_Transition_Type']][1]), "NO")) {
			tmp.cutscore.year <- tail(head(sort(unique(Cutscores$YEAR), na.last=FALSE), -1), 1)
			tmp.cutscore.grade <- Grades[which(Years==Report_Parameters[['Assessment_Transition']][['Year']])]
			if (is.na(tmp.cutscore.grade)) {
				tmp.cutscore.grade <- grade.values[['interp.df']][['GRADE']][which(grade.values[['years']]==Report_Parameters[['Assessment_Transition']][['Year']])]
			}
			if (is.na(tmp.cutscore.year)) {
				Cutscores[is.na(YEAR), CUTSCORES:=Cutscores[GRADE==tmp.cutscore.grade & is.na(YEAR)][['CUTSCORES']]]
			} else {
				Cutscores[is.na(YEAR) | YEAR < Report_Parameters$Assessment_Transition$Year,
					CUTSCORES:=Cutscores[GRADE==tmp.cutscore.grade & YEAR==tmp.cutscore.year][['CUTSCORES']]]
			}
		}
		if (identical(toupper(Report_Parameters[['Assessment_Transition']][['Assessment_Transition_Type']][2]), "NO")) {
			tmp.cutscore.grade <- Grades[which(Years==Report_Parameters[['Assessment_Transition']][['Year']])]
			if (is.na(tmp.cutscore.grade)) {
				tmp.cutscore.grade <- grade.values[['interp.df']][['GRADE']][which(grade.values[['years']]==Report_Parameters[['Assessment_Transition']][['Year']])]
			}
			Cutscores[!is.na(YEAR) & YEAR >= Report_Parameters$Assessment_Transition$Year,
				CUTSCORES:=Cutscores[GRADE==tmp.cutscore.grade & YEAR >= Report_Parameters$Assessment_Transition$Year][['CUTSCORES']]]
		}
	}

	if (grade.values$year_span > 0) {
		low.year <- year.function(Report_Parameters$Current_Year, (1-grade.values$year_span), 1)
		high.year <- year.function(Report_Parameters$Current_Year, studentGrowthPlot.year.span-grade.values$year_span, 1)
		year.text <- c(year.function(Report_Parameters$Current_Year, (1-grade.values$year_span), grade.values$year_span+grade.values$increment_for_projection_current, "character", test.season),
			rep(" ", studentGrowthPlot.year.span))
		year.text <- head(year.text, studentGrowthPlot.year.span)
		content_area.text <- grade.values$interp.df$CONTENT_AREA[match(gsub("-", "_", year.text), grade.values$years)]
		content_area.text[is.na(content_area.text)] <- " "
		for (i in which(content_area.text %in% names(SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["Content_Areas_Labels"]]))) {
			content_area.text[i] <- SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["Content_Areas_Labels"]][[content_area.text[i]]]
		}

		if (grade.values$increment_for_projection_current > 0) {
			grades.text.numbers <- c(Grades[grade.values$year_span:1],
				grades.content_areas.reported.in.state$GRADE[match(Grades[1], grades.content_areas.reported.in.state$GRADE)+seq(grade.values$increment_for_projection_current)])
			tmp.grades.text.numbers <- head(grade.values$interp.df$GRADE[-1], studentGrowthPlot.year.span)
		} else {
			grades.text.numbers <- Grades[grade.values$year_span:1]
			tmp.grades.text.numbers <- head(grade.values$interp.df$GRADE[-1], studentGrowthPlot.year.span)
		}
		if (!is.null(Report_Parameters[['Configuration']][['Zero_to_K']])) {
			grades.text.numbers[which(grades.text.numbers==0)] <- "K"
		}
		grades.text.numbers.missing <- which(is.na(grades.text.numbers))
		grades.text.numbers.non.tested <- which(!tmp.grades.text.numbers %in% grades.content_areas.reported.in.state$GRADE)
		grades.text.eoct <- which(grades.text.numbers=="EOCT")
		grades.text <- c(paste(grade.label, grades.text.numbers), rep(" ", studentGrowthPlot.year.span))
		grades.text[grades.text.numbers.non.tested] <- "Non-tested Grade"
		grades.text[grades.text.numbers.missing] <- missing.data.symbol
		grades.text[grades.text.eoct] <- "EOCT"
		grades.text <- head(grades.text, studentGrowthPlot.year.span)

		if (missing(Plotting_Scale_Scores)) Plotting_Scale_Scores <- Scale_Scores
		scale.scores.values <- c(Plotting_Scale_Scores[grade.values$year_span:1], rep(NA, studentGrowthPlot.year.span))
		scale.scores.values <- head(scale.scores.values, studentGrowthPlot.year.span)

		scale.scores.text <- c(Scale_Scores[grade.values$year_span:1], rep(" ", studentGrowthPlot.year.span))
		scale.scores.text[which(is.na(Scale_Scores[grade.values$year_span:1]))] <- missing.data.symbol
		scale.scores.text <- head(scale.scores.text, studentGrowthPlot.year.span)

		ach.levels.text <- c(ach.level.labels(Achievement_Levels[grade.values$year_span:1]), rep(" ", studentGrowthPlot.year.span))
		ach.levels.text <- head(ach.levels.text, studentGrowthPlot.year.span)

		if (grade.values$year_span > 1) {
			gp.values <- c(SGP[(grade.values$year_span-1):1], rep(NA, studentGrowthPlot.year.span))
			gp.values <- head(gp.values, studentGrowthPlot.year.span-1)

			gp.text <- c(SGP[(grade.values$year_span-1):1], rep(" ", studentGrowthPlot.year.span))
			gp.text[which(is.na(SGP[(grade.values$year_span-1):1]))] <- missing.data.symbol
			gp.text <- head(gp.text, studentGrowthPlot.year.span-1)

			gp.levels.text <- c(sgp.level.labels(SGP_Levels[(grade.values$year_span-1):1]), rep(" ", studentGrowthPlot.year.span))
			gp.levels.text <- head(gp.levels.text, studentGrowthPlot.year.span-1)
		} else {
			gp.values <- rep(NA, studentGrowthPlot.year.span-1)
			gp.text <- rep(" ", studentGrowthPlot.year.span-1)

			gp.levels.text <- rep(" ", studentGrowthPlot.year.span-1)
		}

		if (grade.values$increment_for_projection_current==0) {
			cuts.ny1 <- rep(NA, number.growth.levels+1)
			cuts.ny1.text <- rep(NA, number.growth.levels+1)
		} else {
			cuts.ny1 <- Plotting_Cuts[["NY1"]]
			cuts.ny1.text <- Cuts[["NY1"]]
		}
	}

	if (grade.values$year_span==0) {
		low.year <- year.function(Report_Parameters$Current_Year, 0, 1)
		high.year <- year.function(Report_Parameters$Current_Year, studentGrowthPlot.year.span-1, 1)
		year.text <- rep(" ", studentGrowthPlot.year.span)
		content_area.text <- rep(" ", studentGrowthPlot.year.span)

		grades.text <- rep(" ", studentGrowthPlot.year.span)

		scale.scores.values <- rep(NA, studentGrowthPlot.year.span)
		scale.scores.text <- rep(" ", studentGrowthPlot.year.span)

		ach.levels.text <- rep(" ", studentGrowthPlot.year.span)

		gp.values <- rep(NA, studentGrowthPlot.year.span-1)
		gp.text <- rep(" ", studentGrowthPlot.year.span-1)

		gp.levels.text <- rep(" ", studentGrowthPlot.year.span-1)

		cuts.ny1 <- rep(NA, number.growth.levels+1)
		cuts.ny1.text <- rep(NA, number.growth.levels+1)
	}

	current.year <- year.function(Report_Parameters$Current_Year, 0, 1)
	xscale.range <- range(low.year,high.year) + c(-0.075, 0.1)*diff(range(low.year,high.year))
	if (is.null(Report_Parameters$Assessment_Transition)) {
		tmp.year.cut <- NULL
		xscale.range.list <- list(xscale.range)
	} else {
		tmp.year.cut <- year.function(Report_Parameters$Assessment_Transition$Year, 0, 1)-0.5
		if (tmp.year.cut <= xscale.range[1]) xscale.range[1] <- tmp.year.cut-0.9
		xscale.range.list <- list(c(xscale.range[1], tmp.year.cut-0.025), c(tmp.year.cut+0.025, xscale.range[2]))
	}

	if (Report_Parameters$Content_Area_Title %in% names(SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["Transformed_Achievement_Level_Cutscores"]]) &&
		is.null(Report_Parameters$Assessment_Transition)) {
			tmp.range <-
				range(head(tail(SGP::SGPstateData[[Report_Parameters$State]][["Student_Report_Information"]][["Transformed_Achievement_Level_Cutscores"]][[Report_Parameters$Content_Area]],-1),-1), na.rm=TRUE)
		low.score <- min(cuts.ny1,
			Plotting_Scale_Scores,
			tmp.range,
			na.rm=TRUE)
		high.score <- max(cuts.ny1,
			Plotting_Scale_Scores,
			tmp.range,
			na.rm=TRUE)
		yscale.range <- extendrange(c(low.score, high.score), f=0.15)
	} else {
		low.score <- min(cuts.ny1,
			Plotting_Scale_Scores,
			Cutscores$CUTSCORES[Cutscores$GRADE==grade.values$interp.df$GRADE[1] & Cutscores$CUTLEVEL==1],
			na.rm=TRUE)
		high.score <- max(cuts.ny1,
			Plotting_Scale_Scores,
			Cutscores$CUTSCORES[Cutscores$GRADE==tail(grade.values$interp.df$GRADE, 1) & Cutscores$CUTLEVEL %in% (number.achievement.level.regions-1)],
			na.rm=TRUE)
		yscale.range <- extendrange(c(low.score,high.score), f=0.15)
	}

	if (is.null(Report_Parameters$Assessment_Transition)) {
		subject.report.vp <- viewport(layout = grid.layout(2, 3, widths = unit(c(1.15, 5.4, 1.5)/8.05, rep("npc", 3)),
						heights = unit(c(2.45, 0.9)/3.35, rep("npc", 2))), gp=gpar(fill="transparent"))
	} else {
		subject.report.vp <- viewport(layout = grid.layout(2, 3, widths = unit(c(1.15, 5.65, 1.25)/8.05, rep("npc", 3)),
						heights = unit(c(2.45, 0.9)/3.35, rep("npc", 2))), gp=gpar(fill="transparent"))
	}

	growth.chart.vp <- viewport(name="growth.chart.vp",
		layout.pos.row=1, layout.pos.col=2,
		xscale=xscale.range,
		yscale=yscale.range,
		gp=gpar(fill="transparent"))

	if (is.null(Report_Parameters$Assessment_Transition)) {
		right.vp <- viewport(name="right.vp",
			layout.pos.row=1, layout.pos.col=3,
			xscale=c(0,1),
			yscale=c(0,1),
			gp=gpar(fill="transparent"))
	} else {
		right.vp <- viewport(name="right.vp",
			layout.pos.row=1, layout.pos.col=3,
			xscale=c(0,1),
			yscale=yscale.range,
			gp=gpar(fill="transparent"))
	}

	left.vp <- viewport(name="left.vp",
		layout.pos.row=1, layout.pos.col=1,
		xscale=c(0, 1),
		yscale=yscale.range,
		gp=gpar(fill="transparent"))

	if (is.null(Report_Parameters$Assessment_Transition)) {
		growth.and.margins.vp <- viewport(name="growth.and.margins.vp",
			layout.pos.row=1, layout.pos.col=1:2)
	} else {
		growth.and.margins.vp <- viewport(name="growth.and.margins.vp",
			layout.pos.row=1, layout.pos.col=1:3)
		bottom.right.vp.size <- 1.05
	}

	bottom.vp <- viewport(name="bottom.vp",
		layout.pos.row=2, layout.pos.col=2,
		xscale=xscale.range,
		yscale=c(0,3),
		gp=gpar(fill="transparent"))

	bottom.left.right.vp <- viewport(name="bottom.left.right.vp",
		layout.pos.row=2, layout.pos.col=1:3,
		xscale=c(0,1),
		yscale=c(0,3),
		gp=gpar(fill="transparent"))

	bottom.right.vp <- viewport(name="bottom.right.vp",
		layout.pos.row=2, layout.pos.col=3,
		xscale=c(0,1),
		yscale=c(0,3),
		gp=gpar(fill="transparent"))

	bottom.left.vp <- viewport(name="bottom.left.vp",
		layout.pos.row=2, layout.pos.col=1,
		xscale=c(0,1),
		yscale=c(0,3),
		gp=gpar(fill="transparent"))

	pushViewport(subject.report.vp)


	#############################################
	### Growth.Chart Viewport
	#############################################

	pushViewport(growth.chart.vp)

	for (j in seq(length(Report_Parameters$Assessment_Transition[['Year']])+1)) {
		for (i in seq(number.achievement.level.regions[[j]]-1)) {
			if (is.null(tmp.year.cut)) {
				tmp.year.sequence <- list(seq(length((low.year-1):(high.year+1))))
			} else {
				tmp.year.sequence <- list(match((low.year-1):tmp.year.cut, (low.year-1):(high.year+1)), match(rev((high.year+1):tmp.year.cut), (low.year-1):(high.year+1)))
			}
			temp <- cbind(temp_id=seq_len(nrow(grade.values$interp.df)), grade.values$interp.df, YEAR=grade.values$years)
			temp$YEAR <- sapply(temp$YEAR, function(x) get.my.cutscore.year(Report_Parameters$State, Report_Parameters$Content_Area, as.character(x)))
			temp <- merge(temp, subset(Cutscores, CUTLEVEL==i), all.x=TRUE)
			temp <- temp[order(temp$temp_id),]$CUTSCORES
			if (length(temp[which(!is.na(temp))])==1) {
				temp[which(is.na(temp))] <- temp[which(!is.na(temp))]
			} else {
				temp[which(is.na(temp))] <- approx(temp, xout=which(is.na(temp)), rule=2)$y
			}
			assign(paste("level_", j, "_", i, "_curve", sep=""), splinefun(((low.year-1):(high.year+1))[tmp.year.sequence[[j]]], temp[tmp.year.sequence[[j]]], method="natural"))
		}

		tmp.x.points <- seq(xscale.range.list[[j]][1], xscale.range.list[[j]][2], length=round(diff(xscale.range.list[[j]])/diff(xscale.range)*40))
		x.boundary.values.1 <- c(xscale.range.list[[j]][1], tmp.x.points, xscale.range.list[[j]][2])
		y.boundary.values.1 <- c(yscale.range[1], eval(parse(text=paste("level_", j, "_1_curve(tmp.x.points)", sep=""))), yscale.range[1])
		assign(paste("x.boundary.values.", number.achievement.level.regions[[j]], sep=""),
			c(xscale.range.list[[j]][1], tmp.x.points, xscale.range.list[[j]][2]))
		assign(paste("y.boundary.values.", number.achievement.level.regions[[j]], sep=""),
			c(yscale.range[2], eval(parse(text=paste("level_", j, "_", number.achievement.level.regions[[j]]-1, "_curve(tmp.x.points)", sep=""))), yscale.range[2]))

		if (number.achievement.level.regions[[j]] > 2) {
			for (i in 2:(number.achievement.level.regions[[j]]-1)) {
				assign(paste("x.boundary.values.", i, sep=""),
					c(tmp.x.points, rev(tmp.x.points)))

				assign(paste("y.boundary.values.", i, sep=""),
					eval(parse(text=paste("c(level_", j, "_", i, "_curve(tmp.x.points), level_", j, "_", i-1, "_curve(rev(tmp.x.points)))", sep=""))))
			}
		}

		## Keep achievement level bands from extending outside of chart box
		eval(parse(text=paste("y.boundary.values.", seq(number.achievement.level.regions[[j]]), "[y.boundary.values.", seq(number.achievement.level.regions[[j]]), " > max(yscale.range)] <- max(yscale.range)", sep="")))
		eval(parse(text=paste("y.boundary.values.", seq(number.achievement.level.regions[[j]]), "[y.boundary.values.", seq(number.achievement.level.regions[[j]]), " < min(yscale.range)] <- min(yscale.range)", sep="")))

		for (i in seq(number.achievement.level.regions[[j]])) {
			grid.polygon(x=get(paste("x.boundary.values.", i, sep="")),
				y=get(paste("y.boundary.values.", i, sep="")),
				default.units="native",
				gp=gpar(fill=achievement.level.region.colors[[j]][i], lwd=0.8, col="white"))
		}
	}

	if (!is.null(Report_Parameters$Assessment_Transition)) {
		grid.lines(x=tmp.year.cut, y=yscale.range, default.units="native", gp=gpar(lwd=1.8, col=border.color))
		grid.lines(x=tmp.year.cut+0.028, y=yscale.range, default.units="native", gp=gpar(lwd=0.4, col=border.color))
		grid.lines(x=tmp.year.cut-0.028, y=yscale.range, default.units="native", gp=gpar(lwd=0.4, col=border.color))

		if (grade.values$year_span != 0) {
			for (j in seq(length(Report_Parameters$Assessment_Transition[['Year']])+1)) {
				tmp.transition.names <- names(SGP::SGPstateData[[Report_Parameters$State]][["Assessment_Program_Information"]][["Assessment_Transition"]])
				tmp.test.abbreviation <-
					unlist(SGP::SGPstateData[[Report_Parameters$State]][["Assessment_Program_Information"]][["Assessment_Transition"]][grep('Assessment_Abbreviation', tmp.transition.names)])
				tmp.test.abbreviation.text <- rep(" ", length(year.text))
				tmp.test.abbreviation.text[which((low.year:high.year < Report_Parameters$Assessment_Transition$Year)[which(year.text!=" ")])] <-
					tmp.test.abbreviation[1]
				tmp.test.abbreviation.text[which((low.year:high.year >= Report_Parameters$Assessment_Transition$Year)[which(year.text!=" ")])] <-
					tmp.test.abbreviation[2]
				grid.stext(tmp.test.abbreviation.text, x=unit(low.year:high.year, "native"), y=convertY(unit(0.07, "npc"), "native"), gp=gpar(cex=0.7))
				grid.stext(sapply(content_area.text, capwords), x=unit(low.year:high.year, "native"), y=convertY(unit(0.03, "npc"), "native"), gp=gpar(cex=0.7))
			}
		}
	}

	if (grade.values$year_span == 0) {
		grid.text(x=0.5, y=0.5, paste("No", test.abbreviation, "Data"), gp=gpar(col=border.color, cex=2))
	}

	if (is.null(Report_Parameters$Assessment_Transition) && sgPlot.show.content_area.progression) {
		grid.stext(sapply(content_area.text, capwords), x=unit(low.year:high.year, "native"), y=convertY(unit(0.05, "npc"), "native"), gp=gpar(cex=0.75))
		# grid.text(x=low.year:high.year, y=convertY(unit(0.05, "npc"), "native"), sapply(content_area.text, capwords), gp=gpar(col="white", cex=0.75), default.units="native")
	}

	if (connect.points=="Arrows") {
		growth.arrow.coors.x <- c(.05, .85, .8, 1, .8, .85, .05, .053, .0555, .0575, .0585, .059, .0585, .0575, .0555, .053, .05)
		growth.arrow.coors.y <- c(-.2, -.2, -.5, 0, .5, .2,  seq(.2, -.2, length=11))

		for (i in 1:length(gp.values)) {
			tmp.lag <- which(is.na(rev(scale.scores.values[1:i]))==FALSE)
			if (!is.na(gp.values[i]) & length(tmp.lag) > 0) {
				lag.to.prior.score <- min(tmp.lag, na.rm=TRUE)
				if (lag.to.prior.score == 1) my.lty <- 1 else my.lty <- 2
				arrow.rise <- convertY(unit(scale.scores.values[i+1], "native") - unit(scale.scores.values[i+1-lag.to.prior.score], "native"), "inches")
				arrow.run <- convertX(unit(lag.to.prior.score, "native") - unit(0, "native"), "inches")
				arrow.angle <- atan2(as.numeric(arrow.rise),as.numeric(arrow.run))*180/pi

				## Arrows connecting achievement scores

				pushViewport(viewport(x=unit(low.year+i-lag.to.prior.score, "native"), y=unit(scale.scores.values[i+1-lag.to.prior.score], "native"),
					width=unit((lag.to.prior.score-0.08)/cos(arrow.angle*pi/180), "native"), height=unit(0.05, "npc"), angle=arrow.angle, just=c("left", "center"),
					xscale=c(0, 1), yscale=c(-0.5, 0.5)))
				grid.polygon(x=growth.arrow.coors.x, y=growth.arrow.coors.y, default.units="native",
					gp=gpar(lwd=0.3, lty=my.lty, col=border.color, fill=arrow.color(as.numeric(gp.values[i]))))
				grid.curve(0.05, -0.2, 0.05, 0.2, curvature=0.3, ncp=11, square=FALSE, default.units="native", gp=gpar(lwd=0.3, col=border.color))
				popViewport()
			}
		}
	} ## END Report_Parameters[['Configuration']][['Connect_Points']]=="Arrows"


	if (paste(Grades[1], Content_Areas[1]) != tail(with(grades.content_areas.reported.in.state, paste(GRADE, CONTENT_AREA)), 1) & !is.na(cuts.ny1[1]) & show.fan){
		for (i in seq(number.growth.levels)) {
			grid.polygon(x=c(current.year, rep(current.year+grade.values$increment_for_projection_current, 2), current.year),
				y=c(scale.scores.values[which(current.year==low.year:high.year)], max(yscale.range[1], cuts.ny1[i]),
				min(yscale.range[2], cuts.ny1[i+1]), scale.scores.values[which(current.year==low.year:high.year)]),
				default.units="native", gp=gpar(col=NA, lwd=0, fill=arrow.legend.color[i], alpha=0.45))
			grid.roundrect(x=unit(current.year+grade.values$increment_for_projection_current, "native"),
				y=unit((max(yscale.range[1], cuts.ny1[i])+min(yscale.range[2], cuts.ny1[i+1]))/2, "native"),
				height=unit(min(yscale.range[2], as.numeric(cuts.ny1[i+1])) - max(yscale.range[1], as.numeric(cuts.ny1[i])), "native"),
				width=unit(0.04, "native"), r=unit(0.45, "snpc"), gp=gpar(lwd=0.3, col=border.color, fill=arrow.legend.color[i]))

			if (is.null(Report_Parameters[['SGP_Targets']])) {
				grid.text(x=current.year+grade.values$increment_for_projection_current+0.05,
					y=(max(yscale.range[1], cuts.ny1[i])+min(yscale.range[2], cuts.ny1[i+1]))/2, growth.level.labels[i],
					default.units="native", just="left", gp=gpar(cex=0.4, col=border.color))
			}
			if (i != 1 & show.fan.cutscores) {
				grid.text(x=current.year+grade.values$increment_for_projection_current+0.05, y=cuts.ny1[i], as.character(cuts.ny1.text[i]),
					default.units="native", just="left", gp=gpar(cex=0.4, col=border.color))
			}
		}
	}

	if (!is.null(Report_Parameters[['SGP_Targets']])) {
		for (i in tmp.target.types) {
			if (length(grep("Current", i))==0) {
				current.year.x.coor <- current.year
				current.year.x.coor.lag <- min(which(!is.na(tail(Scale_Scores, -1))), na.rm=TRUE)
				x.coor.label.adjustment <- -0.075; label.position <- "right"; tmp.index <- 1
				if (is.null(Report_Parameters$Assessment_Transition)) {
					tmp.achievement.level <- which(tail(head(Achievement_Levels, current.year.x.coor.lag+1), 1)==achievement.level.labels[[tmp.index]])
				} else {
					if (year.function(Report_Parameters$Assessment_Transition$Year, 0, 1) <= current.year.x.coor-current.year.x.coor.lag) tmp.index <- 2
					tmp.achievement.level <- which(tail(head(Achievement_Levels, current.year.x.coor.lag+1), 1)==achievement.level.labels[[tmp.index]])
				}
				show.targets <- TRUE
			} else {
				current.year.x.coor <- current.year+grade.values$increment_for_projection_current
				current.year.x.coor.lag <- grade.values$increment_for_projection_current
				x.coor.label.adjustment <- 0.075; label.position <- "left"; tmp.index <- 1; show.targets <- show.fan
				if (is.null(Report_Parameters$Assessment_Transition)) {
					tmp.achievement.level <- which(head(Achievement_Levels, 1)==achievement.level.labels[[tmp.index]])
				} else {
					if (year.function(Report_Parameters$Assessment_Transition$Year, 0, 1) <= current.year.x.coor) tmp.index <- 2
					tmp.achievement.level <- which(head(Achievement_Levels, 1)==achievement.level.labels[[tmp.index]])
				}
			}

			if (show.targets) {
				if (length(grep("CUKU", i))>0 & tmp.achievement.level <= level.to.get.cuku[[tmp.index]]) {
					label.position <- c(label.position, "center")
					tmp.target.label <- c(CU.label, target.label)
					y.coordinates <- c(as.numeric(convertY(convertY(unit(Plotting_SGP_Scale_Score_Targets[[i]][['NY1']], "native"), "inches")+unit(0.0375, "inches"), "native")),
							as.numeric(convertY(convertY(unit(Plotting_SGP_Scale_Score_Targets[[i]][['NY1']], "native"), "inches")-unit(0.0375, "inches"), "native")))
				}
				if (length(grep("CUKU", i))>0 & tmp.achievement.level > level.to.get.cuku[[tmp.index]]) {
					label.position <- c(label.position, "top")
					tmp.target.label <- c(KU.label, target.label)
					y.coordinates <- c(Plotting_SGP_Scale_Score_Targets[[i]][['NY1']],
						as.numeric(convertY(convertY(unit(Plotting_SGP_Scale_Score_Targets[[i]][['NY1']], "native"), "inches")-unit(0.1, "inches"), "native")))
				}
				if (length(grep("MUSU", i))>0 & tmp.achievement.level <= level.to.get.musu[[tmp.index]]) {
					label.position <- c(label.position, "bottom")
					tmp.target.label <- c(MU.label, target.label)
					y.coordinates <- c(as.numeric(convertY(convertY(unit(Plotting_SGP_Scale_Score_Targets[[i]][['NY1']], "native"), "inches")+unit(0.1, "inches"), "native")),
						Plotting_SGP_Scale_Score_Targets[[i]][['NY1']])
				}
				if (length(grep("MUSU", i))>0 & tmp.achievement.level > level.to.get.musu[[tmp.index]]) {
					label.position <- c(label.position, "bottom")
					tmp.target.label <- c(SU.label, target.label)
					y.coordinates <- c(as.numeric(convertY(convertY(unit(Plotting_SGP_Scale_Score_Targets[[i]][['NY1']], "native"), "inches")+unit(0.1, "inches"), "native")),
						Plotting_SGP_Scale_Score_Targets[[i]][['NY1']])
				}
				grid.lines(x=c(current.year.x.coor-current.year.x.coor.lag, current.year.x.coor),
					y=c(scale.scores.values[which(current.year.x.coor-current.year.x.coor.lag==low.year:high.year)], Plotting_SGP_Scale_Score_Targets[[i]][['NY1']]),
					gp=gpar(lwd=0.8, col=border.color), default.units="native")
				grid.circle(x=current.year.x.coor, y=Plotting_SGP_Scale_Score_Targets[[i]][['NY1']], r=unit(c(0.05, 0.04, 0.025, 0.0125), "inches"),
					gp=gpar(col=c("black", "white", "black", "white"), lwd=0.01, fill=c("black", "white", "black", "white")), default.units="native")
				grid.text(x=current.year.x.coor+x.coor.label.adjustment, y=y.coordinates, tmp.target.label, default.units="native", just=label.position, gp=gpar(cex=0.5, col=border.color))
			}
		}
	}

	grid.circle(x=low.year:high.year, y=scale.scores.values, r=unit(0.04, "inches"), gp=gpar(col=border.color, lwd=0.7, fill="white"), default.units="native")

	popViewport()


	#################################
	### Left Viewport
	#################################

	pushViewport(left.vp)

	y.boundary.legend.1 <- c(yscale.range[1], yscale.range[1], rep(level_1_1_curve(xscale.range[1]), 2))
	assign(paste("y.boundary.legend.", number.achievement.level.regions[[1]], sep=""),
		c(yscale.range[2], yscale.range[2], rep(eval(parse(text=paste("level_1", "_", number.achievement.level.regions[[1]]-1, "_curve(xscale.range[1])", sep=""))), 2)))

	if (number.achievement.level.regions[[1]] > 2) {
		for (i in 2:(number.achievement.level.regions[[1]]-1)) {
			assign(paste("y.boundary.legend.", i, sep=""),
				eval(parse(text=paste("c(rep(level_1", "_", i-1, "_curve(xscale.range[1]), 2), rep(level_1", "_", i, "_curve(xscale.range[1]), 2))", sep=""))))
		}
	}

	for (i in seq(number.achievement.level.regions[[1]])){
	grid.polygon(x=c(0,1,1,0),
		y=get(paste("y.boundary.legend.", i, sep="")),
		default.units="native",
		gp=gpar(fill=achievement.level.region.colors[[1]][i], lwd=0.5, col=border.color, alpha=0.7))
	}

	grid.text(x=0.94, y=(level_1_1_curve(xscale.range[1]) + yscale.range[1])/2, names(achievement.level.labels[[1]])[1],
		gp=gpar(col=border.color, fontface=2, fontfamily="Helvetica-Narrow", cex=.85), default.units="native", just="right")
	grid.text(x=0.94, y=(eval(parse(text=paste("level_1", "_", number.achievement.level.regions[[1]]-1, "_curve(xscale.range[1])", sep=""))) + yscale.range[2])/2,
		names(achievement.level.labels[[1]])[number.achievement.level.regions[[1]]],
		gp=gpar(col=border.color, fontface=2, fontfamily="Helvetica-Narrow", cex=.85), default.units="native", just="right")

	if (number.achievement.level.regions[[1]] > 2) {
		for (i in 2:(number.achievement.level.regions[[1]]-1)) {
		grid.text(x=.94, y=(eval(parse(text=paste("(level_1", "_", i-1, "_curve(xscale.range[1]) + level_1", "_", i, "_curve(xscale.range[1]))/2", sep="")))),
			names(achievement.level.labels[[1]])[i],
			gp=gpar(col=border.color, fontface=2, fontfamily="Helvetica-Narrow", cex=.85), default.units="native", just="right")
		}
	}

	grid.lines(0, c(yscale.range[1], yscale.range[2]), gp=gpar(lwd=.8, col=border.color), default.units="native")

	popViewport()


	############################################
	### Growth and Margins Viewport
	############################################

	pushViewport(growth.and.margins.vp)
		if (grade.values$year_span == 0) {
			grid.roundrect(r=unit(.01, "snpc"), gp=gpar(lwd=1.8, col=border.color, clip=TRUE, fill=rgb(1, 1, 1, 0.5)))
		} else {
			grid.roundrect(r=unit(.01, "snpc"), gp=gpar(lwd=1.8, col=border.color, clip=TRUE))
		}
	popViewport()


	##############################################
	### Bottom Viewport
	##############################################

	pushViewport(bottom.vp)
	if (is.null(Report_Parameters[['SGP_Targets']])) {
		grid.text(x=low.year:high.year, y=2.67, grades.text, gp=gpar(col=border.color, cex=.75), default.units="native")
		grid.text(x=low.year:high.year, y=2.3, year.text, gp=gpar(col=border.color, cex=.6), default.units="native")
		grid.text(x=low.year:high.year, y=1.7, scale.scores.text, gp=gpar(col=border.color, cex=.65), default.units="native")
		grid.text(x=low.year:high.year, y=1.3, ach.levels.text, gp=gpar(col=border.color, cex=.6), default.units="native")
		grid.text(x=(low.year+1):high.year-0.5, y=0.7, gp.text, gp=gpar(col=border.color, cex=.65), default.units="native")
		grid.text(x=(low.year+1):high.year-0.5, y=0.3, gp.levels.text, gp=gpar(col=border.color, cex=.6), default.units="native")
	} else {
		grid.text(x=low.year:high.year, y=2.75, grades.text, gp=gpar(col=border.color, cex=.6), default.units="native")
		grid.text(x=low.year:high.year, y=2.45, year.text, gp=gpar(col=border.color, cex=.5), default.units="native")
		grid.text(x=low.year:high.year, y=1.95, scale.scores.text, gp=gpar(col=border.color, cex=.5), default.units="native")
		grid.text(x=low.year:high.year, y=1.65, ach.levels.text, gp=gpar(col=border.color, cex=.5), default.units="native")
		grid.text(x=(low.year+1):high.year-0.5, y=0.85, gp.text, gp=gpar(col=border.color, cex=.5), default.units="native")
		grid.text(x=(low.year+1):high.year-0.5, y=0.55, gp.levels.text, gp=gpar(col=border.color, cex=.5), default.units="native")

		tmp.projection.names.list <- list()
		if (length(grep("Current", tmp.target.types)) > 0) {
			tmp.projection.names.list[["Current"]] <- grep("Current", tmp.target.types, value=TRUE)
		}

		if (length(grep("Current", tmp.target.types, invert=TRUE)) > 0) {
			tmp.projection.names.list[["Lagged"]] <- grep("Current", tmp.target.types, value=TRUE, invert=TRUE)
		}

		for (i in seq_along(tmp.projection.names.list)) {
			if (length(grep("Current", tmp.projection.names.list[[i]])) > 0) {
				tmp.projection.names <- tmp.projection.names.list[[i]]
				tmp.projection.year.from <- current.year
				tmp.projection.year.to <- current.year+grade.values$increment_for_projection_current
				if (!is.null(Report_Parameters$Assessment_Transition) && year.function(Report_Parameters$Assessment_Transition$Year, 0, 1) <= tmp.projection.year.from) {
					achievement.level.label.index.from <- 2
				} else {
					achievement.level.label.index.from <- 1
				}
				if (!is.null(Report_Parameters$Assessment_Transition) && year.function(Report_Parameters$Assessment_Transition$Year, 0, 1) <= tmp.projection.year.to) {
					achievement.level.label.index.to <- 2
				} else {
					achievement.level.label.index.to <- 1
				}
				tmp.achievement.level <- which(head(Achievement_Levels, 1)==achievement.level.labels[[achievement.level.label.index.from]])
			} else {
				tmp.projection.names <- tmp.projection.names.list[[i]]
				tmp.projection.year.lag <- min(which(!is.na(tail(Scale_Scores, -1))), na.rm=TRUE)
				tmp.projection.year.from <- current.year-tmp.projection.year.lag
				tmp.projection.year.to <- current.year
				if (!is.null(Report_Parameters$Assessment_Transition) && year.function(Report_Parameters$Assessment_Transition$Year, 0, 1) <= tmp.projection.year.from) {
					achievement.level.label.index.from <- 2
				} else {
					achievement.level.label.index.from <- 1
				}
				if (!is.null(Report_Parameters$Assessment_Transition) && year.function(Report_Parameters$Assessment_Transition$Year, 0, 1) <= tmp.projection.year.to) {
					achievement.level.label.index.to <- 2
				} else {
					achievement.level.label.index.to <- 1
				}
				tmp.achievement.level <- which(tail(head(Achievement_Levels, tmp.projection.year.lag+1), 1)==achievement.level.labels[[achievement.level.label.index.from]])
			}

			if ((length(grep("CUKU", tmp.projection.names)) > 0 & tmp.achievement.level <= level.to.get.cuku[[achievement.level.label.index.from]]) |
				length(grep("MUSU", tmp.projection.names))==0) {
				level.to.get.cuku.label <- names(achievement.level.labels[[achievement.level.label.index.to]])[level.to.get.cuku[[achievement.level.label.index.to]]+1]
				grid.text(x=tmp.projection.year.to, y=1.35,
					paste(level.to.get.cuku.label, " (", SGP_Scale_Score_Targets[[grep("CUKU", tmp.projection.names, value=TRUE)]][['NY1']], ")", sep=""),
					gp=gpar(col=border.color, cex=.4), default.units="native")
				grid.text(x=tmp.projection.year.to, y=0.25,
					paste(CU.label, " (", SGP_Targets[[grep("CUKU", tmp.projection.names, value=TRUE)]], ")", sep=""),
					gp=gpar(col=border.color, cex=.4), default.units="native")
			} else {
				level.to.get.cuku.label <- names(achievement.level.labels[[achievement.level.label.index.to]])[level.to.get.cuku[[achievement.level.label.index.to]]+1]
				level.to.get.musu.label <- names(achievement.level.labels[[achievement.level.label.index.to]])[level.to.get.musu[[achievement.level.label.index.to]]+1]
				grid.text(x=tmp.projection.year.to, y=1.35,
					paste(level.to.get.cuku.label, " (", SGP_Scale_Score_Targets[[grep("CUKU", tmp.projection.names, value=TRUE)]][['NY1']], ")/", level.to.get.musu.label, " (", SGP_Scale_Score_Targets[[grep("MUSU", tmp.projection.names, value=TRUE)]][['NY1']], ")", sep=""),
					gp=gpar(col=border.color, cex=.4), default.units="native")
				if (tmp.achievement.level <= level.to.get.musu[[achievement.level.label.index.from]]) {
					grid.text(x=tmp.projection.year.to, y=0.25,
					paste(KU.label, " (", SGP_Targets[[grep('CUKU', tmp.projection.names, value=TRUE)]], ")/Move Up (", SGP_Targets[[grep('MUSU', tmp.projection.names, value=TRUE)]], ")", sep=""),
					gp=gpar(col=border.color, cex=.4), default.units="native")
				}
				if (tmp.achievement.level > level.to.get.musu[[achievement.level.label.index.from]]) {
					grid.text(x=tmp.projection.year.to, y=0.25,
					paste(KU.label, " (", SGP_Targets[[grep('CUKU', tmp.projection.names, value=TRUE)]], ")/Stay Up (", SGP_Targets[[grep('MUSU', tmp.projection.names, value=TRUE)]], ")", sep=""),
					gp=gpar(col=border.color, cex=.4), default.units="native")
				}
			}
		}
	}
	popViewport()


	##############################
	### Bottom Left.Right Viewport
	##############################

	pushViewport(bottom.left.right.vp)
	if (is.null(Report_Parameters[['SGP_Targets']])) {
		grid.lines(x=c(0,1), y=2, gp=gpar(lwd=1.8, col=border.color), default.units="native")
		grid.lines(x=c(0,1), y=1, gp=gpar(lwd=1, col=border.color), default.units="native")
		grid.lines(x=c(0,1), y=0, gp=gpar(lwd=1.8, col=border.color), default.units="native")
	} else {
		grid.lines(x=c(0,1), y=2.2, gp=gpar(lwd=1.8, col=border.color), default.units="native")
		grid.lines(x=c(0,1), y=1.1, gp=gpar(lwd=1, col=border.color), default.units="native")
		grid.lines(x=c(0,1), y=0, gp=gpar(lwd=1.8, col=border.color), default.units="native")
	}
	popViewport()


	################################
	### Bottom Right Viewport
	################################


	pushViewport(bottom.right.vp)
	if (is.null(Report_Parameters[['SGP_Targets']])) {
		grid.text(x=0.1, y=1.5, achievement.label, gp=gpar(col=border.color, cex=bottom.right.vp.size), just="left", default.units="native")
		grid.text(x=0.1, y=.5, growth.label, gp=gpar(col=border.color, cex=bottom.right.vp.size), just="left", default.units="native")
	} else {
		grid.text(x=0.1, y=1.65, achievement.label, gp=gpar(col=border.color, cex=bottom.right.vp.size+0.1), just="left", default.units="native")
		grid.text(x=0.1, y=0.55, growth.label, gp=gpar(col=border.color, cex=bottom.right.vp.size+0.1), just="left", default.units="native")
	}
	popViewport()


	###########################
	### Bottom Left Viewport
	###########################

	pushViewport(bottom.left.vp)
	if (is.null(Report_Parameters[['SGP_Targets']])) {
		grid.text(x=0.9, y=1.7, scale_score.label, gp=gpar(col=border.color, cex=bottom.left.vp.size), just="right", default.units="native")
		grid.text(x=0.9, y=1.3, achievement_level.label, gp=gpar(col=border.color, cex=bottom.left.vp.size), just="right", default.units="native")
		grid.text(x=0.9, y=0.7, growth_percentile.label, gp=gpar(col=border.color, cex=bottom.left.vp.size), just="right", default.units="native")
		grid.text(x=0.9, y=0.3, growth_level.label, gp=gpar(col=border.color, cex=bottom.left.vp.size), just="right", default.units="native")
	} else {
		grid.text(x=0.9, y=1.95, scale_score.label, gp=gpar(col=border.color, cex=bottom.left.vp.size-0.125), just="right", default.units="native")
		grid.text(x=0.9, y=1.65, achievement_level.label, gp=gpar(col=border.color, cex=bottom.left.vp.size-0.125), just="right", default.units="native")
		grid.text(x=0.9, y=1.35, achievement_target.label, gp=gpar(col=border.color, cex=bottom.left.vp.size-0.125), just="right", default.units="native")
		grid.text(x=0.9, y=0.85, growth_percentile.label, gp=gpar(col=border.color, cex=bottom.left.vp.size-0.125), just="right", default.units="native")
		grid.text(x=0.9, y=0.55, growth_level.label, gp=gpar(col=border.color, cex=bottom.left.vp.size-0.125), just="right", default.units="native")
		grid.text(x=0.9, y=0.25, growth_target.label, gp=gpar(col=border.color, cex=bottom.left.vp.size-0.125), just="right", default.units="native")
	}
	popViewport()


	###################################
	### Right Viewport
	###################################

	if (is.null(Report_Parameters$Assessment_Transition)) { ### right.vp without ASSESSMENT_TRANSITION
		pushViewport(right.vp)

		grid.roundrect(width=unit(0.95, "native"), r=unit(.02, "snpc"), gp=gpar(lwd=1.8, col=border.color, fill=legend.fill.color), just="center")

		if (!split.content_area.tf) {
			grid.text(x=.5, y=0.875, content.area.label, gp=gpar(col=border.color, cex=title.ca.size, fontface=2, fontfamily="Helvetica-Narrow"), default.units="native")
		} else {
			grid.text(x=.5, y=0.92, content_area.label.pieces[1], gp=gpar(col=border.color, cex=title.ca.size-0.25, fontface=2, fontfamily="Helvetica-Narrow"), default.units="native")
			grid.text(x=.5, y=0.83, content_area.label.pieces[2], gp=gpar(col=border.color, cex= title.ca.size-0.25, fontface=2, fontfamily="Helvetica-Narrow"), default.units="native")
		}
		grid.text(x=0.08, y=0.75, achievement.label, gp=gpar(col=border.color, cex=.85, fontface=2, fontfamily="Helvetica-Narrow"), default.units="native", just="left")
		grid.text(x=0.08, y=0.525, growth.label, gp=gpar(col=border.color, cex=.75, fontface=2, fontfamily="Helvetica-Narrow"), default.units="native", just="left")
		grid.text(x=0.275, y=0.455, level.label, gp=gpar(col=border.color, cex=.6), default.units="native", just="center")
		grid.text(x=0.75, y=0.455, percentiles.label, gp=gpar(col=border.color, cex=.6), default.units="native", just="center")

		grid.roundrect(x=unit(0.3, "native"), y=unit(0.64, "native"), width=unit(0.3, "native"), height=unit(0.1, "native"), r=unit(0.06, "char"),
			gp=gpar(col="grey72", lwd=0.4, fill="grey72"))
		grid.circle(x=0.3, y=0.64, r=unit(0.04, "inches"), gp=gpar(col=border.color, lwd=0.7, fill="white"), default.units="native")
		if (!split.content_area.tf) {
			grid.text(x=0.7, y=0.658, paste(test.abbreviation, content.area.label), gp=gpar(col=border.color, cex=legend.size), default.units="native")
			grid.text(x=0.7, y=0.622, scale_score.label, gp=gpar(col=border.color, cex= legend.size), default.units="native")
		} else {
			grid.text(x=0.7, y=0.7, test.abbreviation, gp=gpar(col=border.color, cex=legend.size), default.units="native")
			grid.text(x=0.7, y=0.665, content_area.label.pieces[1], gp=gpar(col=border.color, cex=legend.size), default.units="native")
			grid.text(x=0.7, y=0.63, content_area.label.pieces[2], gp=gpar(col=border.color, cex=legend.size), default.units="native")
			grid.text(x=0.7, y=0.595, scale_score.label, gp=gpar(col=border.color, cex= legend.size), default.units="native")
		}

		y.center <- seq(0.05, 0.4, length=number.growth.levels+1)
		arrow.legend.coors.x <- c(.25, .75, .75, 1, .5, 0, .25)
		arrow.legend.coors.y <- c(0, 0, 1.3, 1.1, 2, 1.1, 1.3)

		for (i in seq(number.growth.levels)) {
			pushViewport(viewport(x=unit(0.3, "native"), y=unit(y.center[i], "native"),
				width=unit(0.07, "native"), height=unit(y.center[2]-y.center[1], "npc"), just=c("center", "bottom"),
				xscale=c(0, 1), yscale=c(0, 2)))
			grid.polygon(x=arrow.legend.coors.x, y=arrow.legend.coors.y, default.units="native",
				gp=gpar(lwd=0.3, col=border.color, fill=arrow.legend.color[i]))

			popViewport()

			pushViewport(viewport(x=unit(0.2, "native"), y=unit(y.center[i], "native"),
				width=unit(0.04, "native"), height=unit(y.center[2]-y.center[1], "npc"), just=c("center", "bottom")))
				grid.roundrect(x=0.5, y=0.5, width=1, height=1, r=unit(.45, "snpc"),
				gp=gpar(lwd=0.3, col=border.color, fill=arrow.legend.color[i]))
			popViewport()

			grid.polygon(x=c(0.05, rep(0.1875, 2)),
				y=c((head(y.center,1)+tail(y.center,1))/2, y.center[i], y.center[i]+y.center[2]-y.center[1]), default.units="native",
				gp=gpar(col=NA, lwd=0, fill=arrow.legend.color[i], alpha=0.45))

			grid.text(x=0.375, y=((y.center[1]+y.center[2])/2)+(i-1)*(y.center[2]-y.center[1]), growth.level.labels[i], default.units="native",
				gp=gpar(col=border.color, cex=.5), just="left")
			grid.text(x=0.925, y=((y.center[1]+y.center[2])/2)+(i-1)*(y.center[2]-y.center[1]), growth.level.cutscores.text[i], default.units="native",
				gp=gpar(col=border.color, cex=.5), just="right")
		}
	} else { ### right.vp with ASSESSMENT_TRANSITION
		pushViewport(right.vp)

		y.boundary.legend.1 <- c(yscale.range[1], yscale.range[1], rep(level_2_1_curve(xscale.range[2]), 2))
		assign(paste("y.boundary.legend.", number.achievement.level.regions[[2]], sep=""),
			c(yscale.range[2], yscale.range[2], rep(eval(parse(text=paste("level_2", "_", number.achievement.level.regions[[2]]-1, "_curve(xscale.range[2])", sep=""))), 2)))

		if (number.achievement.level.regions[[2]] > 2) {
			for (i in 2:(number.achievement.level.regions[[2]]-1)) {
				assign(paste("y.boundary.legend.", i, sep=""),
					eval(parse(text=paste("c(rep(level_2", "_", i-1, "_curve(xscale.range[2]), 2), rep(level_2", "_", i, "_curve(xscale.range[2]), 2))", sep=""))))
			}
		}

		for (i in seq(number.achievement.level.regions[[2]])){
		grid.polygon(x=c(0,1,1,0),
			y=get(paste("y.boundary.legend.", i, sep="")),
			default.units="native",
			gp=gpar(fill=achievement.level.region.colors[[2]][i], lwd=0.5, col=border.color, alpha=0.7))
		}

		grid.text(x=0.06, y=(level_2_1_curve(xscale.range[2]) + yscale.range[1])/2, names(achievement.level.labels[[2]])[1],
			gp=gpar(col=border.color, fontface=2, fontfamily="Helvetica-Narrow", cex=.85), default.units="native", just="left")
		grid.text(x=0.06, y=(eval(parse(text=paste("level_2", "_", number.achievement.level.regions[[2]]-1, "_curve(xscale.range[2])", sep=""))) + yscale.range[2])/2,
			names(achievement.level.labels[[2]])[number.achievement.level.regions[[2]]],
			gp=gpar(col=border.color, fontface=2, fontfamily="Helvetica-Narrow", cex=.85), default.units="native", just="left")

		if (number.achievement.level.regions[[2]] > 2) {
			for (i in 2:(number.achievement.level.regions[[2]]-1)) {
			grid.text(x=.06, y=(eval(parse(text=paste("(level_2", "_", i-1, "_curve(xscale.range[2]) + level_2", "_", i, "_curve(xscale.range[2]))/2", sep="")))),
				names(achievement.level.labels[[2]])[i],
				gp=gpar(col=border.color, fontface=2, fontfamily="Helvetica-Narrow", cex=.85), default.units="native", just="left")
			}
		}

		grid.lines(0, c(yscale.range[1], yscale.range[2]), gp=gpar(lwd=.8, col=border.color), default.units="native")
	}
	popViewport(2)

} ## END studentGrowthPlot Function
