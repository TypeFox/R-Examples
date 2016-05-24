`bubblePlot_Styles` <- 
	function(sgp_object,
		state,
		bPlot.years=NULL,
		bPlot.content_areas=NULL,
		bPlot.districts=NULL,
		bPlot.schools=NULL,
		bPlot.instructors=NULL,
		bPlot.styles=c(1),
		bPlot.levels=NULL,
		bPlot.level.cuts=NULL,
		bPlot.full.academic.year=TRUE,
		bPlot.minimum.n=10,
		bPlot.anonymize=FALSE,
		bPlot.prior.achievement=TRUE, 
		bPlot.draft=FALSE,
		bPlot.demo=FALSE,
		bPlot.output ="PDF",
		bPlot.format="print",
		bPlot.folder="Visualizations/bubblePlots") {


	DISTRICT_NUMBER <- DISTRICT_NAME <- SCHOOL_NUMBER <- SCHOOL_NAME <- SCHOOL_ENROLLMENT_STATUS <- YEAR <- CONTENT_AREA <- MEDIAN_SGP_COUNT <- NULL ## To prevent R CMD check warnings
	ID <- YEAR_INTEGER_TMP <- SCALE_SCORE <- SGP <- GRADE <- NULL ## To prevent R CMD check warnings
	INSTRUCTOR_NUMBER <- INSTRUCTOR_NAME <- INSTRUCTOR_ENROLLMENT_STATUS <- NULL
	SGP_TARGET <- VALID_CASE <- ENROLLMENT_STATUS <- NULL
	SCALE_SCORE_PRIOR <- SGP_PRIOR <- SGP_TARGET_PRIOR <- ACHIEVEMENT_LEVEL_PRIOR <- CONTENT_AREA_PRIOR <- SGP_NORM_GROUP <- NULL

	### Define relevant quantities

	# State stuff

	if (state %in% c(datasets::state.abb, "DEMO")) {
		state.name.label <- c(datasets::state.name, "DEMONSTRATION")[state==c(datasets::state.abb, "DEMO")]
		test.abbreviation.label <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Abbreviation"]]
	} else {
		state.name.label <- test.abbreviation.label <- state
	}
		state.name.file.label <- gsub("_", " ", state.name.label)

	# draft message

	if (bPlot.draft) {
		bPlot.message <- c("grid.text(x=unit(50, 'native'), y=unit(50, 'native'), 'CONFIDENTIAL - DO NOT DISTRIBUTE', rot=-30, gp=gpar(col='grey80', cex=2.9, alpha=0.8, fontface=2))")
	} else {
		bPlot.message <- NULL
	}

	if (is.null(bPlot.level.cuts)) {
		bPlot.level.cuts <- seq(0,100,by=20)
	}

	### Utility functions	

	"%w/o%" <- function(x,y) x[!x %in% y]

	pretty_year <- function(x) sub("_", "-", x)

	create.bPlot.labels <- function(year.iter, y.variable.iter, bubblePlot_LEVEL) {
		pretty_year <- function(x) sub("_", "-", x)
		my.labels <- list()
		my.labels$x.year.label <- pretty_year(year.iter)
		if (length(grep("PRIOR", y.variable.iter)) > 0) {
			if (year.iter=="All Years") {
				my.labels$y.year <- "All Years"
			} else {
				my.labels$y.year <- paste(as.numeric(unlist(strsplit(as.character(year.iter), "_")))-1, collapse="_")
			}
			if (bubblePlot_LEVEL=="Summary") my.labels$y.year.label <- paste(pretty_year(my.labels$y.year), "Prior Percent at/above Proficient")
			if (bubblePlot_LEVEL=="Individual") my.labels$y.year.label <- list(PRIOR=paste(pretty_year(my.labels$y.year), "Achievement Level"), CURRENT=paste(pretty_year(year.iter), "Achievement Level"))
			if (bubblePlot_LEVEL=="Summary") my.labels$main.title <- paste(test.abbreviation.label, "Growth & Prior Achievement")
			if (bubblePlot_LEVEL=="Individual") my.labels$main.title <- paste(test.abbreviation.label, "Growth & Achievement")
			if (bubblePlot_LEVEL=="Summary") my.labels$pdf.title <- "Bubble_Plot_(Prior_Achievement)"
			if (bubblePlot_LEVEL=="Individual") my.labels$pdf.title <- "Student_Bubble_Plot"
		} else {
			my.labels$y.year <- year.iter
			if (bubblePlot_LEVEL=="Summary") my.labels$y.year.label <- paste(pretty_year(my.labels$y.year), "Percent at/above Proficient")
			if (bubblePlot_LEVEL=="Individual") my.labels$y.year.label <- paste(pretty_year(my.labels$y.year), "Achievement Level")
			if (bubblePlot_LEVEL=="Summary") my.labels$main.title <- paste(test.abbreviation.label, "Growth & Achievement")
			if (bubblePlot_LEVEL=="Individual") my.labels$main.title <- paste(test.abbreviation.label, "Growth & Achievement")
			if (bubblePlot_LEVEL=="Summary") my.labels$pdf.title <- "Bubble_Plot_(Current_Achievement)"
			if (bubblePlot_LEVEL=="Individual") my.labels$pdf.title <- "Student_Bubble_Plot_(Current_Achievement)"
		}
		return(my.labels)
	}

	names.merge <- function(tmp.data, bPlot.anonymize) {
		if (!"INSTRUCTOR_NUMBER" %in% names(tmp.data) & !"SCHOOL_NUMBER" %in% names(tmp.data) & "DISTRICT_NUMBER" %in% names(tmp.data)) {
			tmp.names <- unique(data.table(sgp_object@Data[!is.na(DISTRICT_NUMBER), 
				list(DISTRICT_NUMBER, DISTRICT_NAME, SCHOOL_NUMBER, SCHOOL_NAME, YEAR)], key=c("DISTRICT_NUMBER", "YEAR"))) # Keep other institution NUMBER to iterate over in some plots
			if (bPlot.anonymize) {
				tmp.names$DISTRICT_NAME <- paste("District", as.numeric(as.factor(tmp.names$DISTRICT_NUMBER)))
			}
			setkey(tmp.data, DISTRICT_NUMBER, YEAR)
		}
		
		if (!"INSTRUCTOR_NUMBER" %in% names(tmp.data) & "SCHOOL_NUMBER" %in% names(tmp.data) & !"DISTRICT_NUMBER" %in% names(tmp.data)) {
			tmp.names <- unique(data.table(sgp_object@Data[!is.na(SCHOOL_NUMBER), 
				list(DISTRICT_NUMBER, DISTRICT_NAME, SCHOOL_NUMBER, SCHOOL_NAME, YEAR)], key=c("SCHOOL_NUMBER", "YEAR"))) # Keep other institution NUMBER to iterate over in some plots
			if (bPlot.anonymize) {
				tmp.names$SCHOOL_NAME <- paste("School", as.numeric(as.factor(tmp.names$SCHOOL_NUMBER)))
			}
			setkey(tmp.data, SCHOOL_NUMBER, YEAR)
		}
		
		if ("INSTRUCTOR_NUMBER" %in% names(tmp.data)) { #Add both school and district number regardless of 
			if (!"INSTRUCTOR_NAME" %in% names(tmp.data)) {
				tmp.num <- seq(length(grep('INSTRUCTOR_NUMBER', names(sgp_object@Data))))
				eval(parse(text=paste("sgp_object@Data$INSTRUCTOR_NAME_", tmp.num,
					"<- paste('Instructor', as.factor(sgp_object@Data$INSTRUCTOR_NUMBER_", tmp.num, "))", sep="")))
			}
			tmp.names <- data.frame(sgp_object@Data[,c(grep('INSTRUCTOR_NUMBER', names(sgp_object@Data)), grep('INSTRUCTOR_NAME', names(sgp_object@Data)),
				grep('SCHOOL_NUMBER', names(sgp_object@Data)), grep('SCHOOL_NAME', names(sgp_object@Data)),
				grep('DISTRICT_NUMBER', names(sgp_object@Data)), grep('DISTRICT_NAME', names(sgp_object@Data))), with=FALSE])
			inst.id.index <- grep('INSTRUCTOR_NUMBER', names(tmp.names)); inst.name.index <- grep('INSTRUCTOR_NAME', names(tmp.names))
			sch.id.index <- grep('SCHOOL_NUMBER', names(tmp.names)); sch.name.index <- grep('SCHOOL_NAME', names(tmp.names))
			dst.id.index <- grep('DISTRICT_NUMBER', names(tmp.names)); dst.name.index <- grep('DISTRICT_NAME', names(tmp.names))
			tmp.names<- eval(parse(text=paste("unique(data.table(INSTRUCTOR_NUMBER=c(", paste("tmp.names[,", inst.id.index, "]", collapse=","),
				"), INSTRUCTOR_NAME=c(", paste("tmp.names[,", inst.name.index, "]", collapse=","),
				"), SCHOOL_NUMBER=rep(", paste("tmp.names[,", sch.id.index, "],", length(inst.id.index), collapse=","),
				"), SCHOOL_NAME=rep(", paste("tmp.names[,", sch.name.index, "],", length(inst.id.index), collapse=","), 
				"), DISTRICT_NUMBER=rep(", paste("tmp.names[,", dst.id.index, "],", length(inst.id.index), collapse=","),
				"), DISTRICT_NAME=rep(", paste("tmp.names[,", dst.name.index, "],", length(inst.id.index), collapse=","), ")))")))
			if (bPlot.anonymize) {
				tmp.names$INSTRUCTOR_NAME <- paste("Instructor", as.numeric(as.factor(tmp.names$INSTRUCTOR_NUMBER)))
				tmp.names$SCHOOL_NAME <- paste("School", as.numeric(as.factor(tmp.names$SCHOOL_NUMBER)))
				tmp.names$DISTRICT_NAME <- paste("District", as.numeric(as.factor(tmp.names$DISTRICT_NUMBER)))
			}

			if ("INSTRUCTOR_NUMBER" %in% names(tmp.data) & "SCHOOL_NUMBER" %in% names(tmp.data) & !"DISTRICT_NUMBER" %in% names(tmp.data)) {
				setkeyv(tmp.names, c("INSTRUCTOR_NUMBER", "SCHOOL_NUMBER"))
				setkeyv(tmp.data, c("INSTRUCTOR_NUMBER", "SCHOOL_NUMBER"))
			}
			if ("INSTRUCTOR_NUMBER" %in% names(tmp.data) & !"SCHOOL_NUMBER" %in% names(tmp.data) & "DISTRICT_NUMBER" %in% names(tmp.data)) {
				setkeyv(tmp.names, c("INSTRUCTOR_NUMBER", "DISTRICT_NUMBER"))
				setkeyv(tmp.data, c("INSTRUCTOR_NUMBER", "DISTRICT_NUMBER"))
			}
			tmp.names[tmp.data, mult="last"][!is.na(INSTRUCTOR_NUMBER)]
		} else tmp.names[tmp.data, mult="last"]
	}

	get.my.iters <- function(tmp.data, bubblePlot_LEVEL, ...) {
		my.iters <- list()

		# Year Stuff

		if (is.null(bPlot.years)) {
			if ("YEAR" %in% names(tmp.data)) {
				my.iters$tmp.years <- tail(sort(unique(tmp.data$YEAR)), 1)
			} else {
				my.iters$tmp.years <- "All Years"
			}
		} else {
			my.iters$tmp.years <- bPlot.years
		}

		# Content Area Stuff

		if (is.null(bPlot.content_areas)) {
			if ("CONTENT_AREA" %in% names(tmp.data)) {
				my.iters$tmp.content_areas <- unique(tmp.data$CONTENT_AREA) %w/o% NA
			} else {
				my.iters$tmp.content_areas <- "All Content Areas"
			}
		} else {
			my.iters$tmp.content_areas <- bPlot.content_areas
		}

		# Reconcile choice of District and Schools

		if (is.null(bPlot.instructors) & is.null(bPlot.schools) & is.null(bPlot.districts)) {
			if ("DISTRICT_NUMBER" %in% names(tmp.data)) {
				if (identical(my.iters$tmp.years, "All Years")) {
					my.iters$tmp.districts <- sort(unique(tmp.data$DISTRICT_NUMBER)) %w/o% NA
				} else {
					my.iters$tmp.districts <- sort(unique(tmp.data[YEAR %in% my.iters$tmp.years]$DISTRICT_NUMBER)) %w/o% NA
				}
			}
			if ("SCHOOL_NUMBER" %in% names(tmp.data)) {
				if (identical(my.iters$tmp.years, "All Years")) {
					my.iters$tmp.schools <- sort(unique(tmp.data$SCHOOL_NUMBER)) %w/o% NA
				} else {
					my.iters$tmp.schools <- sort(unique(tmp.data[YEAR %in% my.iters$tmp.years]$SCHOOL_NUMBER)) %w/o% NA
				}
			}
			if ("INSTRUCTOR_NUMBER" %in% names(tmp.data)) {
				if (identical(my.iters$tmp.years, "All Years")) {
					my.iters$tmp.instructors <- sort(unique(tmp.data$INSTRUCTOR_NUMBER)) %w/o% NA
				} else {
					my.iters$tmp.instructors <- sort(unique(tmp.data[YEAR %in% my.iters$tmp.years]$INSTRUCTOR_NUMBER)) %w/o% NA
				}
			}
		}

		if (is.null(bPlot.instructors) & is.null(bPlot.schools) & !is.null(bPlot.districts)) {
	 		my.iters$tmp.districts <- bPlot.districts
			if ("SCHOOL_NUMBER" %in% names(tmp.data)) my.iters$tmp.schools <- unique(tmp.data$SCHOOL_NUMBER[tmp.data$DISTRICT_NUMBER %in% my.iters$tmp.districts]) %w/o% NA
			if ("INSTRUCTOR_NUMBER" %in% names(tmp.data)) my.iters$tmp.instructors <- unique(tmp.data$INSTRUCTOR_NUMBER[tmp.data$INSTRUCTOR_NUMBER %in% my.iters$tmp.districts]) %w/o% NA
		}

		if (is.null(bPlot.instructors) & !is.null(bPlot.schools) & is.null(bPlot.districts)) {
	 		my.iters$tmp.schools <- bPlot.schools 
			if ("DISTRICT_NUMBER" %in% names(tmp.data)) my.iters$tmp.districts <- unique(tmp.data$DISTRICT_NUMBER[tmp.data$SCHOOL_NUMBER %in% my.iters$tmp.schools]) %w/o% NA
			if ("INSTRUCTOR_NUMBER" %in% names(tmp.data)) my.iters$tmp.instructors <- unique(tmp.data$INSTRUCTOR_NUMBER[tmp.data$INSTRUCTOR_NUMBER %in% my.iters$tmp.districts]) %w/o% NA
		}

		if (!is.null(bPlot.instructors) & is.null(bPlot.schools) & is.null(bPlot.districts)) {
			my.iters$tmp.instructors <- bPlot.instructors 
			if ("DISTRICT_NUMBER" %in% names(tmp.data)) my.iters$tmp.districts <- unique(tmp.data$DISTRICT_NUMBER[tmp.data$SCHOOL_NUMBER %in% my.iters$tmp.schools]) %w/o% NA
			if ("SCHOOL_NUMBER" %in% names(tmp.data)) my.iters$tmp.schools <- unique(tmp.data$SCHOOL_NUMBER[tmp.data$DISTRICT_NUMBER %in% my.iters$tmp.districts]) %w/o% NA
		}

		if (is.null(bPlot.instructors) & !is.null(bPlot.schools) & !is.null(bPlot.districts)) {
	 		my.iters$tmp.districts <- bPlot.districts
	 		my.iters$tmp.schools <- bPlot.schools 
			if ("INSTRUCTOR_NUMBER" %in% names(tmp.data)) my.iters$tmp.instructors <- unique(tmp.data$INSTRUCTOR_NUMBER[tmp.data$INSTRUCTOR_NUMBER %in% my.iters$tmp.districts]) %w/o% NA
		}

		if (!is.null(bPlot.instructors) & !is.null(bPlot.schools) & is.null(bPlot.districts)) {
	 		my.iters$tmp.schools <- bPlot.schools 
			my.iters$tmp.instructors <- bPlot.instructors 
			if ("DISTRICT_NUMBER" %in% names(tmp.data)) my.iters$tmp.districts <- unique(tmp.data$DISTRICT_NUMBER[tmp.data$SCHOOL_NUMBER %in% my.iters$tmp.schools]) %w/o% NA
		}

		if (!is.null(bPlot.instructors) & is.null(bPlot.schools) & !is.null(bPlot.districts)) {
	 		my.iters$tmp.districts <- bPlot.districts
			my.iters$tmp.instructors <- bPlot.instructors 
			if ("SCHOOL_NUMBER" %in% names(tmp.data)) my.iters$tmp.schools <- unique(tmp.data$SCHOOL_NUMBER[tmp.data$DISTRICT_NUMBER %in% my.iters$tmp.districts]) %w/o% NA
		}

		if (!is.null(bPlot.schools) & !is.null(bPlot.districts)) {
	 		my.iters$tmp.districts <- bPlot.districts
	 		my.iters$tmp.schools <- bPlot.schools 
	 		my.iters$tmp.instructors <- bPlot.instructors 
			my.iters$tmp.instructors <- unique(c(my.iters$tmp.instructors, tmp.data$INSTRUCTOR_NUMBER[tmp.data$SCHOOL_NUMBER %in% my.iters$tmp.schools])) %w/o% NA
			my.iters$tmp.schools <- unique(c(my.iters$tmp.schools, tmp.data$SCHOOL_NUMBER[tmp.data$DISTRICT_NUMBER %in% my.iters$tmp.districts])) %w/o% NA
			my.iters$tmp.districts <- unique(c(my.iters$tmp.districts, tmp.data$DISTRICT_NUMBER[tmp.data$SCHOOL_NUMBER %in% my.iters$tmp.schools])) %w/o% NA
		}

		# y.variable (include/not include prior achievement)

		if (bPlot.prior.achievement & length(grep("PERCENT_AT_ABOVE_PROFICIENT_PRIOR", names(tmp.data))) > 0) {
			if (bubblePlot_LEVEL=="Summary") my.iters$tmp.y.variable <- c("PERCENT_AT_ABOVE_PROFICIENT", "PERCENT_AT_ABOVE_PROFICIENT_PRIOR")
			if (bubblePlot_LEVEL=="Individual") my.iters$tmp.y.variable <- c("SCALE_SCORE", "SCALE_SCORE_PRIOR")
		} else {
			if (bubblePlot_LEVEL=="Summary") my.iters$tmp.y.variable <- "PERCENT_AT_ABOVE_PROFICIENT"
			if (bubblePlot_LEVEL=="Individual") my.iters$tmp.y.variable <- "SCALE_SCORE"
		}
		return(my.iters)
	} ## END get.my.iters

	get.my.level.labels <- function(bPlot.level.cuts) {
		tmp.list <- list()
		tmp.list[[1]] <- paste("Less than", bPlot.level.cuts[2], "percent")
		if (length(bPlot.level.cuts) > 3) {
			for (i in 2:(length(bPlot.level.cuts)-1)) {
				tmp.list[[i]] <- paste(bPlot.level.cuts[i], "to", bPlot.level.cuts[i+1], "percent")
			}
		}
		tmp.list[[length(bPlot.level.cuts)-1]] <- paste("More than", bPlot.level.cuts[length(bPlot.level.cuts)-1], "percent")
	do.call(c, tmp.list)
	}

	get.bPlot.data <- function(tmp.bPlot.data) {
		tmp <- "MEDIAN_SGP_COUNT >= bPlot.minimum.n"
		if (content_area.iter != "All Content Areas") tmp <- paste("CONTENT_AREA==as.character(content_area.iter) &", tmp)
		if (year.iter != "All Years") tmp <- paste("YEAR==year.iter &", tmp)
		subset(tmp.bPlot.data, eval(parse(text=tmp)))
	}


#################################################################################################################
####
#### Summary Level bubblePlots
####
#################################################################################################################

### < 100 are @Summary level bubblePlots

bubblePlot_LEVEL <- "Summary"
 

###################################################################
### BubblePlot Style 1 (State level bubblePlots by Schools)
###################################################################

if (1 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 1", date()))

		### Data sets and relevant quantities used for bubblePlots

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR__SCHOOL_ENROLLMENT_STATUS"]][
				SCHOOL_ENROLLMENT_STATUS=="Enrolled School: Yes"]
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR"]]
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) { ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL) 

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=NULL, 
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=NULL, 
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(50, 100, 250, 500),
			bubble_plot_labels.LEVELS=NULL, #levels(bubblePlot[["subset.factor"]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=bPlot.data[["SCHOOL_NAME"]],
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(state.name.label, "School Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="School Size",
			bubble_plot_titles.LEGEND2_P1=NULL,
			bubble_plot_titles.LEGEND2_P2=NULL,

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.04, 0.11),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
			bubble_plot_configs.BUBBLE_COLOR="deeppink2",
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(state.name.file.label, year.iter, capwords(content_area.iter), 
				"State", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "State", "Style_1"),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 1", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 1


#######################################################################################
### BubblePlot Style 2 (State level bubblePlots with district schools highlighted 
### by supplied bPlot.levels factor or without bPlot.levels factor
#######################################################################################

if (2 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 2", date()))

		### Data sets and relevant quantities used for bubblePlots

		if (!is.null(bPlot.levels)) {
			tmp.bPlot.levels.txt <- parse(text=paste("list(", paste("PCT_", bPlot.levels, "=(100*length(grep('Yes',", 
				bPlot.levels, "))/length(grep('Yes|No',", bPlot.levels, ")))", sep="", collapse=","), ")"))
			bPlot.levels <- as.list(bPlot.levels)
		}

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR__SCHOOL_ENROLLMENT_STATUS"]][
				SCHOOL_ENROLLMENT_STATUS=="Enrolled School: Yes"]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[SCHOOL_ENROLLMENT_STATUS=="Enrolled School: Yes", 
					eval(tmp.bPlot.levels.txt), by=list(SCHOOL_NUMBER, CONTENT_AREA, YEAR)]
				setkeyv(tmp.bPlot.data, c("SCHOOL_NUMBER", "CONTENT_AREA", "YEAR"))
				setkeyv(tmp.bPlot.levels.data, c("SCHOOL_NUMBER", "CONTENT_AREA", "YEAR"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				for (l in seq_along(bPlot.level.cuts)) {
					my.level.labels <- get.my.level.labels(bPlot.level.cuts[[l]])
					eval(parse(text=paste("tmp.bPlot.data$", paste("PCT_", bPlot.levels[[l]], sep=""), "<- cut(tmp.bPlot.data$PCT_", bPlot.levels[[l]], 
						", bPlot.level.cuts[[l]], include.lowest=TRUE, labels=my.level.labels, ordered_result=TRUE)", sep="")))
				}
			}
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR"]]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[, eval(tmp.bPlot.levels.txt), by=list(SCHOOL_NUMBER, CONTENT_AREA, YEAR)]
				setkeyv(tmp.bPlot.data, c("SCHOOL_NUMBER", "CONTENT_AREA", "YEAR"))
				setkeyv(tmp.bPlot.levels.data, c("SCHOOL_NUMBER", "CONTENT_AREA", "YEAR"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				for (l in seq_along(bPlot.level.cuts)) {
					my.level.labels <- get.my.level.labels(bPlot.level.cuts[[l]])
					eval(parse(text=paste("tmp.bPlot.data$", paste("PCT_", bPlot.levels[[l]], sep=""), "<- cut(tmp.bPlot.data$PCT_", bPlot.levels[[l]], 
						", bPlot.level.cuts[[l]], include.lowest=TRUE, labels=my.level.labels, ordered_result=TRUE)", sep="")))
				}
			}
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		if (is.null(bPlot.levels)) bPlot.levels <- list(A=NULL)

		for (bPlot.levels.iter in unlist(bPlot.levels)) {  ### Loop over bPlot.levels
		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)

		# Loop over current & prior and bPlot.levels 
		for (levels.iter in levels(factor(eval(parse(text=paste("bPlot.data$PCT_", bPlot.levels.iter, sep="")))))) {
		for (y.variable.iter in my.iters$tmp.y.variable) { 

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL) 

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]==levels.iter), 
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]], 
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(50, 100, 250, 500),
			bubble_plot_labels.LEVELS=levels(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=bPlot.data[["SCHOOL_NAME"]],
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(state.name.label, "School Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="School Size",
			bubble_plot_titles.LEGEND2_P1="Percentage Students",
			bubble_plot_titles.LEGEND2_P2=paste(sapply(head(unlist(strsplit(bPlot.levels.iter, "_")), -1), capwords), collapse=" "),

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.04, 0.11),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
			bubble_plot_configs.BUBBLE_COLOR=NULL,
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(state.name.file.label, year.iter, capwords(content_area.iter), capwords(levels.iter), "State", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "State", "Style_2", bPlot.levels.iter),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## END loop over levels.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter
		} ## End loop over bPlot.levels.iter

		message(paste("\tFinished bubblePlot Style 2", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 2


#######################################################################################
### BubblePlot Style 3 (State level bubblePlots with instructors highlighted 
### by supplied bPlot.levels factor or without bPlot.levels factor
#######################################################################################

if (3 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 3", date()))

		### Data sets and relevant quantities used for bubblePlots

		if (!is.null(bPlot.levels)) {
			tmp.bPlot.levels.txt <- parse(text=paste("100*length(grep('Yes',", bPlot.levels, "))/length(grep('Yes|No',", bPlot.levels, "))", sep=""))
		}

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["STATE"]][["STATE__INSTRUCTOR_NUMBER__INSTRUCTOR_ENROLLMENT_STATUS"]][INSTRUCTOR_ENROLLMENT_STATUS=="Enrolled Instructor: Yes"]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[INSTRUCTOR_ENROLLMENT_STATUS=="Enrolled Instructor: Yes", 
					eval(tmp.bPlot.levels.txt), by=list(INSTRUCTOR_NUMBER)]
				setkeyv(tmp.bPlot.data, c("INSTRUCTOR_NUMBER"))
				setkeyv(tmp.bPlot.levels.data, c("INSTRUCTOR_NUMBER"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				my.level.labels <- get.my.level.labels(bPlot.level.cuts)
				tmp.bPlot.data$V1 <- cut(tmp.bPlot.data$V1, bPlot.level.cuts, include.lowest=TRUE, labels=my.level.labels)
			}
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["STATE"]][["STATE__INSTRUCTOR_NUMBER"]]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[, eval(tmp.bPlot.levels.txt), by=list(INSTRUCTOR_NUMBER)]
				setkeyv(tmp.bPlot.data, c("INSTRUCTOR_NUMBER"))
				setkeyv(tmp.bPlot.levels.data, c("INSTRUCTOR_NUMBER"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				my.level.labels <- get.my.level.labels(bPlot.level.cuts)
				tmp.bPlot.data$V1 <- cut(tmp.bPlot.data$V1, bPlot.level.cuts, include.lowest=TRUE, labels=my.level.labels)
			}
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {
		for (content_area.iter in my.iters$tmp.content_areas) {

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)

		# Loop over current & prior and bPlot.levels 

		for (levels.iter in levels(factor(bPlot.data$V1))) {
		for (y.variable.iter in my.iters$tmp.y.variable) { 

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL) 

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]==levels.iter), 
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]], 
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=levels(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=bPlot.data[["INSTRUCTOR_NAME"]],
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(state.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Class Size",
			bubble_plot_titles.LEGEND2_P1="Percentage Students",
			bubble_plot_titles.LEGEND2_P2=paste(sapply(head(unlist(strsplit(bPlot.levels, "_")), -1), capwords), collapse=" "),

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.08),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
			bubble_plot_configs.BUBBLE_COLOR=NULL,
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(state.name.file.label, year.iter, capwords(content_area.iter), capwords(levels.iter), "State", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "State", "Style_3", bPlot.levels),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## END loop over levels.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 3", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 3


#######################################################################################
### BubblePlot Style 10 (State level bubblePlots with district schools highlighted)
#######################################################################################

if (10 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 10", date()))

		### Data sets and relevant quantities used for bubblePlots

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR__SCHOOL_ENROLLMENT_STATUS"]][
				SCHOOL_ENROLLMENT_STATUS=="Enrolled School: Yes"]
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR"]]
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)

		# Loop over unique districts

		for (district_number.iter in intersect(my.iters$tmp.districts, bPlot.data$DISTRICT_NUMBER)) { ### Loop over DISTRICT NUMBERS
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL)
		district.name.label <- as.character(bPlot.data[DISTRICT_NUMBER==district_number.iter]$DISTRICT_NAME[1])

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[["DISTRICT_NUMBER"]]==district_number.iter),
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=NULL,
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(50, 100, 250, 500),
			bubble_plot_labels.LEVELS=NULL, #levels(bubblePlot[["subset.factor"]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=bPlot.data[["SCHOOL_NAME"]],
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(district.name.label, "School Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="School Size",
			bubble_plot_titles.LEGEND2_P1=NULL,
			bubble_plot_titles.LEGEND2_P2=NULL,

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.04, 0.11),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.01,
			bubble_plot_configs.BUBBLE_COLOR="deeppink2",
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(district.name.label, year.iter, capwords(content_area.iter), bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "District", "Style_10", paste("District", district_number.iter, sep="_")),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 10", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 10


#######################################################################################
### BubblePlot Style 11 (State level bubblePlots with district schools highlighted 
### by supplied bPlot.levels factor or without bPlot.levels factor
#######################################################################################

	if (11 %in% bPlot.styles) {
		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 11", date()))

		### Data sets and relevant quantities used for bubblePlots

		if (!is.null(bPlot.levels)) {
			tmp.bPlot.levels.txt <- parse(text=paste("list(", paste("PCT_", bPlot.levels, "=(100*length(grep('Yes',", 
				bPlot.levels, "))/length(grep('Yes|No',", bPlot.levels, ")))", sep="", collapse=","), ")"))
			bPlot.levels <- as.list(bPlot.levels)
		} else bPlot.levels <- list(A=NULL)

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR__SCHOOL_ENROLLMENT_STATUS"]][
				SCHOOL_ENROLLMENT_STATUS=="Enrolled School: Yes"]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[SCHOOL_ENROLLMENT_STATUS=="Enrolled School: Yes", 
					eval(tmp.bPlot.levels.txt), by=list(SCHOOL_NUMBER, CONTENT_AREA, YEAR)]
				setkeyv(tmp.bPlot.data, c("SCHOOL_NUMBER", "CONTENT_AREA", "YEAR"))
				setkeyv(tmp.bPlot.levels.data, c("SCHOOL_NUMBER", "CONTENT_AREA", "YEAR"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				for (l in seq_along(bPlot.level.cuts)) {
					my.level.labels <- get.my.level.labels(bPlot.level.cuts[[l]])
					eval(parse(text=paste("tmp.bPlot.data$", paste("PCT_", bPlot.levels[[l]], sep=""), "<- cut(tmp.bPlot.data$PCT_", bPlot.levels[[l]], 
						", bPlot.level.cuts[[l]], include.lowest=TRUE, labels=my.level.labels, ordered_result=TRUE)", sep="")))
				}
			}
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR"]]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[, eval(tmp.bPlot.levels.txt), by=list(SCHOOL_NUMBER, CONTENT_AREA, YEAR)]
				setkeyv(tmp.bPlot.data, c("SCHOOL_NUMBER", "CONTENT_AREA", "YEAR"))
				setkeyv(tmp.bPlot.levels.data, c("SCHOOL_NUMBER", "CONTENT_AREA", "YEAR"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				for (l in seq_along(bPlot.level.cuts)) {
					my.level.labels <- get.my.level.labels(bPlot.level.cuts[[l]])
					eval(parse(text=paste("tmp.bPlot.data$", paste("PCT_", bPlot.levels[[l]], sep=""), "<- cut(tmp.bPlot.data$PCT_", bPlot.levels[[l]], 
						", bPlot.level.cuts[[l]], include.lowest=TRUE, labels=my.level.labels, ordered_result=TRUE)", sep="")))
				}
			}
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (bPlot.levels.iter in unlist(bPlot.levels)) {  ### Loop over bPlot.levels
		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas

		# Loop over unique districts

		for (district_number.iter in my.iters$tmp.districts) { ### Loop over DISTRICT NUMBERS

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)[DISTRICT_NUMBER==district_number.iter & !is.na(eval(parse(text=paste("PCT_", bPlot.levels.iter, sep=""))))]

		# Loop over current & prior and bPlot.levels 

		# for (levels.iter in levels(factor(eval(parse(text=paste("bPlot.data$PCT_", bPlot.levels.iter, sep="")))))) { # Only one district at a time, so show all schools in one flat file.
		for (y.variable.iter in my.iters$tmp.y.variable) { 

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL) 
		district.name.label <- as.character(bPlot.data[DISTRICT_NUMBER==district_number.iter]$DISTRICT_NAME[1])

		if (nrow(bPlot.data[which(bPlot.data[["DISTRICT_NUMBER"]]==district_number.iter),]) > 0) {

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[["DISTRICT_NUMBER"]]==district_number.iter), 
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]], 
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(50, 100, 250, 500),
			bubble_plot_labels.LEVELS=levels(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=bPlot.data[["SCHOOL_NAME"]],
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(district.name.label, "School Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="School Size",
			bubble_plot_titles.LEGEND2_P1="Percentage Students",
			bubble_plot_titles.LEGEND2_P2=paste(sapply(head(unlist(strsplit(bPlot.levels.iter, "_")), -1), capwords), collapse=" "),

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.04, 0.11),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
			bubble_plot_configs.BUBBLE_COLOR=NULL,
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(district.name.label, year.iter, capwords(content_area.iter), 
				capwords(bPlot.levels.iter), "Schools", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "District", "Style_11", paste("District", district_number.iter, sep="_"), bPlot.levels.iter),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)
			} # END if (...)
		} ## END loop over y.variable.iter
		# } ## END loop over levels.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter
		} ## End loop over bPlot.levels.iter

		message(paste("\tFinished bubblePlot Style 11", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 11


#######################################################################################
### BubblePlot Style 20 (State level bubblePlots with district teachers highlighted)
#######################################################################################

if (20 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 20", date()))

		### Data sets and relevant quantities used for bubblePlots

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["DISTRICT_NUMBER"]][["DISTRICT_NUMBER__INSTRUCTOR_NUMBER__INSTRUCTOR_ENROLLMENT_STATUS"]][
				INSTRUCTOR_ENROLLMENT_STATUS=="Enrolled Instructor: Yes"]
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["DISTRICT_NUMBER"]][["DISTRICT_NUMBER__INSTRUCTOR_NUMBER"]]
		}

		# Merge in teacher and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)

		# Loop over unique districts

		for (district_number.iter in intersect(my.iters$tmp.districts, bPlot.data$DISTRICT_NUMBER)) { ### Loop over DISTRICT NUMBERS
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL)
		district.name.label <- as.character(bPlot.data[DISTRICT_NUMBER==district_number.iter]$DISTRICT_NAME[1])

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[["DISTRICT_NUMBER"]]==district_number.iter),
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=NULL,
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=NULL, #levels(bubblePlot[["subset.factor"]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=bPlot.data[["INSTRUCTOR_NAME"]],
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(district.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Classroom Size",
			bubble_plot_titles.LEGEND2_P1=NULL,
			bubble_plot_titles.LEGEND2_P2=NULL,

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.08),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.01,
			bubble_plot_configs.BUBBLE_COLOR="deeppink2",
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(district.name.label, year.iter, capwords(content_area.iter), "District", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "District", "Style_20"),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 20", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 20


#######################################################################################
### BubblePlot Style 21 (State level bubblePlots with district schools highlighted 
### by supplied bPlot.levels factor or without bPlot.levels factor (one plot for showing all levels)
#######################################################################################

if (21 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 21", date()))

		### Data sets and relevant quantities used for bubblePlots

		if (!is.null(bPlot.levels)) {
			tmp.bPlot.levels.txt <- parse(text=paste("100*length(grep('Yes',", bPlot.levels, "))/length(grep('Yes|No',", bPlot.levels, "))", sep=""))
		}

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["DISTRICT_NUMBER"]][["DISTRICT_NUMBER__INSTRUCTOR_NUMBER"]][INSTRUCTOR_ENROLLMENT_STATUS=="Enrolled Instructor: Yes"]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[INSTRUCTOR_ENROLLMENT_STATUS=="Enrolled Instructor: Yes",
					eval(tmp.bPlot.levels.txt), by=list(DISTRICT_NUMBER, INSTRUCTOR_NUMBER)]
				setkeyv(tmp.bPlot.data, c("DISTRICT_NUMBER", "INSTRUCTOR_NUMBER"))
				setkeyv(tmp.bPlot.levels.data, c("DISTRICT_NUMBER", "INSTRUCTOR_NUMBER"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				my.level.labels <- get.my.level.labels(bPlot.level.cuts)
				tmp.bPlot.data$V1 <- cut(tmp.bPlot.data$V1, bPlot.level.cuts, include.lowest=TRUE, labels=my.level.labels)
			}
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["DISTRICT_NUMBER"]][["DISTRICT_NUMBER__INSTRUCTOR_NUMBER"]]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[, eval(tmp.bPlot.levels.txt), by=list(DISTRICT_NUMBER, INSTRUCTOR_NUMBER)]
				setkeyv(tmp.bPlot.data, c("DISTRICT_NUMBER", "INSTRUCTOR_NUMBER"))
				setkeyv(tmp.bPlot.levels.data, c("DISTRICT_NUMBER", "INSTRUCTOR_NUMBER"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				my.level.labels <- get.my.level.labels(bPlot.level.cuts)
				tmp.bPlot.data$V1 <- cut(tmp.bPlot.data$V1, bPlot.level.cuts, include.lowest=TRUE, labels=my.level.labels)
			}
		}


		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)

		# Loop over unique districts

		for (district_number.iter in intersect(my.iters$tmp.districts, bPlot.data$DISTRICT_NUMBER)) { ### Loop over DISTRICT NUMBERS
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL)
		district.name.label <- as.character(bPlot.data[DISTRICT_NUMBER==district_number.iter]$DISTRICT_NAME[1])

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[["DISTRICT_NUMBER"]]==district_number.iter), 
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]], 
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=levels(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=bPlot.data[["INSTRUCTOR_NAME"]],
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(district.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Classroom Size",
			bubble_plot_titles.LEGEND2_P1="Percentage Students",
			bubble_plot_titles.LEGEND2_P2=paste(sapply(head(unlist(strsplit(bPlot.levels, "_")), -1), capwords), collapse=" "),

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.08),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
			bubble_plot_configs.BUBBLE_COLOR=NULL,
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(district.name.label, year.iter, capwords(content_area.iter), "District", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "District", "Style_21", bPlot.levels),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 21", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 21



#######################################################################################
### BubblePlot Style 22 (State level bubblePlots with district schools highlighted 
### by supplied bPlot.levels factor or without bPlot.levels factor (one plot for each level)
#######################################################################################

if (22 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 22", date()))

		### Data sets and relevant quantities used for bubblePlots

		if (!is.null(bPlot.levels)) {
			tmp.bPlot.levels.txt <- parse(text=paste("100*length(grep('Yes',", bPlot.levels, "))/length(grep('Yes|No',", bPlot.levels, "))", sep=""))
		}

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["DISTRICT_NUMBER"]][["DISTRICT_NUMBER__INSTRUCTOR_NUMBER"]][INSTRUCTOR_ENROLLMENT_STATUS=="Enrolled Instructor: Yes"]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[INSTRUCTOR_ENROLLMENT_STATUS=="Enrolled Instructor: Yes",
					eval(tmp.bPlot.levels.txt), by=list(DISTRICT_NUMBER, INSTRUCTOR_NUMBER)]
				setkeyv(tmp.bPlot.data, c("DISTRICT_NUMBER", "INSTRUCTOR_NUMBER"))
				setkeyv(tmp.bPlot.levels.data, c("DISTRICT_NUMBER", "INSTRUCTOR_NUMBER"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				my.level.labels <- get.my.level.labels(bPlot.level.cuts)
				tmp.bPlot.data$V1 <- cut(tmp.bPlot.data$V1, bPlot.level.cuts, include.lowest=TRUE, labels=my.level.labels)
			}
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["DISTRICT_NUMBER"]][["DISTRICT_NUMBER__INSTRUCTOR_NUMBER"]]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- sgp_object@Data[, eval(tmp.bPlot.levels.txt), by=list(DISTRICT_NUMBER, INSTRUCTOR_NUMBER)]
				setkeyv(tmp.bPlot.data, c("DISTRICT_NUMBER", "INSTRUCTOR_NUMBER"))
				setkeyv(tmp.bPlot.levels.data, c("DISTRICT_NUMBER", "INSTRUCTOR_NUMBER"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				my.level.labels <- get.my.level.labels(bPlot.level.cuts)
				tmp.bPlot.data$V1 <- cut(tmp.bPlot.data$V1, bPlot.level.cuts, include.lowest=TRUE, labels=my.level.labels)
			}
		}


		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas

		# Loop over unique districts

		for (district_number.iter in my.iters$tmp.districts) { ### Loop over DISTRICT NUMBERS

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)

		# Loop over current & prior and bPlot.levels 

		for (levels.iter in levels(factor(bPlot.data$V1))) {
		for (y.variable.iter in my.iters$tmp.y.variable) { 

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL) 
		district.name.label <- as.character(bPlot.data[DISTRICT_NUMBER==district_number.iter]$DISTRICT_NAME[1])

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]==levels.iter), 
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]], 
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=levels(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=bPlot.data[["INSTRUCTOR_NAME"]],
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(district.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Classroom Size",
			bubble_plot_titles.LEGEND2_P1="Percentage Students",
			bubble_plot_titles.LEGEND2_P2=paste(sapply(head(unlist(strsplit(bPlot.levels, "_")), -1), capwords), collapse=" "),

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.08),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
			bubble_plot_configs.BUBBLE_COLOR=NULL,
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(district.name.label, year.iter, capwords(content_area.iter), capwords(levels.iter), "District", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "District", "Style_22", bPlot.levels),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## END loop over levels.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 22", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 22


###################################################################
### BubblePlot Style 50 (School level bubblePlots by Instructors)
###################################################################

	if (50 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 50", date()), "\n")

		### Data sets and relevant quantities used for bubblePlots

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR__INSTRUCTOR_ENROLLMENT_STATUS"]][!is.na(INSTRUCTOR_ENROLLMENT_STATUS)]
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR"]]
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (district_number.iter in my.iters$tmp.districts) { ### Loop over DISTRICT NUMBERS
		for (school_number.iter in my.iters$tmp.schools) {  ### Loop over schools
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Subset data

		bPlot.data <- tmp.bPlot.data[YEAR==year.iter & SCHOOL_NUMBER==school_number.iter & CONTENT_AREA==content_area.iter & MEDIAN_SGP_COUNT >= bPlot.minimum.n]

		# Create labels and file path

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL) 
		school.name.label <- as.character(bPlot.data[SCHOOL_NUMBER== school_number.iter]$SCHOOL_NAME[1])

		### Create bubblePlot ###
		if (dim(bPlot.data)[1] > 0) { # some institutions don't teach all content areas (elementary math)

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=NULL, 
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=NULL, 
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=NULL, #levels(bubblePlot[["subset.factor"]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=paste(bPlot.data[["INSTRUCTOR_NAME"]], "-", bPlot.data[["SCHOOL_NAME"]]),
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(state.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Number of Students",
			bubble_plot_titles.LEGEND2_P1=NULL,
			bubble_plot_titles.LEGEND2_P2=NULL,

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.14),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
			bubble_plot_configs.BUBBLE_COLOR="deeppink2",
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(school.name.label, year.iter, capwords(content_area.iter), 
				"Instructors", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "Instructor", "Style_50", paste("District", district_number.iter, sep="_")),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)
		} ## END if 
		} ## END loop over y.variable.iter
		} ## End loop over content_area.iter
		} ## End loop over district_number.iter
		} ## End loop over school_number.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 50", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END bubblePlot style 50


#######################################################################################
### BubblePlot Style 53 (STATE level --or multiple districts-- bubblePlots with District Instructors highlighted)
#######################################################################################

	if (53 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 53", date()), "\n")

		### Data sets and relevant quantities used for bubblePlots

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR__INSTRUCTOR_ENROLLMENT_STATUS"]][!is.na(INSTRUCTOR_ENROLLMENT_STATUS)]
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR"]]
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas

		# Subset data

		bPlot.data <- tmp.bPlot.data[YEAR==year.iter & CONTENT_AREA==content_area.iter & DISTRICT_NUMBER %in% my.iters$tmp.districts & MEDIAN_SGP_COUNT >= bPlot.minimum.n]

		# Loop over unique schools IN DISTRICT ONLY
		for (district_number.iter in intersect(my.iters$tmp.districts, bPlot.data$DISTRICT_NUMBER)) { ### Loop over DISTRICT NUMBERS   
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Create labels and file path

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL)
		district.name.label <- as.character(bPlot.data[DISTRICT_NUMBER==district_number.iter]$DISTRICT_NAME[1])

		### Create bubblePlot ###
		
		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[["DISTRICT_NUMBER"]]==district_number.iter),
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=NULL,
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=NULL, #levels(bubblePlot[["subset.factor"]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=paste(bPlot.data[["INSTRUCTOR_NAME"]], "-", bPlot.data[["SCHOOL_NAME"]]),
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(district.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Number of Students",
			bubble_plot_titles.LEGEND2_P1=NULL,
			bubble_plot_titles.LEGEND2_P2=NULL,

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.14),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.01,
			bubble_plot_configs.BUBBLE_COLOR="deeppink2",
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(district.name.label, year.iter, 
				capwords(content_area.iter), "Instructors", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "Instructor", "Style_53", paste("District", district_number.iter, sep="_")),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 53", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END bubblePlot style 53


#######################################################################################
### BubblePlot Style 55 (STATE level --or multiple districts-- bubblePlots with District Instructors highlighted
### by supplied bPlot.levels factor or without bPlot.levels factor)
#######################################################################################

	if (55 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 55", date()), "\n")

		### Data sets and relevant quantities used for bubblePlots

		if (!is.null(bPlot.levels)) {
			tmp.bPlot.levels.txt <- parse(text=paste("list(", paste("PCT_", bPlot.levels, "=(100*length(grep('Yes',", 
				bPlot.levels, "))/length(grep('Yes|No',", bPlot.levels, ")))", sep="", collapse=","), ")"))
			bPlot.levels <- as.list(bPlot.levels)
		} else bPlot.levels <- list(A=NULL)

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR__INSTRUCTOR_ENROLLMENT_STATUS"]][!is.na(INSTRUCTOR_ENROLLMENT_STATUS)]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- data.frame(eval(parse(text=paste("sgp_object@Data[,c(grep('INSTRUCTOR_NUMBER', names(sgp_object@Data)), 
					intersect(grep('ENROLLMENT_STATUS', names(sgp_object@Data)), grep('INSTRUCTOR', names(sgp_object@Data))),", 
					"grep('CONTENT_AREA', names(sgp_object@Data)), grep('YEAR', names(sgp_object@Data)),",
					paste("grep('", bPlot.levels, "', names(sgp_object@Data))", sep="", collapse=","), "), with=FALSE]"))))
				inst.index <- grep('INSTRUCTOR_NUMBER', names(tmp.bPlot.levels.data))
				enroll.index <- grep('ENROLLMENT_STATUS', names(tmp.bPlot.levels.data))
				tmp.bPlot.levels.data<- eval(parse(text=paste("data.table(INSTRUCTOR_NUMBER=c(", paste("tmp.bPlot.levels.data[,", inst.index, "]", collapse=","),
					"), ENROLLMENT_STATUS=c(", paste("as.character(tmp.bPlot.levels.data[,", enroll.index, "])", collapse=","),
					"), CONTENT_AREA=rep(tmp.bPlot.levels.data[, 'CONTENT_AREA'],", length(inst.index),
					"), YEAR=rep(tmp.bPlot.levels.data[, 'YEAR'],", length(inst.index),
					"), ", paste(bPlot.levels, "=rep(tmp.bPlot.levels.data[, '", bPlot.levels, "'],", length(inst.index), ")", sep="", collapse=","),
					", key=c('INSTRUCTOR_NUMBER', 'ENROLLMENT_STATUS'))[!is.na(INSTRUCTOR_NUMBER) & ENROLLMENT_STATUS=='Enrolled Instructor: Yes']")))
				tmp.bPlot.levels.data <- tmp.bPlot.levels.data[, eval(tmp.bPlot.levels.txt), by=list(INSTRUCTOR_NUMBER, CONTENT_AREA, YEAR)]
				setkeyv(tmp.bPlot.data, c("INSTRUCTOR_NUMBER", "CONTENT_AREA", "YEAR"))
				setkeyv(tmp.bPlot.levels.data, c("INSTRUCTOR_NUMBER", "CONTENT_AREA", "YEAR"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				for (l in seq_along(bPlot.level.cuts)) {
					my.level.labels <- get.my.level.labels(bPlot.level.cuts[[l]])
					eval(parse(text=paste("tmp.bPlot.data$", paste("PCT_", bPlot.levels[[l]], sep=""), "<- cut(tmp.bPlot.data$PCT_", bPlot.levels[[l]], 
						", bPlot.level.cuts[[l]], include.lowest=TRUE, labels=my.level.labels, ordered_result=TRUE)", sep="")))
				}
			}
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR"]]
			if (!is.null(bPlot.levels)) {

				tmp.bPlot.levels.data <- data.frame(eval(parse(text=paste("sgp_object@Data[,c(grep('INSTRUCTOR_NUMBER', names(sgp_object@Data)), 
					grep('CONTENT_AREA', names(sgp_object@Data)), grep('YEAR', names(sgp_object@Data)),",
					paste("grep('", bPlot.levels, "', names(sgp_object@Data))", sep="", collapse=","), "), with=FALSE]"))))
				inst.index <- grep('INSTRUCTOR_NUMBER', names(tmp.bPlot.levels.data))
				tmp.bPlot.levels.data<- eval(parse(text=paste("data.table(INSTRUCTOR_NUMBER=c(", paste("tmp.bPlot.levels.data[,", inst.index, "]", collapse=","),
					"), CONTENT_AREA=rep(tmp.bPlot.levels.data[, 'CONTENT_AREA'],", length(inst.index),
					"), YEAR=rep(tmp.bPlot.levels.data[, 'YEAR'],", length(inst.index),
					"), ", paste(bPlot.levels, "=rep(tmp.bPlot.levels.data[, '", bPlot.levels, "'],", length(inst.index), ")", sep="", collapse=","),
					", key='INSTRUCTOR_NUMBER')[!is.na(INSTRUCTOR_NUMBER)]")))
				tmp.bPlot.levels.data <- tmp.bPlot.levels.data[, eval(tmp.bPlot.levels.txt), by=list(INSTRUCTOR_NUMBER, CONTENT_AREA, YEAR)]

				setkeyv(tmp.bPlot.data, c("INSTRUCTOR_NUMBER", "CONTENT_AREA", "YEAR"))
				setkeyv(tmp.bPlot.levels.data, c("INSTRUCTOR_NUMBER", "CONTENT_AREA", "YEAR"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				for (l in seq_along(bPlot.level.cuts)) {
					my.level.labels <- get.my.level.labels(bPlot.level.cuts[[l]])
					eval(parse(text=paste("tmp.bPlot.data$", paste("PCT_", bPlot.levels[[l]], sep=""), "<- cut(tmp.bPlot.data$PCT_", bPlot.levels[[l]], 
						", bPlot.level.cuts[[l]], include.lowest=TRUE, labels=my.level.labels, ordered_result=TRUE)", sep="")))
				}
			}
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data[!is.na(INSTRUCTOR_NUMBER)], bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (bPlot.levels.iter in unlist(bPlot.levels)) {  ### Loop over bPlot.levels
		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)[DISTRICT_NUMBER %in% my.iters$tmp.districts]

		# Loop over unique schools IN DISTRICT ONLY
		for (district_number.iter in intersect(my.iters$tmp.districts, bPlot.data$DISTRICT_NUMBER)) { ### Loop over DISTRICT NUMBERS   
		for (levels.iter in levels(factor(eval(parse(text=paste("bPlot.data$PCT_", bPlot.levels.iter, sep="")))))) {
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Create labels and file path

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL)
		district.name.label <- as.character(bPlot.data[DISTRICT_NUMBER==district_number.iter]$DISTRICT_NAME[1])

		### Create bubblePlot ###

		if (nrow(bPlot.data[which(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]==levels.iter & bPlot.data[["DISTRICT_NUMBER"]]==district_number.iter),]) > 0) {

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]==levels.iter & bPlot.data[["DISTRICT_NUMBER"]]==district_number.iter),
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]],
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=levels(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=paste(bPlot.data[["INSTRUCTOR_NAME"]], "-", bPlot.data[["SCHOOL_NAME"]]),
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(district.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Number of Students",
			bubble_plot_titles.LEGEND2_P1="Percentage Students",
			bubble_plot_titles.LEGEND2_P2=paste(sapply(head(unlist(strsplit(bPlot.levels.iter, "_")), -1), capwords), collapse=" "),

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.14),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.01,
			bubble_plot_configs.BUBBLE_COLOR=NULL,
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(district.name.label, year.iter, capwords(content_area.iter),
				capwords(levels.iter), "Classrooms", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "Instructor", "Style_55", 
				paste("District", district_number.iter, sep="_"), bPlot.levels.iter),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)
			} # END if (...)
		} ## END loop over y.variable.iter
		} ## END loop over levels.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter
		} ## End loop over bPlot.levels.iter

		message(paste("\tFinished bubblePlot Style 55", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END bubblePlot style 55


#######################################################################################
### BubblePlot Style 57 (District level bubblePlots with Schools' Instructors highlighted)
#######################################################################################


	if (57 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 57", date()), "\n")

		### Data sets and relevant quantities used for bubblePlots

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR__INSTRUCTOR_ENROLLMENT_STATUS"]][!is.na(INSTRUCTOR_ENROLLMENT_STATUS)]
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR"]]
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data, bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas
		for (district_number.iter in my.iters$tmp.districts) {  ### Loop over districts

		# Subset data

		bPlot.data <- tmp.bPlot.data[YEAR==year.iter & CONTENT_AREA==content_area.iter & DISTRICT_NUMBER==district_number.iter & MEDIAN_SGP_COUNT >= bPlot.minimum.n]

		# Loop over unique schools IN DISTRICT ONLY
		for (school_number.iter in intersect(my.iters$tmp.schools, bPlot.data$SCHOOL_NUMBER)) {  ### Loop over schools
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Create labels and file path

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL)
		school.name.label <- as.character(bPlot.data[SCHOOL_NUMBER== school_number.iter]$SCHOOL_NAME[1])

		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[["SCHOOL_NUMBER"]]==school_number.iter),
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=NULL,
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=NULL, #levels(bubblePlot[["subset.factor"]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=paste(bPlot.data[["INSTRUCTOR_NAME"]], "-", bPlot.data[["SCHOOL_NAME"]]),
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(school.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Number of Students",
			bubble_plot_titles.LEGEND2_P1=NULL,
			bubble_plot_titles.LEGEND2_P2=NULL,

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.14),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.01,
			bubble_plot_configs.BUBBLE_COLOR="deeppink2",
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(school.name.label, year.iter, capwords(content_area.iter), 
				"Instructor", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "Instructor", "Style_57", paste("District", district_number.iter, sep="_")),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END loop over y.variable.iter
		} ## End loop over school_number.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter

		message(paste("\tFinished bubblePlot Style 57", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END bubblePlot style 57


#######################################################################################
### BubblePlot Style 59 (District level bubblePlots with Schools' Instructors highlighted, Demographic concentrations distinguished)
#######################################################################################

	if (59 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 59", date()), "\n")

		### Data sets and relevant quantities used for bubblePlots

		if (!is.null(bPlot.levels)) {
			tmp.bPlot.levels.txt <- parse(text=paste("list(", paste("PCT_", bPlot.levels, "=(100*length(grep('Yes',", 
				bPlot.levels, "))/length(grep('Yes|No',", bPlot.levels, ")))", sep="", collapse=","), ")"))
			bPlot.levels <- as.list(bPlot.levels)
		} else bPlot.levels <- list(A=NULL)

		if (bPlot.full.academic.year) {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR__INSTRUCTOR_ENROLLMENT_STATUS"]][!is.na(INSTRUCTOR_ENROLLMENT_STATUS)]
			if (!is.null(bPlot.levels)) {
				tmp.bPlot.levels.data <- data.frame(eval(parse(text=paste("sgp_object@Data[,c(grep('INSTRUCTOR_NUMBER', names(sgp_object@Data)), 
					intersect(grep('ENROLLMENT_STATUS', names(sgp_object@Data)), grep('INSTRUCTOR', names(sgp_object@Data))),", 
					"grep('CONTENT_AREA', names(sgp_object@Data)), grep('YEAR', names(sgp_object@Data)),",
					paste("grep('", bPlot.levels, "', names(sgp_object@Data))", sep="", collapse=","), "), with=FALSE]"))))
				inst.index <- grep('INSTRUCTOR_NUMBER', names(tmp.bPlot.levels.data))
				enroll.index <- grep('ENROLLMENT_STATUS', names(tmp.bPlot.levels.data))
				tmp.bPlot.levels.data<- eval(parse(text=paste("data.table(INSTRUCTOR_NUMBER=c(", paste("tmp.bPlot.levels.data[,", inst.index, "]", collapse=","),
					"), ENROLLMENT_STATUS=c(", paste("as.character(tmp.bPlot.levels.data[,", enroll.index, "])", collapse=","),
					"), CONTENT_AREA=rep(tmp.bPlot.levels.data[, 'CONTENT_AREA'],", length(inst.index),
					"), YEAR=rep(tmp.bPlot.levels.data[, 'YEAR'],", length(inst.index),
					"), ", paste(bPlot.levels, "=rep(tmp.bPlot.levels.data[, '", bPlot.levels, "'],", length(inst.index), ")", sep="", collapse=","),
					", key=c('INSTRUCTOR_NUMBER', 'ENROLLMENT_STATUS'))[!is.na(INSTRUCTOR_NUMBER) & ENROLLMENT_STATUS=='Enrolled Instructor: Yes']")))
				tmp.bPlot.levels.data <- tmp.bPlot.levels.data[, eval(tmp.bPlot.levels.txt), by=list(INSTRUCTOR_NUMBER, CONTENT_AREA, YEAR)]
				setkeyv(tmp.bPlot.data, c("INSTRUCTOR_NUMBER", "CONTENT_AREA", "YEAR"))
				setkeyv(tmp.bPlot.levels.data, c("INSTRUCTOR_NUMBER", "CONTENT_AREA", "YEAR"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				for (l in seq_along(bPlot.level.cuts)) {
					my.level.labels <- get.my.level.labels(bPlot.level.cuts[[l]])
					eval(parse(text=paste("tmp.bPlot.data$", paste("PCT_", bPlot.levels[[l]], sep=""), "<- cut(tmp.bPlot.data$PCT_", bPlot.levels[[l]], 
						", bPlot.level.cuts[[l]], include.lowest=TRUE, labels=my.level.labels, ordered_result=TRUE)", sep="")))
				}
			}
		} else {
			tmp.bPlot.data <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR"]]
			if (!is.null(bPlot.levels)) {

				tmp.bPlot.levels.data <- data.frame(eval(parse(text=paste("sgp_object@Data[,c(grep('INSTRUCTOR_NUMBER', names(sgp_object@Data)), 
					grep('CONTENT_AREA', names(sgp_object@Data)), grep('YEAR', names(sgp_object@Data)),",
					paste("grep('", bPlot.levels, "', names(sgp_object@Data))", sep="", collapse=","), "), with=FALSE]"))))
				inst.index <- grep('INSTRUCTOR_NUMBER', names(tmp.bPlot.levels.data))
				tmp.bPlot.levels.data<- eval(parse(text=paste("data.table(INSTRUCTOR_NUMBER=c(", paste("tmp.bPlot.levels.data[,", inst.index, "]", collapse=","),
					"), CONTENT_AREA=rep(tmp.bPlot.levels.data[, 'CONTENT_AREA'],", length(inst.index),
					"), YEAR=rep(tmp.bPlot.levels.data[, 'YEAR'],", length(inst.index),
					"), ", paste(bPlot.levels, "=rep(tmp.bPlot.levels.data[, '", bPlot.levels, "'],", length(inst.index), ")", sep="", collapse=","),
					", key='INSTRUCTOR_NUMBER')[!is.na(INSTRUCTOR_NUMBER)]")))
				tmp.bPlot.levels.data <- tmp.bPlot.levels.data[, eval(tmp.bPlot.levels.txt), by=list(INSTRUCTOR_NUMBER, CONTENT_AREA, YEAR)]

				setkeyv(tmp.bPlot.data, c("INSTRUCTOR_NUMBER", "CONTENT_AREA", "YEAR"))
				setkeyv(tmp.bPlot.levels.data, c("INSTRUCTOR_NUMBER", "CONTENT_AREA", "YEAR"))
				tmp.bPlot.data <- tmp.bPlot.levels.data[tmp.bPlot.data]
				for (l in seq_along(bPlot.level.cuts)) {
					my.level.labels <- get.my.level.labels(bPlot.level.cuts[[l]])
					eval(parse(text=paste("tmp.bPlot.data$", paste("PCT_", bPlot.levels[[l]], sep=""), "<- cut(tmp.bPlot.data$PCT_", bPlot.levels[[l]], 
						", bPlot.level.cuts[[l]], include.lowest=TRUE, labels=my.level.labels, ordered_result=TRUE)", sep="")))
				}
			}
		}

		# Merge in school and district names and anonymize school names (if requested)

		tmp.bPlot.data <- names.merge(tmp.bPlot.data[!is.na(INSTRUCTOR_NUMBER)], bPlot.anonymize)

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(tmp.bPlot.data, bubblePlot_LEVEL)

		### Start loops for bubblePlots

		for (bPlot.levels.iter in unlist(bPlot.levels)) {  ### Loop over bPlot.levels
		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas
		for (district_number.iter in my.iters$tmp.districts) {  ### Loop over districts

		# Subset data

		bPlot.data <- get.bPlot.data(tmp.bPlot.data)[DISTRICT_NUMBER==district_number.iter]

		# Loop over unique schools IN DISTRICT ONLY
		for (school_number.iter in intersect(my.iters$tmp.schools, bPlot.data$SCHOOL_NUMBER)) {  ### Loop over schools
		for (levels.iter in levels(factor(eval(parse(text=paste("bPlot.data$PCT_", bPlot.levels.iter, sep="")))))) {
		for (y.variable.iter in my.iters$tmp.y.variable) {  ### Loop over CURRENT and PRIOR achievement (if requested)

		# Create labels and file path

		bPlot.labels <- create.bPlot.labels(year.iter, y.variable.iter, bubblePlot_LEVEL)
		school.name.label <- as.character(bPlot.data[SCHOOL_NUMBER== school_number.iter]$SCHOOL_NAME[1])

		### Create bubblePlot ###
			
		bubblePlot(
			bubble_plot_data.X=bPlot.data[["MEDIAN_SGP"]],
			bubble_plot_data.Y=bPlot.data[[y.variable.iter]],
			bubble_plot_data.SUBSET=which(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]==levels.iter), 
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=bPlot.data[["MEDIAN_SGP_COUNT"]],
			bubble_plot_data.LEVELS=bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]],
			bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[["MEDIAN_SGP"]], " (", bPlot.data[["MEDIAN_SGP_COUNT"]], ")", sep=""),
				paste(bPlot.data[[y.variable.iter]], " (", bPlot.data[[paste(y.variable.iter, "_COUNT", sep="")]], ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Median Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
			bubble_plot_labels.SIZE=c(10, 25, 50, 100),
			bubble_plot_labels.LEVELS=levels(bPlot.data[[paste("PCT_", bPlot.levels.iter, sep="")]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Median SGP (Count)"),
				paste(bPlot.labels$y.year.label, " (Count)")),
			bubble_plot_labels.BUBBLE_TITLES=paste(bPlot.data[["INSTRUCTOR_NAME"]], "-", bPlot.data[["SCHOOL_NAME"]]),
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(school.name.label, "Classroom Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, test.abbreviation.label, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="Number of Students",
			bubble_plot_titles.LEGEND2_P1="Percentage Students",
			bubble_plot_titles.LEGEND2_P2=paste(sapply(head(unlist(strsplit(bPlot.levels.iter, "_")), -1), capwords), collapse=" "),

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.14),
			bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.01,
			bubble_plot_configs.BUBBLE_COLOR=NULL,
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(school.name.label, year.iter, capwords(content_area.iter), 
				capwords(levels.iter), capwords(bPlot.levels.iter), "Classrooms", bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "Instructor", "Style_59", 
				paste("District", district_number.iter, sep="_"), bPlot.levels.iter),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)
		} ## END loop over y.variable.iter
		} ## END loop over levels.iter
		} ## End loop over school_number.iter
		} ## End loop over district_number.iter
		} ## End loop over content_area.iter
		} ## End loop over year.iter
		} ## End loop over bPlot.levels.iter

		message(paste("\tFinished bubblePlot Style 59", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ## END bubblePlot style 59


#################################################################################################################
#################################################################################################################
####
#### Individual Level bubblePlots
####
#################################################################################################################
#################################################################################################################

### >= 100 are @Data level bubblePlots

	if (any(bPlot.styles>=100)) {
		bubblePlot_LEVEL <- "Individual"

		if (any(c(150, 153) %in% bPlot.styles)) {
			if (!bPlot.prior.achievement) {
				message("bPlot.prior.achievement must be set to TRUE to use bubblePlot Style 150 &/or 153.  Setting bPlot.prior.achievement = TRUE internally")
				bPlot.prior.achievement <- TRUE
			}
		}

		### Utility functions

		get.my.cutscore.year <- function(state, content_area, year) {
			tmp.cutscore.years <- sapply(strsplit(names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]])[
				grep(content_area, names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]]))], "[.]"), function(x) x[2])
			if (any(!is.na(tmp.cutscore.years))) {
				if (year %in% tmp.cutscore.years) {
					return(paste(content_area, year, sep="."))
				} else {
					if (year==sort(c(year, tmp.cutscore.years))[1]) {
						return(content_area)
					} else {
						return(paste(content_area, sort(tmp.cutscore.years)[which(year==sort(c(year, tmp.cutscore.years)))-1], sep="."))
					}
				}
			} else {
				return(content_area)
			}
		}

		### Get tmp.years, tmp.content_areas, and tmp.y.variable

		my.iters <- get.my.iters(sgp_object@Data["VALID_CASE"], bubblePlot_LEVEL)


		### Copy @Data

		slot.data <- copy(sgp_object@Data)

		### Create PRIOR Scale Score, SGP, SGP_TARGET and CONTENT_AREA

		if (bPlot.prior.achievement) {
			if (!all(c("SCALE_SCORE_PRIOR", "SGP_PRIOR", "SGP_TARGET_PRIOR", "CONTENT_AREA_PRIOR") %in% names(slot.data))) {
				slot.data[, YEAR_INTEGER_TMP:=as.integer(as.factor(slot.data$YEAR))]
				setkeyv(slot.data, c("ID", "CONTENT_AREA", "YEAR_INTEGER_TMP", "VALID_CASE")) ## CRITICAL that VALID_CASE is last in group
				if (!"SCALE_SCORE_PRIOR" %in% names(slot.data)) {
					slot.data[,SCALE_SCORE_PRIOR:=slot.data[SJ(ID, CONTENT_AREA, YEAR_INTEGER_TMP-1), mult="last"][["SCALE_SCORE"]]]
				}
				if (!"SGP_PRIOR" %in% names(slot.data)) {
					slot.data[,SGP_PRIOR:=slot.data[SJ(ID, CONTENT_AREA, YEAR_INTEGER_TMP-1), mult="last"][["SGP"]]]
				}
				if (!"SGP_TARGET_PRIOR" %in% names(slot.data)) {
					slot.data[,SGP_TARGET_PRIOR:=slot.data[SJ(ID, CONTENT_AREA, YEAR_INTEGER_TMP-1), mult="last"][["SGP_TARGET_3_YEAR"]]]
				}
				if (!"CONTENT_AREA_PRIOR" %in% names(slot.data) & "SGP_NORM_GROUP" %in% names(slot.data)) {
					slot.data[,CONTENT_AREA_PRIOR:=SGP_NORM_GROUP]
					levels(slot.data$CONTENT_AREA_PRIOR) <- sapply(strsplit(sapply(strsplit(sapply(strsplit(levels(slot.data$CONTENT_AREA_PRIOR), ";"), function(x) rev(x)[2]), "/"), function(x) rev(x)[1]), "_"), function(x) paste(head(x, -1), collapse=" "))
				}
				slot.data[,YEAR_INTEGER_TMP:=NULL]
			}
		}
	} # END Individual Plot setup
 

###################################################################
### bubblePlot style 100 (Individual Student within Grade Chart)
###################################################################

	if (100 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 100", date()))

		### Key slot.data for fast subsetting

		setkeyv(slot.data, c("VALID_CASE", "YEAR", "CONTENT_AREA", "DISTRICT_NUMBER"))

		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas
		for (district.iter in seq_along(my.iters$tmp.districts)) { ### Loop over districts (seq_along to get integer for anonymize)
		
		# Subset data

		tmp.bPlot.data.1 <- slot.data[SJ("VALID_CASE", year.iter, content_area.iter, my.iters$tmp.districts[district.iter])]

		tmp.unique.schools <- my.iters$tmp.schools[my.iters$tmp.schools %in% unique(tmp.bPlot.data.1$SCHOOL_NUMBER)]
		for (school.iter in seq_along(tmp.unique.schools)) { ### Loop over schools (seq_along to get integer for anonymize)

		# Subset data

		tmp.bPlot.data <- tmp.bPlot.data.1[SCHOOL_NUMBER==tmp.unique.schools[school.iter] & !is.na(SGP) & !is.na(SCALE_SCORE)]

		for (grade.iter in intersect(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area.iter]],
			sort(unique(tmp.bPlot.data$GRADE)))) { 
				
		bPlot.data <- subset(tmp.bPlot.data, GRADE==grade.iter)

		if (dim(bPlot.data)[1] > 0) {

		# Anonymize district, school and student names (if requested)

		if (bPlot.anonymize) {
			bPlot.data$FIRST_NAME <- "Student"; bPlot.data$LAST_NAME <- seq(dim(bPlot.data)[1])
			bPlot.data$SCHOOL_NAME <- paste("Sample School", school.iter); bPlot.data$DISTRICT_NAME <- paste("Sample District", district.iter)
		}

		# Create labels

		bPlot.labels <- create.bPlot.labels(year.iter, "SCALE_SCORE_PRIOR", bubblePlot_LEVEL)

		# Create cutscore ranges

		my.content_area <- get.my.cutscore.year(state, content_area.iter, as.character(year.iter))
		tmp.y.range <- extendrange(c(bPlot.data[["SCALE_SCORE"]],
			SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[my.content_area]][[paste("GRADE", grade.iter, sep="_")]]), f=0.1)
		tmp.loss.hoss <- SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[content_area.iter]][[paste("loss.hoss", grade.iter, sep="_")]]
		tmp.y.ticks <- sort(c(max(tmp.loss.hoss[1], tmp.y.range[1]), min(tmp.loss.hoss[2], tmp.y.range[2]),
			SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[my.content_area]][[paste("GRADE", grade.iter, sep="_")]])) 

		# Get median SGP for grade, school, content area combination

		if (bPlot.full.academic.year) {
			school.content_area.grade.median <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR__GRADE__SCHOOL_ENROLLMENT_STATUS"]][
				SCHOOL_ENROLLMENT_STATUS=="Enrolled School: Yes" & SCHOOL_NUMBER==tmp.unique.schools[school.iter] & 
				CONTENT_AREA==content_area.iter & YEAR==year.iter & GRADE==grade.iter][["MEDIAN_SGP"]]
			if (length(school.content_area.grade.median)==0) school.content_area.grade.median <- NA
		} else {
			school.content_area.grade.median <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__CONTENT_AREA__YEAR__GRADE"]][
				SCHOOL_NUMBER==tmp.unique.schools[school.iter] & CONTENT_AREA==content_area.iter & YEAR==year.iter & GRADE==grade.iter][["MEDIAN_SGP"]]
			if (length(school.content_area.grade.median)==0) school.content_area.grade.median <- NA
		}


		### Custom draft message with two median SGP lines

		if (!is.na(school.content_area.grade.median)) {
			school.content_area.grade.median.line <- 
				c(paste("grid.lines(x=unit(", school.content_area.grade.median, ", 'native'), y=c(0.03,0.97), gp=gpar(col='blue', lwd=1.75, lty=2, alpha=0.75))", sep=""),
				paste("grid.text('Grade ", grade.iter, " Median = ", school.content_area.grade.median, "', x=unit(", school.content_area.grade.median,
					", 'native'), y=0.005, gp=gpar(col='blue', cex=0.85))", sep=""))
		} else {
			school.content_area.grade.median.line <- NULL
		}
		bPlot.message.style.100 <- c("grid.text(x=unit(50, 'native'), y=unit(mean(bubble_plot_data.Y), 'native'), 'CONFIDENTIAL \n STUDENT DATA -\n DO NOT DISTRIBUTE', 
				rot=-30, gp=gpar(col='grey80', cex=2, alpha=0.8, fontface=2))",
			school.content_area.grade.median.line)


		### Create bubblePlot ###

		bubblePlot(
			bubble_plot_data.X=bPlot.data[["SGP"]],
			bubble_plot_data.Y=bPlot.data[["SCALE_SCORE"]],
			bubble_plot_data.SUBSET=NULL,
			bubble_plot_data.INDICATE=NULL,
			bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
			bubble_plot_data.SIZE=rep(50, length(bPlot.data[["SGP"]])),
			bubble_plot_data.LEVELS=NULL, 
			bubble_plot_data.BUBBLE_TIPS_LINES=list(
				paste(bPlot.data$SGP, " (", bPlot.data$SGP_TARGET, ")", sep=""),
				paste(bPlot.data$ACHIEVEMENT_LEVEL_PRIOR, " (", bPlot.data$SCALE_SCORE_PRIOR, ")", sep=""),
				paste(bPlot.data$ACHIEVEMENT_LEVEL, " (", bPlot.data$SCALE_SCORE, ")", sep="")),
			bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Student Growth Percentile")),
			bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label$CURRENT),
			bubble_plot_labels.SIZE=NULL,
			bubble_plot_labels.LEVELS=NULL, #levels(bubblePlot[["subset.factor"]]),
			bubble_plot_labels.BUBBLE_TIPS_LINES=list(
				paste(bPlot.labels$x.year.label, "Student Growth Percentile (Target)"),
				paste(bPlot.labels$y.year.label$PRIOR, " (Scale Score)", sep=""),
				paste(bPlot.labels$y.year.label$CURRENT, " (Scale Score)", sep="")),
			bubble_plot_labels.BUBBLE_TITLES=paste(bPlot.data$FIRST_NAME, bPlot.data$LAST_NAME),
			bubble_plot_titles.MAIN=bPlot.labels$main.title,
			bubble_plot_titles.SUB1=paste(bPlot.data$SCHOOL_NAME[1], "Student Performance"),
			bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, "Grade", grade.iter, capwords(content_area.iter)),
			bubble_plot_titles.LEGEND1="",
			bubble_plot_titles.LEGEND2_P1=NULL,
			bubble_plot_titles.LEGEND2_P2=NULL,

			bubble_plot_configs.BUBBLE_MIN_MAX=c(0.07, 0.07),
			bubble_plot_configs.BUBBLE_X_TICKS=c(1, seq(10,90,10), 99),
			bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.7, 5), 1, rep(0.7, 5)),
			bubble_plot_configs.BUBBLE_Y_TICKS=tmp.y.ticks,
			bubble_plot_configs.BUBBLE_Y_BANDS=tmp.y.ticks,
			bubble_plot_configs.BUBBLE_Y_BAND_LABELS=SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]],
			bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
			bubble_plot_configs.BUBBLE_COLOR="blue",
			bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
			bubble_plot_configs.BUBBLE_TIPS="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
			bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
			bubble_plot_configs.BUBBLE_PLOT_LEGEND="FALSE",
			bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
			bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS=NULL,
			bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message.style.100,
			bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(gsub(" ", "_", bPlot.data$SCHOOL_NAME[1]), "Grade", grade.iter,
				year.iter, capwords(content_area.iter), bPlot.labels$pdf.title, sep="_"), ".pdf", sep=""),
			bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "Individual", "Style_100", gsub(" ", "_", bPlot.data$DISTRICT_NAME[1])),
			bubble_plot_pdftk.CREATE_CATALOG=FALSE)

		} ## END if dim(bPlot.data)[1] > 0
		} ## END grade.iter loop
		} ## END school.iter loop
		} ## END district.iter loop
		} ## END content_area.iter loop
		} ## END year.iter loop

		message(paste("\tFinished bubblePlot Style 100", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END if bubblePlot style 100

######################################################################
### bubblePlot style 150 (Individual Student within INSTRUCTOR Chart)
######################################################################

	if (150 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 150", date()), "\n")

		if (bPlot.demo) {
			bPlot.anonymize <- TRUE
			setkeyv(slot.data, c("VALID_CASE", "YEAR", "CONTENT_AREA", "GRADE"))
			my.iters$tmp.districts <- '-999'; my.iters$tmp.schools <- '-99'
		} else setkeyv(slot.data, c("VALID_CASE", "YEAR", "CONTENT_AREA", "DISTRICT_NUMBER"))


		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas
		for (district.iter in seq_along(my.iters$tmp.districts)) { ### Loop over districts (seq_along to get integer for anonymize)

			# Subset data

			if (bPlot.demo) {
				tmp.ids <- list()
				tmp.grades.reported <- SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area.iter]][-1]
				tmp.grades.reported <- as.character(tmp.grades.reported[tmp.grades.reported %in% unique(slot.data)["VALID_CASE"][["GRADE"]]])
				for (i in seq_along(tmp.grades.reported)) {
					tmp.ids[[i]] <- as.character(sample(unique(slot.data[SJ("VALID_CASE", year.iter, content_area.iter, tmp.grades.reported[i])][['ID']]), 30))
				}
				
				tmp.bPlot.data.1.long <- slot.data[SJ("VALID_CASE", year.iter, content_area.iter)][ID %in% unlist(tmp.ids)]
	
				tmp.bPlot.data.1.long[['INSTRUCTOR_NUMBER']] <- factor(paste("Grade_", tmp.bPlot.data.1.long[['GRADE']], sep=""))
				tmp.bPlot.data.1.long[['INSTRUCTOR_NAME']] <- factor('Psuedo-Instructor')
				tmp.bPlot.data.1.long[['SCHOOL_NUMBER']] <- factor('-99')
				tmp.bPlot.data.1.long[['SCHOOL_NAME']] <- factor('Psuedo School')
				tmp.bPlot.data.1.long[['DISTRICT_NUMBER']] <- factor('-999')
				tmp.bPlot.data.1.long[['DISTRICT_NAME']] <- factor('Psuedo District')
				setkeyv(tmp.bPlot.data.1.long, "INSTRUCTOR_NUMBER")
			} else {
				if (!"INSTRUCTOR_NUMBER" %in% names(sgp_object@Data_Supplementary)) {
					stop("\tNOTE: Indvidividual level Instructor bubble plots require an INSTRUCTOR_NUMBER lookup table embedded in @Data_Supplementary")
				}
				tmp.bPlot.data.1.long <- data.table(slot.data[SJ("VALID_CASE", year.iter, content_area.iter, my.iters$tmp.districts[district.iter])],
                                                                key=c("ID", "CONTENT_AREA", "YEAR"))[sgp_object@Data_Supplementary$INSTRUCTOR_NUMBER, nomatch=0]
			}

			setkeyv(tmp.bPlot.data.1.long, "INSTRUCTOR_NUMBER")
			tmp.unique.schools <- my.iters$tmp.schools[my.iters$tmp.schools %in% unique(tmp.bPlot.data.1.long$SCHOOL_NUMBER)]
			for (school.iter in seq_along(tmp.unique.schools)) { ### Loop over schools (seq_along to get integer for anonymize)
	
			# Subset data
	
			tmp.bPlot.data <- tmp.bPlot.data.1.long[SCHOOL_NUMBER==tmp.unique.schools[school.iter] & !is.na(SGP)]
			
			for (instructor.iter in sort(unique(tmp.bPlot.data$INSTRUCTOR_NUMBER))) { ### Loop over unique teachers in school
				bPlot.data <- tmp.bPlot.data[SJ(instructor.iter)]
		
				if (dim(bPlot.data)[1] > 0) {
		
				for (grade.iter in intersect(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area.iter]],
					sort(unique(bPlot.data$GRADE)))) { ### Loop over unique grades levels for instructor (usually only one)
			
			# Anonymize district, school and student names (if requested)
	
				if (bPlot.anonymize) {
					bPlot.data$FIRST_NAME <- "Student"; bPlot.data$LAST_NAME <- seq(dim(bPlot.data)[1])
					bPlot.data$SCHOOL_NAME <- paste("Psuedo School", school.iter); bPlot.data$DISTRICT_NAME <- paste("Psuedo District", district.iter)
				}
	
			# Create labels and file path
	
				bPlot.labels <- create.bPlot.labels(year.iter, "SCALE_SCORE", bubblePlot_LEVEL) # Only produce "Current Year" plots
				prior.year <- paste(as.numeric(unlist(strsplit(as.character(year.iter), "_")))-1, collapse="_")
	
			# Create cutscore ranges
	
				my.content_area <- get.my.cutscore.year(state, content_area.iter, as.character(bPlot.labels$y.year)) 
				tmp.y.range <- extendrange(c(bPlot.data[["SCALE_SCORE"]], 
					SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[my.content_area]][[paste("GRADE", grade.iter, sep="_")]]), f=0.1)
				tmp.loss.hoss <- SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[content_area.iter]][[paste("loss.hoss", grade.iter, sep="_")]]
				tmp.y.ticks <- sort(c(max(tmp.loss.hoss[1], tmp.y.range[1]), min(tmp.loss.hoss[2], tmp.y.range[2]),
					SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[my.content_area]][[paste("GRADE", grade.iter, sep="_")]])) 
	
			# Get median SGP for grade, school, content area combination
	
				instructor.content_area.grade.median <- median(bPlot.data$SGP, na.rm=TRUE)
	
			### Custom draft message with two median SGP lines
	
				bPlot.message.style.150 <- c("grid.text(x=unit(50, 'native'), y=unit(mean(bubble_plot_data.Y), 'native'), 'CONFIDENTIAL \n STUDENT DATA -\n DO NOT DISTRIBUTE', 
						rot=-30, gp=gpar(col='grey80', cex=2, alpha=0.8, fontface=2))", 
					paste("grid.lines(x=unit(", instructor.content_area.grade.median, 
						", 'native'), y=c(0.03,0.97), gp=gpar(col='blue', lwd=1.75, lty=2, alpha=0.75))", sep=""),
					paste("grid.text('Classroom Median = ", instructor.content_area.grade.median, "', x=unit(", 
						instructor.content_area.grade.median, ", 'native'), y=0.005, gp=gpar(col='blue', cex=0.85))", sep=""))
	
			### Create bubblePlot ###
	
				bubblePlot(
					bubble_plot_data.X=bPlot.data[['SGP']],
					bubble_plot_data.Y=bPlot.data[['SCALE_SCORE']],
					bubble_plot_data.SUBSET=NULL,
					bubble_plot_data.INDICATE=NULL,
					bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
					bubble_plot_data.SIZE=rep(50, length(bPlot.data[['SGP']])),
					bubble_plot_data.LEVELS=NULL, 
	
					bubble_plot_data.BUBBLE_TIPS_LINES=list(paste(bPlot.data[['SGP']], " (", bPlot.data[['SGP_TARGET']], ")", sep=""),
						bPlot.data[['SCALE_SCORE']],
						paste(bPlot.data[['SGP_PRIOR']], " (", bPlot.data[['SGP_TARGET_PRIOR']], ")", sep=""),
						paste(bPlot.data[['SCALE_SCORE_PRIOR']], " (", bPlot.data[['ACHIEVEMENT_LEVEL_PRIOR']], ")", sep="")),
					bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Student Growth Percentile (Target)"),
						paste(bPlot.labels$x.year.label, "Scale Score"),
						paste(prior.year, "Prior Student Growth Percentile (Target)"),
						paste(prior.year, "Prior Scale Score (Achievement)")),
	
					bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Student Growth Percentile")),
					bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
					bubble_plot_labels.SIZE=NULL,
					bubble_plot_labels.LEVELS=NULL,
					bubble_plot_labels.BUBBLE_TITLES=paste(bPlot.data$FIRST_NAME, bPlot.data$LAST_NAME),
					bubble_plot_titles.MAIN=bPlot.labels$main.title,
					bubble_plot_titles.SUB1=paste(bPlot.data$SCHOOL_NAME[1], "Student Performance"),
					bubble_plot_titles.SUB2=paste(bPlot.labels$x.year.label, "Instructor", instructor.iter, capwords(content_area.iter)),
					bubble_plot_titles.LEGEND1="",
					bubble_plot_titles.LEGEND2_P1=NULL,
					bubble_plot_titles.LEGEND2_P2=NULL,
		
					bubble_plot_configs.BUBBLE_MIN_MAX=c(0.07, 0.07),
					bubble_plot_configs.BUBBLE_X_TICKS=c(1, seq(10,90,10), 99),
					bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.7, 5), 1, rep(0.7, 5)),
					bubble_plot_configs.BUBBLE_Y_TICKS=tmp.y.ticks,
					bubble_plot_configs.BUBBLE_Y_BANDS=tmp.y.ticks,
					bubble_plot_configs.BUBBLE_Y_BAND_LABELS=SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]],
					bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
					bubble_plot_configs.BUBBLE_COLOR="blue",
					bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
					bubble_plot_configs.BUBBLE_TIPS="TRUE",
					bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
					bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
					bubble_plot_configs.BUBBLE_PLOT_LEGEND="FALSE",
					bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
					bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS=NULL,
					bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message.style.150,
					bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste("Instructor", instructor.iter, year.iter, "Grade", grade.iter, 
						capwords(content_area.iter), "Student_Plot", sep="_"), ".pdf", sep=""),
					bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "Instructor", "Style_150", 
						gsub(" ", "_", bPlot.data$DISTRICT_NAME[1]), gsub(" ", "_", bPlot.data$SCHOOL_NAME[1])),
					bubble_plot_pdftk.CREATE_CATALOG=FALSE)
	
				} ## END grade.iter loop
				} ## END if dim(bPlot.data)[1] > 0
				} ## END instructor.iter loop
				} ## END school.iter loop
			} ## END district.iter loop
			} ## END content_area.iter loop
			} ## END year.iter loop

		message(paste("\tFinished bubblePlot Style 150", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END if bubblePlot style 150


######################################################################
### bubblePlot style 153 (Individual Student within INSTRUCTOR Chart) - course sequence specific
######################################################################


	if (153 %in% bPlot.styles) {

		started.at <- proc.time()
		message(paste("\tStarted bubblePlot Style 153", date()), "\n")

		if (bPlot.demo) {
			bPlot.anonymize <- TRUE
			setkeyv(slot.data, c("VALID_CASE", "YEAR", "CONTENT_AREA", "GRADE"))
			my.iters$tmp.districts <- '-999'; my.iters$tmp.schools <- '-99'
		} else setkeyv(slot.data, c("VALID_CASE", "YEAR", "CONTENT_AREA", "DISTRICT_NUMBER"))


		### Start loops for bubblePlots

		for (year.iter in my.iters$tmp.years) {  ### Loop over year
		for (content_area.iter in my.iters$tmp.content_areas) { ### Loop over content areas
		for (district.iter in seq_along(my.iters$tmp.districts)) { ### Loop over districts (seq_along to get integer for anonymize)

			# Subset data

			if (bPlot.demo) {
				tmp.ids <- list()
				tmp.grades.reported <- SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area.iter]][-1]
				tmp.grades.reported <- as.character(tmp.grades.reported[tmp.grades.reported %in% unique(slot.data)["VALID_CASE"][["GRADE"]]])
				for (i in seq_along(tmp.grades.reported)) {
					tmp.ids[[i]] <- as.character(sample(unique(slot.data[SJ("VALID_CASE", year.iter, content_area.iter, tmp.grades.reported[i])][['ID']]), 30))
				}
				
				tmp.bPlot.data.1.long <- slot.data[SJ("VALID_CASE", year.iter, content_area.iter)][ID %in% unlist(tmp.ids)]
	
				tmp.bPlot.data.1.long[['INSTRUCTOR_NUMBER']] <- factor(paste("Grade_", tmp.bPlot.data.1.long[['GRADE']], sep=""))
				tmp.bPlot.data.1.long[['INSTRUCTOR_NAME']] <- factor('Psuedo-Instructor')
				tmp.bPlot.data.1.long[['SCHOOL_NUMBER']] <- factor('-99')
				tmp.bPlot.data.1.long[['SCHOOL_NAME']] <- factor('Pseudo School')
				tmp.bPlot.data.1.long[['DISTRICT_NUMBER']] <- factor('-999')
				tmp.bPlot.data.1.long[['DISTRICT_NAME']] <- factor('Psuedo District')
			} else {
				if (!"INSTRUCTOR_NUMBER" %in% names(sgp_object@Data_Supplementary)) {
					stop("\tNOTE: Indvidividual level Instructor bubble plots require an INSTRUCTOR_NUMBER lookup table embedded in @Data_Supplementary")
				}
				tmp.bPlot.data.1.long <- data.table(slot.data[SJ("VALID_CASE", year.iter, content_area.iter, my.iters$tmp.districts[district.iter])],
								key=c("ID", "CONTENT_AREA", "YEAR"))[sgp_object@Data_Supplementary$INSTRUCTOR_NUMBER, nomatch=0]
			}
		
			tmp.unique.schools <- my.iters$tmp.schools[my.iters$tmp.schools %in% unique(tmp.bPlot.data.1.long$SCHOOL_NUMBER)]
			for (school.iter in seq_along(tmp.unique.schools)) { ### Loop over schools (seq_along to get integer for anonymize)
	
			# Subset data
	
			tmp.bPlot.data <- tmp.bPlot.data.1.long[SCHOOL_NUMBER==tmp.unique.schools[school.iter] & !is.na(SGP)]
	
			for (instructor.iter in sort(unique(tmp.bPlot.data$INSTRUCTOR_NUMBER))) { ### Loop over unique teachers in school
			bPlot.data <- subset(tmp.bPlot.data, INSTRUCTOR_NUMBER==instructor.iter)
			if ("INSTRUCTOR_LAST_NAME" %in% names(bPlot.data)) {
				tmp.instructor.reference <- paste(bPlot.data$INSTRUCTOR_FIRST_NAME[1], bPlot.data$INSTRUCTOR_LAST_NAME[1])
			} else {
				tmp.instructor.reference <- instructor.iter
			}
	
	
			if (dim(bPlot.data)[1] > 1) { # had error when only 1 kid per teacher
	
			for (grade.iter in sort(unique(bPlot.data$GRADE))) { ### Loop over grades levels. Only one in course sequences
		
		# Anonymize district, school and student names (if requested)

			if (bPlot.anonymize) {
				bPlot.data$FIRST_NAME <- "Student"; bPlot.data$LAST_NAME <- seq(dim(bPlot.data)[1])
				bPlot.data$SCHOOL_NAME <- paste("Psuedo School", school.iter); bPlot.data$DISTRICT_NAME <- paste("Psuedo District", district.iter)
			}

		# Create labels

			bPlot.labels <- create.bPlot.labels(year.iter, "SCALE_SCORE", bubblePlot_LEVEL) # Only produce "Current Year" plots
			prior.year <- paste(as.numeric(unlist(strsplit(as.character(year.iter), "_")))-1, collapse="_")

		# Create cutscore ranges

			my.content_area <- get.my.cutscore.year(state, content_area.iter, as.character(bPlot.labels$y.year)) 
			tmp.y.range <- extendrange(c(bPlot.data[["SCALE_SCORE"]], 
				SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[my.content_area]][[paste("GRADE", grade.iter, sep="_")]]), f=0.1)
			tmp.loss.hoss <- SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[content_area.iter]][[paste("loss.hoss", grade.iter, sep="_")]]
			tmp.y.ticks <- sort(c(max(tmp.loss.hoss[1], tmp.y.range[1]), min(tmp.loss.hoss[2], tmp.y.range[2]),
				SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[my.content_area]][[paste("GRADE", grade.iter, sep="_")]])) 

		# Get median SGP for grade, school, content area combination
		# Report 'Official' Median.  Should be the same as median(bPlot.data$SGP, na.rm=TRUE).  Use that if NULL for some reason (prevents error)

			if (bPlot.full.academic.year) {
				instructor.content_area.grade.median <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR__GRADE__INSTRUCTOR_ENROLLMENT_STATUS"]][
					SCHOOL_NUMBER==tmp.unique.schools[school.iter] & INSTRUCTOR_NUMBER==instructor.iter & CONTENT_AREA==content_area.iter & YEAR==year.iter & GRADE==grade.iter][["MEDIAN_SGP"]]
			} else {
				instructor.content_area.grade.median <- sgp_object@Summary[["SCHOOL_NUMBER"]][["SCHOOL_NUMBER__INSTRUCTOR_NUMBER__CONTENT_AREA__YEAR__GRADE"]][
					SCHOOL_NUMBER==tmp.unique.schools[school.iter] & INSTRUCTOR_NUMBER==instructor.iter & CONTENT_AREA==content_area.iter & YEAR==year.iter & GRADE==grade.iter][["MEDIAN_SGP"]]
			}			
			if (is.null(instructor.content_area.grade.median)) instructor.content_area.grade.median <- median(bPlot.data$SGP, na.rm=TRUE)
			if (bPlot.demo) instructor.content_area.grade.median <- median(bPlot.data$SGP, na.rm=TRUE)

		### Custom message with two median SGP lines

		bPlot.message.style.153 <- c("grid.text(x=unit(50, 'native'), y=unit(mean(bubble_plot_data.Y), 'native'), 'CONFIDENTIAL \n STUDENT DATA -\n DO NOT DISTRIBUTE', 
				rot=-30, gp=gpar(col='grey80', cex=2, alpha=0.8, fontface=2))", 
			paste("grid.lines(x=unit(", instructor.content_area.grade.median, ", 'native'), y=c(0.03,0.97), gp=gpar(col='blue', lwd=1.75, lty=2, alpha=0.75))", sep=""),
			paste("grid.text('Classroom Median = ", instructor.content_area.grade.median, "', x=unit(", instructor.content_area.grade.median, 
				", 'native'), y=0.005, gp=gpar(col='blue', cex=0.85))", sep=""))

		### Create bubblePlot ###

			bubblePlot(
				bubble_plot_data.X=bPlot.data[['SGP']],
				bubble_plot_data.Y=bPlot.data[['SCALE_SCORE']],
				bubble_plot_data.SUBSET=NULL,
				bubble_plot_data.INDICATE=NULL,
				bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
				bubble_plot_data.SIZE=rep(50, length(bPlot.data[['SGP']])),
				bubble_plot_data.LEVELS=NULL, 

				bubble_plot_data.BUBBLE_TIPS_LINES=list(bPlot.data[['SGP']], bPlot.data[['SCALE_SCORE']],
					sapply(gsub("_", " ", as.character(bPlot.data[["CONTENT_AREA_PRIOR"]])), capwords),  # capwords not working with factors with "_" in them...
					bPlot.data[['SGP_PRIOR']],
					paste(bPlot.data[['SCALE_SCORE_PRIOR']], " (", bPlot.data[['ACHIEVEMENT_LEVEL_PRIOR']], ")", sep="")),
				bubble_plot_labels.BUBBLE_TIPS_LINES=list(paste(bPlot.labels$x.year.label, "Student Growth Percentile"),
					paste(bPlot.labels$x.year.label, "Scale Score"),
					paste(prior.year, "Prior Subject / Norm Group"),
					paste(prior.year, "Prior Student Growth Percentile"),
					paste(prior.year, "Prior Scale Score (Achievement)")),

				bubble_plot_labels.X=c("Growth", paste(bPlot.labels$x.year.label, "Student Growth Percentile")),
				bubble_plot_labels.Y=c("Achievement", bPlot.labels$y.year.label),
				bubble_plot_labels.SIZE=NULL,
				bubble_plot_labels.LEVELS=NULL, #levels(bubblePlot[["subset.factor"]]),
				bubble_plot_labels.BUBBLE_TITLES=paste(bPlot.data$FIRST_NAME, bPlot.data$LAST_NAME),
				bubble_plot_titles.MAIN=bPlot.labels$main.title,
				bubble_plot_titles.SUB1=paste(bPlot.labels$x.year.label, bPlot.data$SCHOOL_NAME[1], "Student", capwords(content_area.iter), "Performance"),
				bubble_plot_titles.SUB2=paste("Instructor:", tmp.instructor.reference),
				bubble_plot_titles.LEGEND1="",
				bubble_plot_titles.LEGEND2_P1=NULL,
				bubble_plot_titles.LEGEND2_P2=NULL,
	
				bubble_plot_configs.BUBBLE_MIN_MAX=c(0.07, 0.07),
				bubble_plot_configs.BUBBLE_X_TICKS=c(1, seq(10,90,10), 99),
				bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.7, 5), 1, rep(0.7, 5)),
				bubble_plot_configs.BUBBLE_Y_TICKS=tmp.y.ticks,
				bubble_plot_configs.BUBBLE_Y_BANDS=tmp.y.ticks,
				bubble_plot_configs.BUBBLE_Y_BAND_LABELS=SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]],
				bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0.00,
				bubble_plot_configs.BUBBLE_COLOR="blue",
				bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.9),
				bubble_plot_configs.BUBBLE_TIPS="TRUE",
				bubble_plot_configs.BUBBLE_PLOT_DEVICE=bPlot.output,
				bubble_plot_configs.BUBBLE_PLOT_FORMAT=bPlot.format,
				bubble_plot_configs.BUBBLE_PLOT_LEGEND="FALSE",
				bubble_plot_configs.BUBBLE_PLOT_TITLE="TRUE",
				bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS=NULL,
				bubble_plot_configs.BUBBLE_PLOT_EXTRAS=bPlot.message.style.153,
				bubble_plot_configs.BUBBLE_PLOT_NAME=paste(paste(gsub(" ", "_", bPlot.data$SCHOOL_NAME[1]), "Instructor", 
					instructor.iter, year.iter, capwords(content_area.iter), "Student_Plot", sep="_"), ".pdf", sep=""),
				bubble_plot_configs.BUBBLE_PLOT_PATH=file.path(bPlot.folder, year.iter, "Instructor", "Style_153", 
					gsub(" ", "_", bPlot.data$DISTRICT_NAME[1]), gsub(" ", "_", bPlot.data$SCHOOL_NAME[1])),
				bubble_plot_pdftk.CREATE_CATALOG=FALSE)

			} ## END grade.iter loop
			} ## END if dim(bPlot.data)[1] > 0
			} ## END instructor.iter loop
			} ## END school.iter loop
		} ## END district.iter loop
		} ## END content_area.iter loop
		} ## END year.iter loop

		message(paste("\tFinished bubblePlot Style 153", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END if bubblePlot style 153

####
#### END bubblePlot_Styles
####

} ## END bubblePlot_Styles function
