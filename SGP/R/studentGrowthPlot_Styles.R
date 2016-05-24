`studentGrowthPlot_Styles` <- 
	function(
		sgPlot.data,
		sgPlot.sgp_object,
		sgPlot.cutscores,
		state,
		last.year,
		content_areas,
		districts,
		schools,
		reports.by.student,
		reports.by.instructor,
		reports.by.school,
		sgPlot.years,
		sgPlot.demo.report,
		sgPlot.folder,
		sgPlot.folder.names,
		sgPlot.anonymize,
		sgPlot.front.page,
		sgPlot.header.footer.color,
		sgPlot.fan,
		sgPlot.sgp.targets,
		sgPlot.cleanup,
		sgPlot.baseline,
		sgPlot.sgp.targets.timeframe,
		sgPlot.zip,
		sgPlot.output.format,
		sgPlot.linkages) {

	CUTLEVEL <- ID <- CONTENT_AREA <- GRADE <- CUTSCORES <- YEAR <- NULL ## To prevent R CMD check warnings

	### Utility functions

	"%w/o%" <- function(x,y) x[!x %in% y]

	pretty_year <- function(x) sub("_", "-", x)


	### Check/adjust supplied arguments

	if (any(c("PNG", "PDF_PIECES") %in% sgPlot.output.format | reports.by.student)) {
		reports.by.school <- TRUE
		reports.by.instructor <- FALSE
	}


	### Define quantities/variables related to state

	if (state %in% objects(SGP::SGPstateData)) {
		tmp.abbreviation <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Abbreviation"]]
		tmp.state <- paste(datasets::state.name[state==datasets::state.abb], tmp.abbreviation)
		tmp.organization <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Organization"]]
		number.achievement.level.regions <- length(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Achievement_Level_Labels"]])
		if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["sgp.projections.max.forward.progression.grade"]])) {
			trajectory.cuts <- sort(c(SGP::SGPstateData[[state]][["Growth"]][["Cutscores"]][['Cuts']], SGP::SGPstateData[[state]][["Student_Report_Information"]][["Projection_Fan_Limits"]]))
			trajectory.cuts <- paste(paste("P", trajectory.cuts, sep=""), collapse="|")
		} else {
			trajectory.cuts <- c(1, SGP::SGPstateData[[state]][["Growth"]][["Cutscores"]][['Cuts']], 99)
			trajectory.cuts <- paste(paste("P", trajectory.cuts, sep=""), collapse="|")
		}
		if (!is.null(SGP::SGPstateData[[state]][["Custom_Student_Report"]])) {
			custom.isr.tf <- TRUE
			custom.isr <- SGP::SGPstateData[[state]][["Custom_Student_Report"]]
		} else custom.isr.tf <- FALSE
	} else {
		stop("Construction of student growth plots requires state meta-data to be included in the embedded SGPstateData set.\nPlease augment the SGPstateData set with your data or contact the SGP package maintainer to have your data added to the SGP package.")
	}	

	year_folder <- as.character(last.year) 

	## To Baseline or not to Baseline

	if (sgPlot.baseline) {
		my.sgp <- "SGP_BASELINE"
		my.sgp.level <- "SGP_LEVEL_BASELINE"
		my.sgp.target.label <- c("SGP_TARGET_BASELINE", paste(sgPlot.sgp.targets.timeframe, "YEAR", sep="_"))
	} else {
		my.sgp <- "SGP"
		my.sgp.level <- "SGP_LEVEL"
		my.sgp.target.label <- c("SGP_TARGET", paste(sgPlot.sgp.targets.timeframe, "YEAR", sep="_"))
	}


#############################################################################################################################
###
### Reports assuming DISTRICT, SCHOOLS, GRADES, and then STUDENTS
###
#############################################################################################################################

if (reports.by.school) {

	### Extend sgPlot.folder

	sgPlot.folder <- file.path(sgPlot.folder, "School")

	if (identical(sgPlot.output.format, "PNG")) {
		sgPlot.folder <- file.path(sgPlot.folder, "PNG")
	}


	### Loop over unique DISTRICTS, SCHOOLS, GRADES and then STUDENTS

	tmp.keys <- paste(c("DISTRICT_NUMBER", "SCHOOL_NUMBER", "GRADE", "LAST_NAME", "FIRST_NAME"), last.year, sep=".")

	## Districts

	setkeyv(sgPlot.data, tmp.keys[1])

	for (i in districts) {
		if (reports.by.student) {
			tmp_district_name <- as.character(sgPlot.data[list(i)][[paste("DISTRICT_NAME", last.year, sep=".")]][1])
			district_folder <- "Uncollated_Student_Reports"
		} else {
			if (sgPlot.demo.report | identical(i, -999L)) {
				tmp_district_name <- "Sample District"
				district_folder <- "Sample_District"
			} else {
				if (sgPlot.folder.names=="name") {
					tmp_district_name <- as.character(sgPlot.data[list(i)][[paste("DISTRICT_NAME", last.year, sep=".")]][1])
						district_folder <- gsub(" ", "_", paste(tmp_district_name, i))
				} else {
					tmp_district_name <- as.character(sgPlot.data[list(i)][[paste("DISTRICT_NAME", last.year, sep=".")]][1])
						district_folder <- as.character(i)
				}
			}
		}
		tmp_district_ids <- unique(sgPlot.data[list(i)]$ID)
		tmp_district_data <- subset(sgPlot.data, ID %in% tmp_district_ids)
		if (is.factor(tmp_district_data[[paste("LAST_NAME", last.year, sep=".")]]) | is.factor(tmp_district_data[[paste("FIRST_NAME", last.year, sep=".")]])) {
			tmp_district_data[, paste("LAST_NAME", last.year, sep=".") := as.character(tmp_district_data[[paste("LAST_NAME", last.year, sep=".")]])]
			tmp_district_data[, paste("FIRST_NAME", last.year, sep=".") := as.character(tmp_district_data[[paste("FIRST_NAME", last.year, sep=".")]])]			
		}

	## Schools

	setkeyv(tmp_district_data, tmp.keys[2])

	for (j in intersect(tmp_district_data[[paste("SCHOOL_NUMBER", last.year, sep=".")]], schools)) {
	
		started.at <- proc.time()
		started.date <- date()
	
		if (!is.null(sgPlot.front.page)) {
			sgPlot.front.page.ij <- paste(paste(unlist(strsplit(sgPlot.front.page, "[.]"))[1], i, j, sep="_"), ".pdf", sep="")
			file.copy(sgPlot.front.page, sgPlot.front.page.ij)
		}
	
		if (reports.by.student) {
			tmp_school_name <- as.character(tmp_district_data[list(j)][[paste("SCHOOL_NAME", last.year, sep=".")]][1])
			school_folder <- NULL
		} else {
			if (sgPlot.demo.report | identical(j, -99L)) {
				tmp_school_name <- "Sample School"
				school_folder <- "Sample_School"
			} else {
				if (sgPlot.folder.names=="name") {
					tmp_school_name <- as.character(tmp_district_data[list(j)][[paste("SCHOOL_NAME", last.year, sep=".")]][1])
					school_folder <- gsub(" ", "_", paste(tmp_school_name, j))
				} else {
					tmp_school_name <- as.character(tmp_district_data[list(j)][[paste("SCHOOL_NAME", last.year, sep=".")]][1])
					school_folder <- as.character(j)
				}
			}
	
			if ("PDF" %in% sgPlot.output.format) {
				######################## SCHOOL Report Catalog LaTeX Header #################################################################################
				if (.Platform$OS.type == "windows" & length(unique(tmp_district_data[list(j)]$ID)) > 500) {
					cat("\\let\\mypdfximage\\pdfximage\n\\def\\pdfximage{\\immediate\\mypdfximage}\n\\documentclass[pdftex]{book}\n\\usepackage{hyperref,pdfpages}\n\\hypersetup{%\n",
						file=paste("school_catalog_", i, "_", j, ".tex", sep=""))	
				} else cat("\\documentclass[pdftex]{book}\n\\usepackage{hyperref,pdfpages}\n\\hypersetup{%\n", file=paste("school_catalog_", i, "_", j, ".tex", sep=""))
				cat(paste("pdftitle={", tmp_school_name, ": ", pretty_year(last.year), " ", tmp.state, " Growth and Achievement Reports},\n", sep=""), 
					file=paste("school_catalog_", i, "_", j, ".tex", sep=""), append=TRUE)
				cat(paste("pdfauthor={", tmp.organization$Name, "/Center for Assessment Inc.},\n", sep=""), "pdfcreator={pdfLaTeX},\n", 
				paste("pdfproducer={", tmp.organization$Name, "/Center for Assessment Inc.}}\n", sep=""), 
					"\\begin{document}\n", sep="", file=paste("school_catalog_", i, "_", j, ".tex", sep=""), append=TRUE)
				cat(paste("\\pdfbookmark[-1]{", tmp_district_name, "}{", i, "}\n", sep=""), file=paste("school_catalog_", i, "_", j, ".tex", sep=""), append=TRUE)
				cat(paste("\\pdfbookmark[0]{", tmp_school_name, "}{", j, "}\n", sep=""), file=paste("school_catalog_", i, "_", j, ".tex", sep=""), append=TRUE)
				##############################################################################################################################################
			}
		}
		tmp_school_ids <- unique(tmp_district_data[list(j)]$ID)
		tmp_school_data <- subset(tmp_district_data, ID %in% tmp_school_ids)
	
		## Grades
	
		setkeyv(tmp_school_data, tmp.keys[2])
		grades <- as.character(sort(type.convert(unique(unlist(tmp_school_data[list(j), tmp.keys[3], with=FALSE])), as.is=TRUE)) %w/o% NA)
		setkeyv(tmp_school_data, tmp.keys[3])
	
		for (k in grades) {
		if (reports.by.student) {
			grade_folder <- NULL
		} else {
			if (sgPlot.folder.names=="name") {
				grade_folder <- paste("Grade", k, sep="_")
			} else {
				if (k=="EOCT") {
					grade_folder <- "EOCT"
				} else {
					grade_folder <- substr(paste("0", as.character(k), sep=""), nchar(k), nchar(k)+1)
				}
			}
	
			if ("PDF" %in% sgPlot.output.format) {
				################################ SCHOOL Report Catalog LaTeX Code #########################################################
				cat(paste("\\pdfbookmark[1]{Grade ", k, "}{", j, k, "}\n", sep=""), 
					file=paste("school_catalog_", i, "_", j, ".tex", sep=""), append=TRUE) ## NOTE: j, k included in anchor for uniqueness
				###########################################################################################################################
			}
		}
		tmp_grade_ids <- unique(tmp_school_data[list(k)]$ID)
		tmp_grade_data <- subset(tmp_school_data, ID %in% tmp_grade_ids)
	
		### Create path to pdf files
	
		path.to.pdfs <- paste(c(sgPlot.folder, year_folder, district_folder, school_folder, grade_folder), collapse=.Platform$file.sep)
		dir.create(path.to.pdfs, recursive=TRUE, showWarnings=FALSE)
	
		## Students
	
		setkeyv(tmp_grade_data, tmp.keys[4:5])
	
		for (n in unique(tmp_grade_data[["ID"]])) {
			FIRST_NAME <- sort(tmp_grade_data[ID==n][[tmp.keys[5]]])[1]
			LAST_NAME <- sort(tmp_grade_data[ID==n][[tmp.keys[4]]])[1]
			if (sgPlot.anonymize) {
				student_number <- 1234567890
			} else {
				student_number <- n
			}
			if (sgPlot.folder.names=="name" | sgPlot.anonymize) {
				file_name <- paste(paste(gsub(" |/", "-", FIRST_NAME), gsub(" |/", "-", LAST_NAME), student_number, year_folder, "REPORT", sep="_"), ".pdf", sep="")
				file_name_json <- paste(paste(gsub(" |/", "-", FIRST_NAME), gsub(" |/", "-", LAST_NAME), student_number, year_folder, sep="_"), sep="")
			} else {
				file_name <- paste(n, "_REPORT.pdf", sep="")
				file_name_json <- n
			}
	
		tmp_all_student_data <- tmp_grade_data[ID==n]
		if (length(content_areas) < dim(tmp_all_student_data)[1] | any(duplicated(tmp_all_student_data[["CONTENT_AREA"]]))) {
			tmp_content_areas <- tmp_all_student_data[["CONTENT_AREA"]]
			tmp_index <- unlist(sapply(unique(tmp_content_areas), function(x) which(tmp_content_areas==x), USE.NAMES=FALSE))
			tmp_all_student_data[["CONTENT_AREA"]] <- tmp_content_areas <- unlist(sapply(unique(tmp_content_areas), function(x) paste(x, which(tmp_content_areas==x), sep="."), USE.NAMES=FALSE))[tmp_index]
		} else tmp_content_areas <- content_areas
		num.charts <- length(tmp_content_areas)

		if ("JSON" %in% sgPlot.output.format) {
	
			for (vp in seq_along(tmp_content_areas)) {
				tmp_student_data <- as.data.frame(tmp_all_student_data[CONTENT_AREA==tmp_content_areas[vp]])
				if (dim(tmp_student_data)[1]==0) tmp_student_data <- as.data.frame(data.table(tmp_student_data)[NA][,CONTENT_AREA := tmp_content_areas[vp]])
				tmp.list <- list(
					Scale_Scores_TEXT=as.numeric(subset(tmp_student_data, select=paste("SCALE_SCORE", rev(sgPlot.years), sep="."))),
					Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("TRANSFORMED_SCALE_SCORE", rev(sgPlot.years), sep="."))),
					Achievement_Levels=as.character(unlist(subset(tmp_student_data, select=paste("ACHIEVEMENT_LEVEL", rev(sgPlot.years), sep=".")))),
					SGP=as.numeric(subset(tmp_student_data, select=paste(my.sgp, rev(sgPlot.years), sep="."))),
					SGP_Levels=as.character(unlist(subset(tmp_student_data, select=paste(my.sgp.level, rev(sgPlot.years), sep=".")))),
					Grades=as.character(subset(tmp_student_data, select=paste("GRADE", rev(sgPlot.years), sep="."))),
					Content_Areas=as.character(subset(tmp_student_data, select=paste("CONTENT_AREA_LABELS", rev(sgPlot.years), sep="."))),
					Cuts_TEXT=list(
						NY1=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_1", names(tmp_student_data))), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY2=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_2", names(tmp_student_data))), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY3=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_3", names(tmp_student_data))), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
				Cuts=list(NY1=as.numeric(subset(tmp_student_data, 
						select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY2=as.numeric(subset(tmp_student_data, 
						select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY3=as.numeric(subset(tmp_student_data, 
						select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
				SGP_Targets=list(CUKU=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
					CUKU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]], 
					MUSU=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
					MUSU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]]),
				SGP_Scale_Score_Targets_TEXT=list(CUKU=list(
					NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
					NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
					NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
							MUSU=list(
					NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
					NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
					NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
							CUKU_Current=list(
					NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
					NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
					NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]])),
							MUSU_Current=list(
					NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
					NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
					NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]]))),
				SGP_Scale_Score_Targets=list(CUKU=list(
					NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
					NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
					NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
								MUSU=list(
					NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
					NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
					NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
								CUKU_Current=list(
					NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
					NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
					NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]])),
								MUSU_Current=list(
					NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
					NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
					NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]]))),
				Cutscores=sgPlot.cutscores[[strsplit(tmp_content_areas[vp], "[.]")[[1]][1]]],
				Years=rev(sgPlot.years),
				Report_Parameters=list(Current_Year=last.year, Content_Area=strsplit(tmp_content_areas[vp], "[.]")[[1]][1], Content_Area_Title=tmp_student_data[[paste("CONTENT_AREA_LABELS", last.year, sep=".")]], 
					State=state, SGP_Targets=sgPlot.sgp.targets, Assessment_Transition=sgPlot.linkages, Fan=SGP::SGPstateData[[state]][["SGP_Configuration"]][['sgPlot.fan.condition']]))

			tmp_student_data_JSON <- getJSON(
							tmp.data=tmp.list,
							sgPlot.sgp_object=sgPlot.sgp_object,
							state=state, 
							content_area=strsplit(tmp_content_areas[vp], "[.]")[[1]][1],
							year=last.year,
							years.for.percentile.trajectories=yearIncrement(last.year, c(0,-1)),
							baseline=sgPlot.baseline, 
							data.type="studentGrowthPlot")
			cat(toJSON(tmp_student_data_JSON, pretty=TRUE, na="null", auto_unbox=TRUE), file=file.path(path.to.pdfs, paste(file_name_json, "_", vp, ".json", sep="")))
		}
	} ### END if ("JSON" %in% sgPlot.output.format)


	if ("PDF" %in% sgPlot.output.format) {
		if (!reports.by.student) {
			################################ SCHOOL Report Catalog LaTeX Code ###################################################################################
			if (is.null(sgPlot.front.page) | !is.null(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Include_Front_Page_in_School_Catalog"]])) {
			cat(paste("\\pdfbookmark[2]{", paste(LAST_NAME, ", ", FIRST_NAME, " (", student_number, ")", sep=""), "}{", n , "}\n\\includepdf[pages={1-}, fitpaper=true]{", 
				path.to.pdfs, "/", file_name, "}\n", sep=""), file=paste("school_catalog_", i, "_", j, ".tex", sep=""), append=TRUE)
			} else {
				cat(paste("\\includepdf[pages={1-},fitpaper=true]{", sgPlot.front.page.ij, "}\n\\pdfbookmark[2]{", paste(LAST_NAME, ", ", FIRST_NAME, " (", 
				student_number, ")", sep=""), "}{", n , "}\n\\includepdf[pages={1-},fitpaper=true]{", path.to.pdfs, "/", file_name, "}\n", sep=""), 
					file=paste("school_catalog_", i, "_", j, ".tex", sep=""), append=TRUE)
			}
			#####################################################################################################################################################
		}

		################################ STUDENT Report Front Page Attach LaTeX Code ########################################################################
		cat("\\documentclass[pdftex]{article}\n\\usepackage{hyperref,pdfpages}\n\\hypersetup{%\n", file=paste("student_report_", i, "_", j, ".tex", sep=""))
		cat(paste("pdftitle={", FIRST_NAME, " ", LAST_NAME, " (", student_number, ")", ": ", pretty_year(last.year), " ", 
			tmp.state, " Growth and Achievement Report},\n", sep=""), file=paste("student_report_", i, "_", j, ".tex", sep=""), append=TRUE)
		cat(paste("pdfauthor={", tmp.organization$Name, "/Center for Assessment Inc.},\n", sep=""), file=paste("student_report_", i, "_", j, ".tex", sep=""), append=TRUE)
		cat("pdfcreator={pdfLaTeX},\n", file=paste("student_report_", i, "_", j, ".tex", sep=""), append=TRUE)
		cat(paste("pdfproducer={", tmp.organization$Name, "/Center for Assessment Inc.}}\n", sep=""), file=paste("student_report_", i, "_", j, ".tex", sep=""), append=TRUE) 
		cat("\\begin{document}\n", file=paste("student_report_", i, "_", j, ".tex", sep=""), append=TRUE)
		if (!is.null(sgPlot.front.page)) {
			cat(paste("\\includepdf[pages={1-}, fitpaper=true]{", sgPlot.front.page.ij, "}\n", sep=""), file=paste("student_report_", i, "_", j, ".tex", sep=""), append=TRUE)
		}
		cat(paste("\\includepdf[pages={1-}, fitpaper=true]{", path.to.pdfs, "/", file_name, "}\n", sep=""), 
			file=paste("student_report_", i, "_", j, ".tex", sep=""), append=TRUE)

		cat("\\end{document}", file=paste("student_report_", i, "_", j, ".tex", sep=""), append=TRUE)
		####################################################################################################################################################

		if (!custom.isr.tf) { # boiler plate student reports:

		## Start pdf device
		
		if (num.charts %in% c(1,2)) {
			report.width=11
			report.height=8.5
		}
		if (num.charts==3) {
			report.width=8.5
			report.height=11
		}
		if (num.charts==4) {
			report.width=8.5
			report.height=11
			# report.width=11
			# report.height=8.5
			# report.width=17
			# report.height=11
		}
		if (num.charts==5) {
			report.width=8.5
			report.height=11
			# report.width=11
			# report.height=17
		}
		if (!num.charts %in% 1:5) {
			stop("Individual Student Report Templates currently only available for situations with 1, 2, 3, 4 or 5 content areas.")
		}

		pdf(file.path(path.to.pdfs, file_name), 
			width=report.width, height=report.height, version="1.4", onefile=TRUE)

		########################################################################################################
		###
		### Overall Report viewport creation
		###
		########################################################################################################

		page.count <- 1; page.content_areas <- list(1:num.charts); legend.tf <- FALSE # redefined below if necessary

		if (num.charts==1) {
			report.vp <- viewport(layout = grid.layout(5, 4, widths = unit(c(2.5, 0.1, 8.3, 0.1), rep("inches", 4)), 
				heights = unit(c(0.55, 0.2, 7, 0.25, 0.5), rep("inches", 5))))

			content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=3)
			top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:4)
			bottom.border.vp <- viewport(layout.pos.row=5, layout.pos.col=1:4)
			left.legend.vp <- viewport(layout.pos.row=2:4, layout.pos.col=1)
		}

		if (num.charts==2) {
			report.vp <- viewport(layout = grid.layout(7, 4, widths = unit(c(2.5, 0.1, 8.3, 0.1), rep("inches", 4)), 
				heights = unit(c(0.35, 0.2, 3.55, 0.25, 3.55, 0.2, 0.4), rep("inches", 7))))

			content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=3)
			content_area_2.vp <- viewport(layout.pos.row=5, layout.pos.col=3)
			top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:4)
			bottom.border.vp <- viewport(layout.pos.row=7, layout.pos.col=1:4)
			left.legend.vp <- viewport(layout.pos.row=2:6, layout.pos.col=1)
		}

		if (num.charts==3) {
			report.vp <- viewport(layout = grid.layout(9, 3, widths = unit(c(0.125, 8.3, 0.075), rep("inches", 3)), 
				heights = unit(c(0.35, 0.1, 3.256, 0.14, 3.256, 0.14, 3.256, 0.1, 0.4), rep("inches", 9))))

			content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=2)
			content_area_2.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
			content_area_3.vp <- viewport(layout.pos.row=7, layout.pos.col=2)
			top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:3)
			bottom.border.vp <- viewport(layout.pos.row=9, layout.pos.col=1:3)
		}

		if (num.charts==4) {
			# report.vp <- viewport(layout = grid.layout(7, 3, widths = unit(c(0.125, 8.3, 0.075), rep("inches", 3)), 
			# 	heights = unit(c(0.4, 0.75, 4, 0.5, 4, 0.85, 0.5), rep("inches", 7))))

			# page.content_areas <- list(1:2, 3:4)

			# content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=2)
			# content_area_2.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
			# content_area_3.vp <- viewport(layout.pos.row=3, layout.pos.col=2)
			# content_area_4.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
			# top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:3)
			# bottom.border.vp <- viewport(layout.pos.row=7, layout.pos.col=1:3)

			if (!is.null(Four_Charts <- SGP::SGPstateData[[state]][["Student_Report_Information"]][["Two_Page_Layout"]][["Four_Charts"]])) {
				page.count <- 2
				report.vp <- Four_Charts[["report.vp"]]
				page.content_areas <- Four_Charts[["page.content_areas"]]
				content_area_1.vp <- Four_Charts[["content_area_1.vp"]]
				content_area_2.vp <- Four_Charts[["content_area_2.vp"]]
				content_area_3.vp <- Four_Charts[["content_area_3.vp"]]
				content_area_4.vp <- Four_Charts[["content_area_4.vp"]]
				top.border.vp <- Four_Charts[["top.border.vp"]]
				bottom.border.vp <- Four_Charts[["bottom.border.vp"]]
				if (!is.null(Four_Charts[["left.legend.vp"]])) {left.legend.vp <- Four_Charts[["left.legend.vp"]]; legend.tf <- TRUE}
			} else { ###  17 x 11 Default
				report.vp <- viewport(layout = grid.layout(7, 6, widths = unit(c(2.5, 0.15, 7, 0.2, 7, 0.15), rep("inches", 6)), 
					heights = unit((11/8.5)*c(0.35, 0.2, 3.55, 0.25, 3.55, 0.2, 0.4), rep("inches", 7))))

				content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=3)
				content_area_2.vp <- viewport(layout.pos.row=3, layout.pos.col=5)
				content_area_3.vp <- viewport(layout.pos.row=5, layout.pos.col=3)
				content_area_4.vp <- viewport(layout.pos.row=5, layout.pos.col=5)
				top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:6)
				bottom.border.vp <- viewport(layout.pos.row=7, layout.pos.col=1:6)
				left.legend.vp <- viewport(layout.pos.row=2:6, layout.pos.col=1)
				legend.tf <- TRUE
			}
			###  11 x 8.5 - Panels generally are too small
			# report.vp <- viewport(layout = grid.layout(7, 5, widths = unit(c(0.15, 5.3, 0.1, 5.3, 0.15), rep("inches", 5)), 
			# 	heights = unit(c(0.35, 0.2, 3.55, 0.25, 3.55, 0.2, 0.4), rep("inches", 7))))

			# content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=2)
			# content_area_2.vp <- viewport(layout.pos.row=3, layout.pos.col=4)
			# content_area_3.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
			# content_area_4.vp <- viewport(layout.pos.row=5, layout.pos.col=4)
			# top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:5)
			# bottom.border.vp <- viewport(layout.pos.row=7, layout.pos.col=1:5)
			# #  No Legend
		}

		if (num.charts==5) {
			# report.vp <- viewport(layout = grid.layout(9, 3, widths = unit(c(0.125, 8.3, 0.075), rep("inches", 3)), 
			# 	heights = unit(c(0.35, 0.1, 3.256, 0.14, 3.256, 0.14, 3.256, 0.1, 0.4), rep("inches", 9))))

			# page.content_areas <- list(1:3, 4:5)

			# content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=2)  # 1st Page
			# content_area_2.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
			# content_area_3.vp <- viewport(layout.pos.row=7, layout.pos.col=2)
			# content_area_4.vp <- viewport(layout.pos.row=3, layout.pos.col=2)  # 2nd Page
			# content_area_5.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
			# top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:3)
			# bottom.border.vp <- viewport(layout.pos.row=9, layout.pos.col=1:3)

			if (!is.null(Five_Charts <- SGP::SGPstateData[[state]][["Student_Report_Information"]][["Two_Page_Layout"]][["Five_Charts"]])) {
				page.count <- 2
				report.vp <- Five_Charts[["report.vp"]]
				page.content_areas <- Five_Charts[["page.content_areas"]]
				content_area_1.vp <- Five_Charts[["content_area_1.vp"]]
				content_area_2.vp <- Five_Charts[["content_area_2.vp"]]
				content_area_3.vp <- Five_Charts[["content_area_3.vp"]]
				content_area_4.vp <- Five_Charts[["content_area_4.vp"]]
				content_area_5.vp <- Five_Charts[["content_area_5.vp"]]
				top.border.vp <- Five_Charts[["top.border.vp"]]
				bottom.border.vp <- Five_Charts[["bottom.border.vp"]]
				if (!is.null(Five_Charts[["left.legend.vp"]])) {left.legend.vp <- Five_Charts[["left.legend.vp"]]; legend.tf <- TRUE}
			} else { ###  17 x 11 Default
				report.vp <- viewport(layout = grid.layout(13, 3, widths = unit((11/8.5)*c(0.125, 8.3, 0.075), rep("inches", 3)), 
					heights = unit(c(0.3, 0.10, 3.25, 0.15, 3.25, 0.15, 3.35, 0.15, 3.25, 0.15, 3.25, 0.10, 0.3), rep("inches", 9))))

				content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=2)
				content_area_2.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
				content_area_3.vp <- viewport(layout.pos.row=7, layout.pos.col=2)
				content_area_4.vp <- viewport(layout.pos.row=9, layout.pos.col=2)
				content_area_5.vp <- viewport(layout.pos.row=11, layout.pos.col=2)
				top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:3)
				bottom.border.vp <- viewport(layout.pos.row=13, layout.pos.col=1:3)
			}
		}

		for (pages in 1:page.count) {
		pushViewport(report.vp)

		for (vp in page.content_areas[[pages]]) {
			tmp_student_data <- as.data.frame(tmp_all_student_data[CONTENT_AREA==tmp_content_areas[vp]])
			if (dim(tmp_student_data)[1]==0) tmp_student_data <- as.data.frame(data.table(tmp_student_data)[NA][,c("CONTENT_AREA", paste("CONTENT_AREA_LABELS", last.year, sep=".")) := tmp_content_areas[vp]])
			pushViewport(get(paste("content_area_", vp, ".vp", sep="")))
			studentGrowthPlot(
				Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("SCALE_SCORE", rev(sgPlot.years), sep="."))),
				Plotting_Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("TRANSFORMED_SCALE_SCORE", rev(sgPlot.years), sep="."))),
				Achievement_Levels=as.character(unlist(subset(tmp_student_data, select=paste("ACHIEVEMENT_LEVEL", rev(sgPlot.years), sep=".")))),
				SGP=as.numeric(subset(tmp_student_data, select=paste(my.sgp, rev(sgPlot.years), sep="."))),
				SGP_Levels=as.character(unlist(subset(tmp_student_data, select=paste(my.sgp.level, rev(sgPlot.years), sep=".")))),
				Grades=as.character(subset(tmp_student_data, select=paste("GRADE", rev(sgPlot.years), sep="."))),
				Content_Areas=as.character(subset(tmp_student_data, select=paste("CONTENT_AREA_LABELS", rev(sgPlot.years), sep="."))),
				Cuts=list(
					NY1=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
						grep("YEAR_1", names(tmp_student_data))), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY2=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
						grep("YEAR_2", names(tmp_student_data))), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY3=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
						grep("YEAR_3", names(tmp_student_data))), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
				Plotting_Cuts=list(
					NY1=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY2=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY3=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
				SGP_Targets=list(
					CUKU=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
					CUKU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]], 
					MUSU=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
					MUSU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]]),
				SGP_Scale_Score_Targets=list(
					CUKU=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
					MUSU=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
					CUKU_Current=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]])),
					MUSU_Current=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]]))),
				Plotting_SGP_Scale_Score_Targets=list(
					CUKU=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
					MUSU=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
					CUKU_Current=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]])),
					MUSU_Current=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]]))),
				Cutscores=sgPlot.cutscores[[strsplit(tmp_content_areas[vp], "[.]")[[1]][1]]],
				Years=rev(sgPlot.years),
				Report_Parameters=list(Current_Year=last.year, Content_Area=strsplit(tmp_content_areas[vp], "[.]")[[1]][1], Content_Area_Title=tmp_student_data[[paste("CONTENT_AREA_LABELS", last.year, sep=".")]], 
					State=state, SGP_Targets=sgPlot.sgp.targets, Assessment_Transition=sgPlot.linkages, Fan=SGP::SGPstateData[[state]][["SGP_Configuration"]][['sgPlot.fan.condition']]))
	
			popViewport()
		} ## END loop over tmp_content_areas

		## Top Legend
		if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["sgPlot.use.student.school.name"]])) {
			student_school_name <- sort(unique(as.character(tmp_grade_data[ID==n][[paste("SCHOOL_NAME", last.year, sep=".")]])))[1] # sort to get rid of potential NA values
			if (is.na(student_school_name)) student_school_name <- tmp_school_name
		} else student_school_name <- tmp_school_name

		pushViewport(top.border.vp)
		grid.rect(gp=gpar(fill=sgPlot.header.footer.color, col=sgPlot.header.footer.color))
		grid.text(x=0.025, y=0.5, paste(FIRST_NAME, " ", LAST_NAME, sep=""),
			gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", col="white", cex=1.65), just="left", default.units="native")
		if (page.count > 1) grid.text(x=0.475, y=0.5, paste("Page", pages), gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", col="white", cex=1.65), just="center", default.units="native")
		grid.text(x=0.975, y=0.5, student_school_name, gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", col="white", cex=1.65), just="right", default.units="native")
		popViewport()

		## Bottom Legend
		pushViewport(bottom.border.vp)
		grid.rect(gp=gpar(fill=sgPlot.header.footer.color, col=sgPlot.header.footer.color))
		grid.text(x=0.02, y=0.70, paste("For more information please visit the", tmp.organization$Name, paste("(", tmp.organization$Abbreviation, ")", sep=""),
			"at", tmp.organization$URL, "or contact", tmp.organization$Contact), gp=gpar(cex=0.8, col="white"), default.units="native", just="left")
		copyright.text <- paste("Cooperatively developed by the ", tmp.organization$Abbreviation, " & the Center for Assessment, Inc.", sep="")
		grid.text(x=0.02, y=0.30, paste(copyright.text, " Distributed by the ", tmp.organization$Abbreviation, ".", sep=""), 
			gp=gpar(cex=0.8, col="white"), default.units="native", just="left")

		# grid.text(x=0.995, y=0.18, copyright.text, gp=gpar(col="white", cex=0.45), default.units="native", just="right")
		# grid.text(x=unit(0.992, "native")-convertWidth(grobWidth(textGrob(copyright.text, gp=gpar(cex=0.45))), "native"), y=0.19, "\\co", 
		#	 gp=gpar(col="white", cex=0.55, fontfamily="HersheySymbol"), default.units="native", just="right")
		popViewport()

		## Left Legend (Only with one or two content areas depicted)

		if (num.charts %in% c(1,2) | legend.tf) { #,4

			pushViewport(left.legend.vp)
	
			# Interpretation
			interpretation.y <- 0.93
			achievement.level.region.colors <- paste("grey", round(seq(62, 91, length=number.achievement.level.regions)), sep="")
	
			grid.roundrect(x=unit(0.5, "native"), y=unit(interpretation.y, "native"), width=unit(0.9, "native"), height=unit(0.06, "native"), 
				gp=gpar(fill=sgPlot.header.footer.color, col="black"))
			grid.text(x=0.5, y=interpretation.y+0.012, "How to interpret this student", gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", cex=1.05, col="white"))
			grid.text(x=0.5, y=interpretation.y-0.012, "growth & achievement report", gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", cex=1.05, col="white"))
	
			grid.roundrect(x=unit(0.2, "native"), y=unit(interpretation.y-0.08, "native"), width=unit(0.1, "native"), height=unit(0.05, "native"), r=unit(0.02, "inches"), 
				gp=gpar(fill=achievement.level.region.colors[1], lwd=1))
			grid.circle(x=0.2, y=interpretation.y-0.08, r=0.02, default.units="native", gp=gpar(fill="white"))
			grid.text(x=0.325, y=interpretation.y-0.0675, tmp.abbreviation, gp=gpar(cex=0.9), default.units="native", just="left")
			grid.text(x=0.325, y=interpretation.y-0.0925, "Scale Score", gp=gpar(cex=0.9), default.units="native", just="left")
	
			tmp.rect.height <- 0.125/number.achievement.level.regions
			for (l in seq(number.achievement.level.regions)) {
				grid.rect(x=unit(0.2, "native"), y=unit(interpretation.y-0.125-(l-1)*tmp.rect.height, "native"), width=unit(0.1, "native"),
					height=unit(tmp.rect.height, "native"), gp=gpar(fill=rev(achievement.level.region.colors)[l], col="white", lwd=1), just=c("center", "top"))
			} 
			grid.roundrect(x=unit(0.2, "native"), y=interpretation.y-0.125, width=unit(0.1, "native"), height=unit(0.125, "native"),
				 r=unit(0.02, "inches"),gp=gpar(col="black", lwd=1.5), just=c("center", "top"))
			grid.text(x=0.325, y=interpretation.y-0.1625, tmp.abbreviation, default.units="native", just="left")
			grid.text(x=0.325, y=interpretation.y-0.1875, "Achievement", default.units="native", just="left")
			grid.text(x=0.325, y=interpretation.y-0.2125, "Levels", default.units="native", just="left")
	
			grid.polygon(x=c(0.1875, 0.1875, 0.17, 0.2, 0.23, 0.2125, 0.2125), y=interpretation.y-c(0.35, 0.30, 0.31, 0.27, 0.31, 0.30, 0.35), 
				default.units="native",gp=gpar(fill="grey50"))
			grid.text(x=0.325, y=interpretation.y-0.285, "Student", gp=gpar(cex=0.9), default.units="native", just="left")
			grid.text(x=0.325, y=interpretation.y-0.31, "Growth", gp=gpar(cex=0.9), default.units="native", just="left")
			grid.text(x=0.325, y=interpretation.y-0.335, "Percentile", gp=gpar(cex=0.9), default.units="native", just="left")
	
			if (!is.null(sgPlot.sgp.targets)) {
				grid.circle(x=0.2, y=interpretation.y-0.41, r=unit(c(0.08, 0.065, 0.045, 0.0225), "inches"),
					gp=gpar(col=c("black", "white", "black", "white"), lwd=0.01, fill=c("black", "white", "black", "white")), default.units="native")
				grid.text(x=0.325, y=interpretation.y-0.385, "Catch Up/Keep Up", gp=gpar(cex=0.9), default.units="native", just="left")
				grid.text(x=0.325, y=interpretation.y-0.41, "Move Up/Stay Up", gp=gpar(cex=0.9), default.units="native", just="left")
				grid.text(x=0.325, y=interpretation.y-0.435, "Targets", gp=gpar(cex=0.9), default.units="native", just="left")
			}
	
			# Suggested uses
			if (is.null(sgPlot.sgp.targets)) suggested.y <- 0.52 else suggested.y <- 0.42
	
			grid.roundrect(x=unit(0.5, "native"), y=unit(suggested.y, "native"), width=unit(0.9, "native"), height=unit(0.06, "native"), 
				gp=gpar(fill=sgPlot.header.footer.color, col="black"))
			grid.text(x=0.5, y=suggested.y, "Suggested Uses", gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", cex=1.05, col="white"))
	
			grid.circle(x=0.075, y=suggested.y-0.07, r=0.01, gp=gpar(fill="black"), default.units="native")
			grid.text(x=0.12, y=suggested.y-0.07, "Review past growth to assess", gp=gpar(cex=0.8), default.units="native", just="left")
			grid.text(x=0.12, y=suggested.y-0.09, "student academic progress toward", gp=gpar(cex=0.8), default.units="native", just="left")
			grid.text(x=0.12, y=suggested.y-0.11, paste(tmp.abbreviation, "achievement goals."), gp=gpar(cex=0.8), default.units="native", just="left")
	
			grid.circle(x=0.075, y=suggested.y-0.14, r=0.01, gp=gpar(fill="black"), default.units="native")
			grid.text(x=0.12, y=suggested.y-0.14, "Develop remediation or enrich-", gp=gpar(cex=0.8), default.units="native", just="left")
			grid.text(x=0.12, y=suggested.y-0.16, "ment plans based on rate of", gp=gpar(cex=0.8), default.units="native", just="left")
			grid.text(x=0.12, y=suggested.y-0.18, "growth needed to reach higher", gp=gpar(cex=0.8), default.units="native", just="left")
			grid.text(x=0.12, y=suggested.y-0.20, paste(tmp.abbreviation, "achievement levels."), gp=gpar(cex=0.8), default.units="native", just="left")
	
			if (sgPlot.fan) {
				grid.circle(x=0.075, y=suggested.y-0.23, r=0.01, gp=gpar(fill="black"), default.units="native")
				grid.text(x=0.12, y=suggested.y-0.23, "Identify the rate of progress", gp=gpar(cex=0.8), default.units="native", just="left")
				grid.text(x=0.12, y=suggested.y-0.25, "needed in order to reach or", gp=gpar(cex=0.8), default.units="native", just="left")
				grid.text(x=0.12, y=suggested.y-0.27, "maintain proficient status", gp=gpar(cex=0.8), default.units="native", just="left")
				grid.text(x=0.12, y=suggested.y-0.29, paste("on the", tmp.abbreviation, "next year."), gp=gpar(cex=0.8), default.units="native", just="left")
			}
	
			# Extra stuff	
			grid.lines(x=1.0, y=c(0.025,0.975), gp=gpar(lwd=1.8), default.units="native")	
			popViewport()
			
		} ## END Left legend

		if (page.count > 1 & pages < max(page.count))	grid.newpage()
		} ## END loop over pages

		## Turn pdf device off
		dev.off()
		} else { 
		#############################
		### Custom Student reports
		#############################
		
			make.custom.isr <- as.function(custom.isr$Custom_ISR_Function$value)
			environment(make.custom.isr) <- environment()
	
			## Turn pdf device ON and OFF within make.custom.isr
			make.custom.isr()
		
		} # End Custom Student reports

		## Code to LaTeX document attaching first page/adding meta-data

		system(paste("pdflatex -interaction=batchmode student_report_", i, "_", j, ".tex", sep=""), ignore.stdout = TRUE)
		file.rename(paste("student_report_", i, "_", j, ".pdf", sep=""), paste(path.to.pdfs, "/", paste(head(unlist(strsplit(file_name, "_")), -1), collapse="_"), ".pdf", sep=""))
	} ### END "PDF" %in% sgPlot.output.format
	
	if (any(c("PNG", "PDF_PIECES") %in% sgPlot.output.format)) {
		report.png.vp <- viewport(width = unit(8.1, "inches"), height = unit(3.64, "inches"))

		for (vp in seq_along(tmp_content_areas)) {
			if ("PNG" %in% sgPlot.output.format) {
				Cairo(file.path(path.to.pdfs, paste(paste(n, vp, sep="_"), "png", sep=".")), width=8.2, height=3.74, units="in", dpi=144, pointsize=10.5, bg="transparent")
			} else {
				pdf(file.path(path.to.pdfs, paste(paste(n, vp, sep="_"), "pdf", sep=".")), width=8.2, height=3.74, version="1.4") 
			}

			tmp_student_data <- as.data.frame(tmp_all_student_data[CONTENT_AREA==tmp_content_areas[vp]])
			if (dim(tmp_student_data)[1]==0) tmp_student_data <- as.data.frame(data.table(tmp_student_data)[NA][,CONTENT_AREA := tmp_content_areas[vp]])
			pushViewport(report.png.vp)
			studentGrowthPlot(
				Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("SCALE_SCORE", rev(sgPlot.years), sep="."))),
				Plotting_Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("TRANSFORMED_SCALE_SCORE", rev(sgPlot.years), sep="."))),
				Achievement_Levels=as.character(unlist(subset(tmp_student_data, select=paste("ACHIEVEMENT_LEVEL", rev(sgPlot.years), sep=".")))),
				SGP=as.numeric(subset(tmp_student_data, select=paste(my.sgp, rev(sgPlot.years), sep="."))),
				SGP_Levels=as.character(unlist(subset(tmp_student_data, select=paste(my.sgp.level, rev(sgPlot.years), sep=".")))),
				Grades=as.character(subset(tmp_student_data, select=paste("GRADE", rev(sgPlot.years), sep="."))),
				Content_Areas=as.character(subset(tmp_student_data, select=paste("CONTENT_AREA_LABELS", rev(sgPlot.years), sep="."))),
				Cuts=list(
					NY1=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
						grep("YEAR_1", names(tmp_student_data))), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY2=as.numeric(subset(tmp_student_data,  select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
						grep("YEAR_2", names(tmp_student_data))), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY3=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
						grep("YEAR_3", names(tmp_student_data))), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
				Plotting_Cuts=list(
					NY1=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY2=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
					NY3=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
				SGP_Targets=list(
					CUKU=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
					CUKU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]], 
					MUSU=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
					MUSU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]]),
				SGP_Scale_Score_Targets=list(
					CUKU=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
					MUSU=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
					CUKU_Current=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]])),
					MUSU_Current=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]]))),
				Plotting_SGP_Scale_Score_Targets=list(
					CUKU=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
					MUSU=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
					CUKU_Current=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]])),
					MUSU_Current=list(
						NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
						NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
						NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]]))),
				Cutscores=sgPlot.cutscores[[strsplit(tmp_content_areas[vp], "[.]")[[1]][1]]],
				Years=rev(sgPlot.years),
				Report_Parameters=list(Current_Year=last.year, Content_Area=strsplit(tmp_content_areas[vp], "[.]")[[1]][1], Content_Area_Title=tmp_student_data[[paste("CONTENT_AREA_LABELS", last.year, sep=".")]], State=state, SGP_Targets=sgPlot.sgp.targets,
					Assessment_Transition=sgPlot.linkages, Fan=SGP::SGPstateData[[state]][["SGP_Configuration"]][['sgPlot.fan.condition']]))
			popViewport()
			dev.off()
		} ## END loop over tmp_content_areas
	} ### END if ("PNG" %in% sgPlot.output.format)
	
	} ## END for loop for STUDENTS (n)
	} ## END for loop for GRADES (k)
	if (!reports.by.student & "PDF" %in% sgPlot.output.format) {
		cat("\\end{document}", file=paste("school_catalog_", i, "_", j, ".tex", sep=""), append=TRUE)
		system(paste("pdflatex -interaction=batchmode school_catalog_", i, "_", j, ".tex", sep=""), ignore.stdout = TRUE)
		system(paste("pdflatex -interaction=batchmode school_catalog_", i, "_", j, ".tex", sep=""), ignore.stdout = TRUE)
		file.rename(paste("school_catalog_", i, "_", j, ".pdf", sep=""), file.path(sgPlot.folder, year_folder, district_folder, school_folder,
			paste(year_folder, "_", district_folder, "_", school_folder, "_Individual_SGP_Report_Catalog.pdf", sep="")))
	}
	if (sgPlot.cleanup & "PDF" %in% sgPlot.output.format) {
		files.to.remove <- list.files(pattern=paste(i, j, sep="_"), all.files=TRUE)
		lapply(files.to.remove, file.remove)
		files.to.remove <- list.files(path=paste(c(sgPlot.folder, year_folder, district_folder, school_folder), collapse=.Platform$file.sep), 
			pattern="REPORT", all.files=TRUE, full.names=TRUE, recursive=TRUE)
		lapply(files.to.remove, file.remove)
	}

	if (identical(.Platform$OS.type, "unix") & sgPlot.zip & !reports.by.student) {
		if (1000*as.numeric(unlist(strsplit(system(paste("du -s", file.path(sgPlot.folder, year_folder, district_folder, school_folder)), intern=TRUE), "\t"))[1]) > 4000000000) {
			tmp.working.directory <- getwd()
			setwd(file.path(sgPlot.folder, year_folder, district_folder))
			if (paste(school_folder, ".tar.gz", sep="") %in% list.files()) file.remove(paste(school_folder, ".tar.gz", sep=""))
			system(paste("tar cfz", paste(school_folder, ".tar.gz", sep=""), school_folder, sep=" "))
			setwd(tmp.working.directory)
		} else {
			tmp.working.directory <- getwd()
			setwd(file.path(sgPlot.folder, year_folder, district_folder))
			if (paste(school_folder, ".zip", sep="") %in% list.files()) file.remove(paste(school_folder, ".zip", sep=""))
			suppressMessages(
				zip(paste(school_folder, ".zip", sep=""), school_folder, flags="-rmq")
			)
			setwd(tmp.working.directory)
		}
		unlink(file.path(sgPlot.folder, year_folder, district_folder, school_folder), recursive=TRUE)
	}

	message(paste("\tStarted", last.year, tmp_school_name, "student growth plots:", started.date))
	message(paste("\tFinished", last.year, tmp_school_name, "student growth plots:", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END for loop for SCHOOLS (j)
	} ## END for loop for DISTRICTS (i)

	return("DONE")
} ### END if (!reports.by.instructor) 


#############################################################################################################################
###
### Reports assuming DISTRICT, SCHOOLS, INSTRUCTORS, GRADES, and then STUDENTS
###
#############################################################################################################################

if (reports.by.instructor) {

	### Extend sgPlot.folder
	sgPlot.folder <- file.path(sgPlot.folder, "Instructor")

	### Loop over unique DISTRICTS, SCHOOLS, INSTRUCTORS, GRADES and then STUDENTS

	tmp.keys <- paste(c("DISTRICT_NUMBER", "SCHOOL_NUMBER", "INSTRUCTOR_NUMBER", "GRADE", "LAST_NAME", "FIRST_NAME"), last.year, sep=".")

	## Districts
	setkeyv(sgPlot.data, tmp.keys[1])

	for (i in districts) {
		if (sgPlot.demo.report | identical(i, -999L)) {
			tmp_district_name <- "Sample District"
			district_folder <- "Sample_District"
		} else {
			if (sgPlot.folder.names=="name") {
				tmp_district_name <- as.character(sgPlot.data[list(i)][[paste("DISTRICT_NAME", last.year, sep=".")]][1])
					district_folder <- gsub(" ", "_", paste(tmp_district_name, i))
			} else {
				tmp_district_name <- as.character(sgPlot.data[list(i)][[paste("DISTRICT_NAME", last.year, sep=".")]][1])
					district_folder <- as.character(i)
			}
		}
		tmp_district_ids <- unique(sgPlot.data[list(i)]$ID)
		tmp_district_data <- subset(sgPlot.data, ID %in% tmp_district_ids)

	## Schools
	setkeyv(tmp_district_data, tmp.keys[2])

	for (j in intersect(tmp_district_data[[paste("SCHOOL_NUMBER", last.year, sep=".")]], schools)) {

		if (sgPlot.demo.report | identical(j, -99L)) {
			tmp_school_name <- "Sample School"
			school_folder <- "Sample_School"
		} else {
			if (sgPlot.folder.names=="name") {
				tmp_school_name <- as.character(tmp_district_data[list(j)][[paste("SCHOOL_NAME", last.year, sep=".")]][1])
				school_folder <- gsub(" ", "_", paste(tmp_school_name, j))
			} else {
				tmp_school_name <- as.character(tmp_district_data[list(j)][[paste("SCHOOL_NAME", last.year, sep=".")]][1])
				school_folder <- as.character(j)
			}
		}
		tmp_school_ids <- unique(tmp_district_data[list(j)]$ID)
		tmp_school_data <- subset(tmp_district_data, ID %in% tmp_school_ids)

	## Instructor
	setkeyv(tmp_school_data, tmp.keys[2])
	instructors <- sort(unique(unlist(tmp_school_data[list(j), tmp.keys[3], with=FALSE]))) %w/o% NA
	setkeyv(tmp_school_data, tmp.keys[3])

	for (k in instructors) {

	started.at <- proc.time()
	started.date <- date()

	if (!is.null(sgPlot.front.page)) {
		sgPlot.front.page.ijk <- paste(paste(unlist(strsplit(sgPlot.front.page, "[.]"))[1], i, j, k, sep="_"), ".pdf", sep="")
		file.copy(sgPlot.front.page, sgPlot.front.page.ijk)
	}

	if (sgPlot.demo.report | identical(k, -9L)) {
		tmp_instructor_name <- "Sample Instructor"
		instructor_folder <- "Sample Instructor"
		FIRST_NAME_TEACHER <- "Sample"
		LAST_NAME_TEACHER <- "Teacher"
	} else {
		FIRST_NAME_TEACHER <- unlist(strsplit(as.character(tmp_school_data[list(k)][[paste("INSTRUCTOR_NAME", last.year, sep=".")]][1]), ", "))[2]
		LAST_NAME_TEACHER <- unlist(strsplit(as.character(tmp_school_data[list(k)][[paste("INSTRUCTOR_NAME", last.year, sep=".")]][1]), ", "))[1]
		if (sgPlot.folder.names=="name") {
			tmp_instructor_name <- as.character(tmp_school_data[list(k)][[paste("INSTRUCTOR_NAME", last.year, sep=".")]][1])
			instructor_folder <- gsub(" ", "_", paste(tmp_instructor_name, k))
		} else {
			tmp_instructor_name <- as.character(tmp_school_data[list(k)][[paste("INSTRUCTOR_NUMBER", last.year, sep=".")]][1])
			instructor_folder <- as.character(k)
		}
	}
	if ("PDF" %in% sgPlot.output.format) {
		######################## INSTRUCTOR Report Catalog LaTeX Header #################################################################################
		cat("\\documentclass[pdftex]{book}\n\\usepackage{hyperref,pdfpages}\n\\hypersetup{%\n", file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""))
		cat(paste("pdftitle={", tmp_instructor_name, ": ", pretty_year(last.year), " ", tmp.state, " Growth and Achievement Reports},\n", sep=""), 
			file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
		cat(paste("pdfauthor={", tmp.organization$Name, "/Center for Assessment Inc.},\n", sep=""), "pdfcreator={pdfLaTeX},\n", 
		paste("pdfproducer={", tmp.organization$Name, "/Center for Assessment Inc.}}\n", sep=""), 
		"\\begin{document}\n", sep="", file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
		cat(paste("\\pdfbookmark[-1]{", tmp_district_name, "}{", i, "}\n", sep=""), file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
		cat(paste("\\pdfbookmark[0]{", tmp_school_name, "}{", j, "}\n", sep=""), file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
		cat(paste("\\pdfbookmark[1]{", tmp_instructor_name, "}{", k, "}\n", sep=""), file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
		##############################################################################################################################################
	}
	tmp_instructor_ids <- unique(tmp_school_data[list(k)]$ID)
	tmp_instructor_data <- subset(tmp_school_data[list(k)], ID %in% tmp_instructor_ids) ### NOTE: possible that students have multiple instructors so need list(k)

	## Grades
	setkeyv(tmp_instructor_data, tmp.keys[3])
	grades <- sort(unique(unlist(tmp_instructor_data[list(k), tmp.keys[4], with=FALSE]))) %w/o% NA
	setkeyv(tmp_instructor_data, tmp.keys[4])

	for (l in grades) {
		if (sgPlot.folder.names=="name") {
			grade_folder <- paste("Grade", l, sep="_")
		} else {
			grade_folder <- substr(paste("0", as.character(l), sep=""), nchar(l), nchar(l)+1)
		}
	
		if ("PDF" %in% sgPlot.output.format) {
			################################ INSTRUCTOR Report Catalog LaTeX Code #########################################################
			cat(paste("\\pdfbookmark[2]{Grade ", l, "}{", j, k, l, "}\n", sep=""), 
				file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE) ## NOTE: j, k, l included in anchor for uniqueness
			###########################################################################################################################
		}
	
		tmp_grade_ids <- unique(tmp_instructor_data[list(l)]$ID)
		tmp_grade_data <- subset(tmp_instructor_data, ID %in% tmp_grade_ids)
	
		### Create path to pdf files
	
		path.to.pdfs <- file.path(sgPlot.folder, year_folder, district_folder, school_folder, instructor_folder, grade_folder)
		dir.create(path.to.pdfs, recursive=TRUE, showWarnings=FALSE)
	
		## Students
	
		setkeyv(tmp_grade_data, tmp.keys[5:6])
	
		for (n in unique(tmp_grade_data[["ID"]])) {
			FIRST_NAME <- gsub(" |/", "-", sort(tmp_grade_data[ID==n][[tmp.keys[6]]])[1]) 
			LAST_NAME <- gsub(" |/", "-", sort(tmp_grade_data[ID==n][[tmp.keys[5]]])[1])
			if (sgPlot.anonymize) {
				student_number <- 1234567890
			} else {
				student_number <- n
			}
			if (sgPlot.folder.names=="name" | sgPlot.anonymize) {
				file_name <- paste(paste(FIRST_NAME, LAST_NAME, student_number, year_folder, "REPORT", sep="_"), ".pdf", sep="")
			} else {
				file_name <- paste(n, "_REPORT.pdf", sep="")
			}
	
		if ("PDF" %in% sgPlot.output.format) {
			################################ INSTRUCTOR Report Catalog LaTeX Code ###################################################################################
			if (is.null(sgPlot.front.page)) {
				cat(paste("\\pdfbookmark[3]{", paste(LAST_NAME, ", ", FIRST_NAME, " (", student_number, ")", sep=""), "}{", n , "}\n\\includepdf[pages={1-}, fitpaper=true]{", 
					path.to.pdfs, "/", file_name, "}\n", sep=""), file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
			} else {
				cat(paste("\\includepdf[pages={1-}, fitpaper=true]{", sgPlot.front.page.ijk, "}\n\\pdfbookmark[3]{", paste(LAST_NAME, ", ", FIRST_NAME, " (", 
				student_number, ")", sep=""), "}{", n , "}\n\\includepdf[pages={1-}, fitpaper=true]{", path.to.pdfs, "/", file_name, "}\n", sep=""), 
					file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
			}
			#####################################################################################################################################################
	
			################################ STUDENT Report Front Page Attach LaTeX Code ########################################################################
			cat("\\documentclass[pdftex]{article}\n\\usepackage{hyperref,pdfpages}\n\\hypersetup{%\n", file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""))
			cat(paste("pdftitle={", FIRST_NAME, " ", LAST_NAME, " (", student_number, ")", ": ", pretty_year(last.year), " ", 
				tmp.state, " Growth and Achievement Report},\n", sep=""), file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
			cat(paste("pdfauthor={", tmp.organization$Name, "/Center for Assessment Inc.},\n", sep=""), file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
			cat("pdfcreator={pdfLaTeX},\n", file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
			cat(paste("pdfproducer={", tmp.organization$Name, "/Center for Assessment Inc.}}\n", sep=""), file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE) 
			cat("\\begin{document}\n", file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
			if (!is.null(sgPlot.front.page)) {
				cat(paste("\\includepdf[pages={1-}, fitpaper=true]{", sgPlot.front.page.ijk, "}\n", sep=""), file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
			}
			cat(paste("\\includepdf[pages={1-}, fitpaper=true]{", path.to.pdfs, "/", file_name, "}\n", sep=""), 
				file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
	
			cat("\\end{document}", file=paste("student_report_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
			####################################################################################################################################################
	
			## Start pdf device		
			if (num.charts %in% c(1,2)) {
				report.width=11
				report.height=8.5
			}
			if (num.charts==3) {
				report.width=8.5
				report.height=11
			}
			if (num.charts==4) {
				report.width=17
				report.height=11
			}
			if (num.charts==5) {
				report.width=11
				report.height=17
			}
			if (!num.charts %in% 1:5) {
				stop("Individual Student Report Templates currently only available for situations with 1, 2, 3, 4 or 5 content areas.")
			}
	
			pdf(file.path(path.to.pdfs, file_name), 
				  width=report.width, height=report.height, version="1.4")
	
			########################################################################################################
			###
			### Overall Report viewport creation
			###
			########################################################################################################
	
			if (num.charts==1) {
				report.vp <- viewport(layout = grid.layout(5, 4, widths = unit(c(2.5, 0.1, 8.3, 0.1), rep("inches", 4)), 
					heights = unit(c(0.55, 0.2, 7, 0.25, 0.5), rep("inches", 5))))
	
				content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=3)
				top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:4)
				bottom.border.vp <- viewport(layout.pos.row=5, layout.pos.col=1:4)
				left.legend.vp <- viewport(layout.pos.row=2:4, layout.pos.col=1)
			}
	
			if (num.charts==2) {
				report.vp <- viewport(layout = grid.layout(7, 4, widths = unit(c(2.5, 0.1, 8.3, 0.1), rep("inches", 4)), 
					heights = unit(c(0.35, 0.2, 3.55, 0.25, 3.55, 0.2, 0.4), rep("inches", 7))))
	
				content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=3)
				content_area_2.vp <- viewport(layout.pos.row=5, layout.pos.col=3)
				top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:4)
				bottom.border.vp <- viewport(layout.pos.row=7, layout.pos.col=1:4)
				left.legend.vp <- viewport(layout.pos.row=2:6, layout.pos.col=1)
			}
	
			if (num.charts==3) {
				report.vp <- viewport(layout = grid.layout(9, 3, widths = unit(c(0.125, 8.3, 0.075), rep("inches", 3)), 
					heights = unit(c(0.35, 0.1, 3.256, 0.14, 3.256, 0.14, 3.256, 0.1, 0.4), rep("inches", 9))))
	
				content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=2)
				content_area_2.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
				content_area_3.vp <- viewport(layout.pos.row=7, layout.pos.col=2)
				top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:3)
				bottom.border.vp <- viewport(layout.pos.row=9, layout.pos.col=1:3)
			}
	
			if (num.charts==4) {
				report.vp <- viewport(layout = grid.layout(7, 6, widths = unit(c(2.5, 0.15, 7, 0.2, 7, 0.15), rep("inches", 6)), 
					heights = unit((11/8.5)*c(0.35, 0.2, 3.55, 0.25, 3.55, 0.2, 0.4), rep("inches", 7))))
	
				content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=3)
				content_area_2.vp <- viewport(layout.pos.row=3, layout.pos.col=5)
				content_area_3.vp <- viewport(layout.pos.row=5, layout.pos.col=3)
				content_area_4.vp <- viewport(layout.pos.row=5, layout.pos.col=5)
				top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:6)
				bottom.border.vp <- viewport(layout.pos.row=7, layout.pos.col=1:6)
				left.legend.vp <- viewport(layout.pos.row=2:6, layout.pos.col=1)
			}
	
			if (num.charts==5) {
				report.vp <- viewport(layout = grid.layout(13, 3, widths = unit((11/8.5)*c(0.125, 8.3, 0.075), rep("inches", 3)), 
					heights = unit(c(0.3, 0.10, 3.25, 0.15, 3.25, 0.15, 3.35, 0.15, 3.25, 0.15, 3.25, 0.10, 0.3), rep("inches", 9))))
	
				content_area_1.vp <- viewport(layout.pos.row=3, layout.pos.col=2)
				content_area_2.vp <- viewport(layout.pos.row=5, layout.pos.col=2)
				content_area_3.vp <- viewport(layout.pos.row=7, layout.pos.col=2)
				content_area_4.vp <- viewport(layout.pos.row=9, layout.pos.col=2)
				content_area_5.vp <- viewport(layout.pos.row=11, layout.pos.col=2)
				top.border.vp <- viewport(layout.pos.row=1, layout.pos.col=1:3)
				bottom.border.vp <- viewport(layout.pos.row=13, layout.pos.col=1:3)
			}
	
			pushViewport(report.vp)
	
			for (vp in seq_along(tmp_content_areas)) {
				tmp_student_data <- as.data.frame(tmp_all_student_data[CONTENT_AREA==tmp_content_areas[vp]])
				if (dim(tmp_student_data)[1]==0) tmp_student_data <- as.data.frame(data.table(tmp_student_data)[NA][,CONTENT_AREA := tmp_content_areas[vp]])
				pushViewport(get(paste("content_area_", vp, ".vp", sep="")))
				studentGrowthPlot(
					Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("SCALE_SCORE", rev(sgPlot.years), sep="."))),
					Plotting_Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("TRANSFORMED_SCALE_SCORE", rev(sgPlot.years), sep="."))),
					Achievement_Levels=as.character(unlist(subset(tmp_student_data, select=paste("ACHIEVEMENT_LEVEL", rev(sgPlot.years), sep=".")))),
					SGP=as.numeric(subset(tmp_student_data, select=paste(my.sgp, rev(sgPlot.years), sep="."))),
					SGP_Levels=as.character(unlist(subset(tmp_student_data, select=paste(my.sgp.level, rev(sgPlot.years), sep=".")))),
					Grades=as.character(subset(tmp_student_data, select=paste("GRADE", rev(sgPlot.years), sep="."))),
					Content_Areas=as.character(subset(tmp_student_data, select=paste("CONTENT_AREA_LABELS", rev(sgPlot.years), sep="."))),
					Cuts=list(
						NY1=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_1", names(tmp_student_data))), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY2=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_2", names(tmp_student_data))), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY3=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_3", names(tmp_student_data))), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
					Plotting_Cuts=list(
						NY1=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY2=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY3=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
					SGP_Targets=list(
						CUKU=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
						CUKU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]], 
						MUSU=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
						MUSU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]]),
					SGP_Scale_Score_Targets=list(
						CUKU=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
						MUSU=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
						CUKU_Current=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]])),
						MUSU_Current=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]]))),
					Plotting_SGP_Scale_Score_Targets=list(
						CUKU=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
						MUSU=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
						CUKU_Current=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]])),
						MUSU_Current=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]]))),
					Cutscores=sgPlot.cutscores[[strsplit(tmp_content_areas[vp], "[.]")[[1]][1]]],
					Years=rev(sgPlot.years),
					Report_Parameters=list(Current_Year=last.year, Content_Area=strsplit(tmp_content_areas[vp], "[.]")[[1]][1],
						Content_Area_Title=tmp_student_data[[paste("CONTENT_AREA_LABELS", last.year, sep=".")]],
						State=state, Denote_Content_Area=tmp_student_data[['CONTENT_AREA_RESPONSIBILITY']]=="Content Area Responsibility: Yes", SGP_Targets=sgPlot.sgp.targets,
						Assessment_Transition=sgPlot.linkages, Fan=SGP::SGPstateData[[state]][["SGP_Configuration"]][['sgPlot.fan.condition']]))
			popViewport()
	
			} ## END loop over tmp_content_areas
	
			## Top Legend
			if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["sgPlot.use.student.school.name"]])) {
				student_school_name <- sort(unique(tmp_student_data[[paste("SCHOOL_NAME", last.year, sep=".")]]))[1] # sort to get rid of potential NA values
			} else student_school_name <- tmp_school_name

			pushViewport(top.border.vp)
			grid.rect(gp=gpar(fill=sgPlot.header.footer.color, col=sgPlot.header.footer.color))
			grid.text(x=0.025, y=0.5, paste(FIRST_NAME, " ", LAST_NAME, sep=""), 
				gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", col="white", cex=1.65), just="left", default.units="native")
			grid.text(x=0.975, y=0.5, paste(student_school_name, ": ", FIRST_NAME_TEACHER, " ", LAST_NAME_TEACHER, sep=""), 
				gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", col="white", cex=1.65), just="right", default.units="native")
			popViewport()
	
	
			## Bottom Legend
			pushViewport(bottom.border.vp)
			grid.rect(gp=gpar(fill=sgPlot.header.footer.color, col=sgPlot.header.footer.color))
			grid.text(x=0.02, y=0.70, paste("For more information please visit the", tmp.organization$Name, paste("(", tmp.organization$Abbreviation, ")", sep=""), 
				"at", tmp.organization$URL, "or contact", tmp.organization$Contact), gp=gpar(cex=0.8, col="white"), default.units="native", just="left")
			copyright.text <- paste("Cooperatively developed by the ", tmp.organization$Abbreviation, " & the Center for Assessment, Inc.", sep="")
			grid.text(x=0.02, y=0.30, paste(copyright.text, " Distributed by the ", tmp.organization$Abbreviation, ".", sep=""), 
				gp=gpar(cex=0.8, col="white"), default.units="native", just="left")
	
			popViewport()
	
			## Left Legend (Only with one or two content areas depicted)
	
			if (num.charts %in% c(1,2,4)) {
				pushViewport(left.legend.vp)
		
				# Interpretation
				interpretation.y <- 0.93
				achievement.level.region.colors <- paste("grey", round(seq(62, 91, length=number.achievement.level.regions)), sep="")
	
				grid.roundrect(x=unit(0.5, "native"), y=unit(interpretation.y, "native"), width=unit(0.9, "native"), height=unit(0.06, "native"), 
					gp=gpar(fill=sgPlot.header.footer.color, col="black"))
				grid.text(x=0.5, y=interpretation.y+0.012, "How to interpret this student", gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", cex=1.05, col="white"))
				grid.text(x=0.5, y=interpretation.y-0.012, "growth & achievement report", gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", cex=1.05, col="white"))
	
				grid.roundrect(x=unit(0.2, "native"), y=unit(interpretation.y-0.08, "native"), width=unit(0.1, "native"), height=unit(0.05, "native"), r=unit(0.02, "inches"), 
					gp=gpar(fill=achievement.level.region.colors[1], lwd=1))
				grid.circle(x=0.2, y=interpretation.y-0.08, r=0.02, default.units="native", gp=gpar(fill="white"))
				grid.text(x=0.325, y=interpretation.y-0.0675, tmp.abbreviation, gp=gpar(cex=0.9), default.units="native", just="left")
				grid.text(x=0.325, y=interpretation.y-0.0925, "Scale Score", gp=gpar(cex=0.9), default.units="native", just="left")
	
				tmp.rect.height <- 0.125/number.achievement.level.regions
				for (p in seq(number.achievement.level.regions)) {
					grid.rect(x=unit(0.2, "native"), y=unit(interpretation.y-0.125-(p-1)*tmp.rect.height, "native"), width=unit(0.1, "native"), height=unit(tmp.rect.height, "native"),
					gp=gpar(fill=rev(achievement.level.region.colors)[p], col="white", lwd=1), just=c("center", "top"))
				} 
				grid.roundrect(x=unit(0.2, "native"), y=interpretation.y-0.125, width=unit(0.1, "native"), height=unit(0.125, "native"), r=unit(0.02, "inches"),
					gp=gpar(col="black", lwd=1.5), just=c("center", "top"))
				grid.text(x=0.325, y=interpretation.y-0.1625, tmp.abbreviation, default.units="native", just="left")
				grid.text(x=0.325, y=interpretation.y-0.1875, "Achievement", default.units="native", just="left")
				grid.text(x=0.325, y=interpretation.y-0.2125, "Levels", default.units="native", just="left")
	
				grid.polygon(x=c(0.1875, 0.1875, 0.17, 0.2, 0.23, 0.2125, 0.2125), y=interpretation.y-c(0.35, 0.30, 0.31, 0.27, 0.31, 0.30, 0.35), default.units="native",
					gp=gpar(fill="grey50"))
				grid.text(x=0.325, y=interpretation.y-0.285, "Student", gp=gpar(cex=0.9), default.units="native", just="left")
				grid.text(x=0.325, y=interpretation.y-0.31, "Growth", gp=gpar(cex=0.9), default.units="native", just="left")
				grid.text(x=0.325, y=interpretation.y-0.335, "Percentile", gp=gpar(cex=0.9), default.units="native", just="left")
	
				# Suggested uses
				suggested.y <- 0.52
	
				grid.roundrect(x=unit(0.5, "native"), y=unit(suggested.y, "native"), width=unit(0.9, "native"), height=unit(0.06, "native"), 
					gp=gpar(fill=sgPlot.header.footer.color, col="black"))
				grid.text(x=0.5, y=suggested.y, "Suggested Uses", gp=gpar(fontface="bold", fontfamily="Helvetica-Narrow", cex=1.05, col="white"))
	
				grid.circle(x=0.075, y=suggested.y-0.07, r=0.01, gp=gpar(fill="black"), default.units="native")
				grid.text(x=0.12, y=suggested.y-0.07, "Review past growth to assess", gp=gpar(cex=0.8), default.units="native", just="left")
				grid.text(x=0.12, y=suggested.y-0.09, "student academic progress toward", gp=gpar(cex=0.8), default.units="native", just="left")
				grid.text(x=0.12, y=suggested.y-0.11, paste(tmp.abbreviation, "achievement goals."), gp=gpar(cex=0.8), default.units="native", just="left")
	
				grid.circle(x=0.075, y=suggested.y-0.14, r=0.01, gp=gpar(fill="black"), default.units="native")
				grid.text(x=0.12, y=suggested.y-0.14, "Develop remediation or enrich-", gp=gpar(cex=0.8), default.units="native", just="left")
				grid.text(x=0.12, y=suggested.y-0.16, "ment plans based on rate of", gp=gpar(cex=0.8), default.units="native", just="left")
				grid.text(x=0.12, y=suggested.y-0.18, "growth needed to reach higher", gp=gpar(cex=0.8), default.units="native", just="left")
				grid.text(x=0.12, y=suggested.y-0.20, paste(tmp.abbreviation, "achievement levels."), gp=gpar(cex=0.8), default.units="native", just="left")
	
				if (sgPlot.fan) {
					grid.circle(x=0.075, y=suggested.y-0.23, r=0.01, gp=gpar(fill="black"), default.units="native")
					grid.text(x=0.12, y=suggested.y-0.23, "Identify the rate of progress", gp=gpar(cex=0.8), default.units="native", just="left")
					grid.text(x=0.12, y=suggested.y-0.25, "needed in order to reach or", gp=gpar(cex=0.8), default.units="native", just="left")
					grid.text(x=0.12, y=suggested.y-0.27, "maintain proficient status", gp=gpar(cex=0.8), default.units="native", just="left")
					grid.text(x=0.12, y=suggested.y-0.29, paste("on the", tmp.abbreviation, "next year."), gp=gpar(cex=0.8), default.units="native", just="left")
				}
	
				# Extra stuff
				grid.lines(x=1.0, y=c(0.025,0.975), gp=gpar(lwd=1.8), default.units="native")
				popViewport()				
			} ## END Left legend
	
			## Turn pdf device off
			dev.off()
	
			## Code to LaTeX document attaching first page/adding meta-data
			system(paste("pdflatex -interaction=batchmode student_report_", i, "_", j, "_", k, ".tex", sep=""), ignore.stdout = TRUE)
			file.rename(paste("student_report_", i, "_", j, "_", k, ".pdf", sep=""), paste(path.to.pdfs, "/", paste(head(unlist(strsplit(file_name, "_")), -1), collapse="_"), ".pdf", sep=""))
		} ### END "PDF" %in% sgPlot.output.format
	
		if (any(c("PNG", "PDF_PIECES") %in% sgPlot.output.format)) {
			report.png.vp <- viewport(width = unit(8.1, "inches"), height = unit(3.64, "inches"))
	
			for (vp in seq_along(tmp_content_areas)) {
				if ("PNG" %in% sgPlot.output.format) {
					Cairo(file.path(path.to.pdfs, paste(paste(n, vp, sep="_"), "png", sep=".")), width=8.2, height=3.74, units="in", dpi=144, pointsize=10.5, bg="transparent")
				} else {
					pdf(file.path(path.to.pdfs, paste(paste(n, vp, sep="_"), "pdf", sep=".")), width=8.2, height=3.74, version="1.4") 
				}
	
				tmp_student_data <- as.data.frame(tmp_all_student_data[CONTENT_AREA==tmp_content_areas[vp]])
				if (dim(tmp_student_data)[1]==0) tmp_student_data <- as.data.frame(data.table(tmp_student_data)[NA][,CONTENT_AREA := tmp_content_areas[vp]])
				pushViewport(report.png.vp)
				studentGrowthPlot(
					Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("SCALE_SCORE", rev(sgPlot.years), sep="."))),
					Plotting_Scale_Scores=as.numeric(subset(tmp_student_data, select=paste("TRANSFORMED_SCALE_SCORE", rev(sgPlot.years), sep="."))),
					Achievement_Levels=as.character(unlist(subset(tmp_student_data, select=paste("ACHIEVEMENT_LEVEL", rev(sgPlot.years), sep=".")))),
					SGP=as.numeric(subset(tmp_student_data, select=paste(my.sgp, rev(sgPlot.years), sep="."))),
					SGP_Levels=as.character(unlist(subset(tmp_student_data, select=paste(my.sgp.level, rev(sgPlot.years), sep=".")))),
					Grades=as.character(subset(tmp_student_data, select=paste("GRADE", rev(sgPlot.years), sep="."))),
					Content_Areas=as.character(subset(tmp_student_data, select=paste("CONTENT_AREA_LABELS", rev(sgPlot.years), sep="."))),
					Cuts=list(
						NY1=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_1", names(tmp_student_data))), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY2=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_2", names(tmp_student_data))), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY3=as.numeric(subset(tmp_student_data, select=setdiff(intersect(grep(trajectory.cuts, names(tmp_student_data)), 
							grep("YEAR_3", names(tmp_student_data))), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
					Plotting_Cuts=list(
						NY1=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_1_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY2=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_2_CURRENT_TRANSFORMED", names(tmp_student_data))))),
						NY3=as.numeric(subset(tmp_student_data, select=intersect(grep(trajectory.cuts, names(tmp_student_data)), grep("YEAR_3_CURRENT_TRANSFORMED", names(tmp_student_data)))))),
					SGP_Targets=list(
						CUKU=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
						CUKU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]], 
						MUSU=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], sep="_"), last.year, sep=".")]], 
						MUSU_Current=tmp_student_data[[paste(paste(my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "CURRENT", sep="_"), last.year, sep=".")]]),
					SGP_Scale_Score_Targets=list(
						CUKU=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
						MUSU=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3", sep="_")]])),
						CUKU_Current=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]])),
						MUSU_Current=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT", sep="_")]]))),
					Plotting_SGP_Scale_Score_Targets=list(
						CUKU=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
						MUSU=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_TRANSFORMED", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_TRANSFORMED", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_TRANSFORMED", sep="_")]])),
						CUKU_Current=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]])),
						MUSU_Current=list(
							NY1=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_1_CURRENT_TRANSFORMED", sep="_")]]),
							NY2=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_2_CURRENT_TRANSFORMED", sep="_")]]),
							NY3=as.numeric(tmp_student_data[[paste('SCALE_SCORE', my.sgp.target.label[1], "MOVE_UP_STAY_UP", my.sgp.target.label[2], "PROJ_YEAR_3_CURRENT_TRANSFORMED", sep="_")]]))),
					Cutscores=sgPlot.cutscores[[strsplit(tmp_content_areas[vp], "[.]")[[1]][1]]],
					Years=rev(sgPlot.years),
					Report_Parameters=list(Current_Year=last.year, Content_Area=strsplit(tmp_content_areas[vp], "[.]")[[1]][1],
						Content_Area_Title=tmp_student_data[[paste("CONTENT_AREA_LABELS", last.year, sep=".")]],
						State=state, Denote_Content_Area=tmp_student_data[['CONTENT_AREA_RESPONSIBILITY']]=="Content Area Responsibility: Yes", SGP_Targets=sgPlot.sgp.targets,
						Assessment_Transition=sgPlot.linkages, Fan=SGP::SGPstateData[[state]][["SGP_Configuration"]][['sgPlot.fan.condition']]))
				popViewport()
				dev.off()
				} ## END loop over tmp_content_areas
			} ### END if ("PNG" %in% sgPlot.output.format)
			} ## END for loop for STUDENTS (n)
		} ## END for loop for GRADES (l)

		cat("\\end{document}", file=paste("instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), append=TRUE)
		system(paste("pdflatex -interaction=batchmode instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), ignore.stdout = TRUE)
		system(paste("pdflatex -interaction=batchmode instructor_catalog_", i, "_", j, "_", k, ".tex", sep=""), ignore.stdout = TRUE)
		file.rename(paste("instructor_catalog_", i, "_", j, "_", k, ".pdf", sep=""), file.path(sgPlot.folder, year_folder, district_folder, school_folder, instructor_folder,
			paste(year_folder, "_", district_folder, "_", school_folder, "_Individual_SGP_Report_Catalog.pdf", sep="")))
		if (sgPlot.cleanup & "PDF" %in% sgPlot.output.format) {
			files.to.remove <- list.files(pattern=paste(i, j, k, sep="_"), all.files=TRUE)
			lapply(files.to.remove, file.remove)
			files.to.remove <- list.files(path=file.path(sgPlot.folder, year_folder, district_folder, school_folder, instructor_folder), pattern="REPORT", all.files=TRUE, full.names=TRUE, recursive=TRUE)
			lapply(files.to.remove, file.remove)
		}
	
		if (identical(.Platform$OS.type, "unix") & sgPlot.zip) {
			if (1000*as.numeric(unlist(strsplit(system(paste("du -s", file.path(sgPlot.folder, year_folder, district_folder, school_folder, instructor_folder)), intern=TRUE), "\t"))[1]) > 4000000000) {
				tmp.working.directory <- getwd()
				setwd(file.path(sgPlot.folder, year_folder, district_folder, school_folder))
				if (paste(instructor_folder, ".tar.gz", sep="") %in% list.files()) file.remove(paste(instructor_folder, ".tar.gz", sep=""))
				system(paste("tar cfz", paste(instructor_folder, ".tar.gz", sep=""), instructor_folder, sep=" "))
				setwd(tmp.working.directory)
			} else {
				tmp.working.directory <- getwd()
				setwd(file.path(sgPlot.folder, year_folder, district_folder, school_folder))
				if (paste(instructor_folder, ".zip", sep="") %in% list.files()) file.remove(paste(instructor_folder, ".zip", sep=""))
				suppressMessages(
					zip(paste(instructor_folder, ".zip", sep=""), instructor_folder, flags="-rmq")
				)
				setwd(tmp.working.directory)
			}
			unlink(file.path(sgPlot.folder, year_folder, district_folder, school_folder, instructor_folder), recursive=TRUE)
		}

		message(paste("\tStarted", last.year, tmp_school_name, "Instructor:", tmp_instructor_name, "student growth plots:", started.date))
		message(paste("\tFinished", last.year, tmp_school_name, "Instructor:", tmp_instructor_name, "student growth plots:", date(), "in", convertTime(timetaken(started.at)), "\n"))

	} ## END for loop for INSTRUCTORS (k)
	} ## END for loop for SCHOOLS (j)
	} ## END for loop for DISTRICTS (i)

	return("DONE")

} ### END if (reports.by.instructor)

} ## END studentGrowthPlot_Styles function
