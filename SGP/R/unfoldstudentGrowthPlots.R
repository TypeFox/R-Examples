`unfoldstudentGrowthPlots` <- 
function(directory="Visualizations"){

	tmp.directory <- file.path(directory, "studentGrowthPlots_FLAT")
	dir.create(tmp.directory, recursive=TRUE, showWarnings=FALSE)

	### Get list of zip files

	list.of.zip.files <- list.files(file.path(directory, "studentGrowthPlots"), recursive=TRUE, pattern="zip")
	list.of.zip.files <- paste(getwd(), file.path(directory, "studentGrowthPlots"), list.of.zip.files, sep="/")

	### unzip the files

	setwd(tmp.directory)

	for (i in list.of.zip.files) {
		system(paste("unzip -oq", i))
	}

	### Get new list of individual files

	list.of.directories <- list.files()
	list.of.pdf.files <- list.files(recursive=TRUE)
	list.of.pdf.files <- paste(getwd(), list.of.pdf.files, sep="/")
	list.of.studentGrowthPlots <- list.of.pdf.files[-grep("Catalog", list.of.pdf.files)]
	list.of.pdf.catalogs <- grep("Catalog", list.of.pdf.files, value=TRUE)

	### Create directory for School Catalogs

	dir.create("School_Catalogs", recursive=TRUE, showWarnings=FALSE)
	dir.create("studentGrowthPlots", recursive=TRUE, showWarnings=FALSE)

	for (i in list.of.pdf.catalogs) {
		file.rename(i, file.path("School_Catalogs", sapply(strsplit(i, "/"), tail, 1)))
	}

	for (i in list.of.studentGrowthPlots) {
		file.rename(i, file.path("studentGrowthPlots", sapply(strsplit(i, "/"), tail, 1)))
	}

	### Remove empty directories

	for (i in list.of.directories) {
		unlink(file.path(getwd(), i), recursive=TRUE)
	}
} ### END unfoldstudentGrowthPlots function
