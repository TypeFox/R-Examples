####
# Get the desired distribution of iRefIndex from the web repository:
####
get_irefindex = function(tax_id="All", iref_version="current", data_folder=getwd()) {

	# 1. Get data folder:
	if (data_folder == "data") {
		datafolder = system.file("data", package = "iRefR")
	} else if (data_folder == "home") {
		datafolder = R.home()
	} else {
		datafolder = data_folder
	}

	# 2. Get release dates and full URLs:
	if (iref_version == "current") {
		iref_version = "13.0"
		release_date = "08122013"
		url = paste("http://irefindex.org/download/irefindex/data/current/psi_mitab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
	} else {
		if (iref_version == "7.0") {
			release_date = "11042010"
			url = paste("http://irefindex.org/download/irefindex/data/archive/release_", iref_version, "/psimi_tab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		}
		if (iref_version == "8.0") {
			release_date = "01192011"
			url = paste("http://irefindex.org/download/irefindex/data/archive/release_", iref_version, "/psimi_tab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		}
		if (iref_version == "9.0") {
			release_date = "10182011"
			url = paste("http://irefindex.org/download/irefindex/data/archive/release_", iref_version, "/psimi_tab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		}
		if (iref_version == "10.0") {
			release_date = "03022013"
			url = paste("http://irefindex.org/download/irefindex/data/archive/release_", iref_version, "/psi_mitab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		}
		if (iref_version == "11.0") {
			release_date = "05192013"
			url = paste("http://irefindex.org/download/irefindex/data/archive/release_", iref_version, "/psi_mitab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		}
		if (iref_version == "12.0") {
			release_date = "06062013"
			url = paste("http://irefindex.org/download/irefindex/data/archive/release_", iref_version, "/psi_mitab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		}
		if (iref_version == "13.0") {
			release_date = "08122013"
			url = paste("http://irefindex.org/download/irefindex/data/archive/release_", iref_version, "/psi_mitab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		}
	}

	# 3. Check if file already exists. Otherwise, download and save:
	file_location = paste(datafolder, "/", tax_id,".mitab.",release_date,".txt",sep="")
	if (file.exists(file_location) == TRUE) {
		cat("Reading available iRefIndex file...\n")
		irefindex_tab = unique(read.table(file_location, header=TRUE, comment.char="", sep='\t', quote=""))
	} else {
		cat("Downloading iRefIndex file...\n")
		zipfile = paste(datafolder, "/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		download.file(url, destfile=zipfile)
		unzip(zipfile, exdir=datafolder)
		file.remove(zipfile)
		cat("Reading downloaded file...\n")
		txtfile = paste(datafolder, "/", tax_id,".mitab.",release_date,".txt", sep="")
		irefindex_tab = unique(read.table(txtfile, header=TRUE, comment.char="", sep='\t', quote=""))
		cat("File has been saved as:\n")
		cat(paste(txtfile, "\n"))
		save(file = paste(datafolder, "/", tax_id,".mitab.",release_date,".RData",sep=""), list = "irefindex_tab")
	}

	irefindex_tab
}
