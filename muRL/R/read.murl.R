read.murl <- function(file = "murljobs.csv", header = TRUE, stringsAsFactors = 								FALSE, field.title = "title", field.fname = "fname" , 
								field.lname = "lname", fields.address = "address", 
								field.city = "city", field.state = "state", 
								field.zip = "zipcode", field.position = "position", 
								field.subfield = "subfield",  
								field.dept = "dept", field.institution = "institution",
								colClasses = c("character"), ...){
	
	## read data from local file or webpage
	if(is.data.frame(file)){
		data <- file
	}else{
	data <- read.csv(file, header = header, stringsAsFactors = stringsAsFactors, colClasses = colClasses, ...)
	}

	## Rename field titles to those expected by write.murl
	names(data)[which(names(data) == field.title)] <- "title"
	names(data)[which(names(data) == field.fname)] <- "fname"
	names(data)[which(names(data) == field.lname)] <- "lname"
	names(data)[grep(fields.address, names(data))] <- paste("address", 1:length(grep(fields.address, names(data))), sep="")

	names(data)[which(names(data) == field.city)] <- "city"
	names(data)[which(names(data) == field.state)] <- "state"
	names(data)[which(names(data) == field.zip)] <- "zip"		
	names(data)[which(names(data) == field.position)] <- "position"
	names(data)[which(names(data) == field.subfield)] <- "subfield"
	names(data)[which(names(data) == field.dept)] <- "dept"
	names(data)[which(names(data) == field.institution)] <- "institution"
	
	## Which are not == 5, 10?
	n510 <- which(!(nchar(data$zip) %in% c(5,10)))
	## Which are == 10, but lack hyphen?
	zl10 <- which(nchar(data$zip) == 10)
	nhyp <- zl10[-grep("-", data$zip[zl10])]
	## Of == 10 with hyphen, are they 5 then 4?
	yhyp <- zl10[grep("-", data$zip[zl10])]
	spl54 <- function(ddd){
		nchar(ddd) == c(5,4)
	}
	if(length(zl10)>0){
			yhypn <- yhyp[which(unlist(lapply(lapply(strsplit(data$zip[zl10[grep("-", data$zip[zl10])]], "-"), spl54), sum)) !=2)]
	}else{yhypn <- NULL}
	
	allns <- sort(unique(c(n510, nhyp, yhypn)))
	
	if(length(allns)>0){
		warning("At least one postal code is of non-standard (US) form.  Please check postal code(s) of unit(s) ", paste(allns, collapse=" "), " and respecify if needed.")
	}
		return(data)	
}
