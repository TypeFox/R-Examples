get.file.name <- function(dir, prefix="", key="", ending="", verbal=TRUE, onename=FALSE) {
# From given directory dir, returns all files that start from given prefix,
# and that end with given ending, and those that contain the keyword key. 
# param:
# onename - flag if TRUE, then will return only the LAST file name that fits the description, not all the files.
# Example:
#
# 

def.return <- character(0)

# Obtain all files in the directory
all.files <- dir(dir)

if(length(all.files) == 0) {
	if(verbal)
		print(paste("No files exist in directory: ", dir, sep=""))
	return(def.return)
}

# 1. Extract files with desired prefix
file.id <- grep(paste("^", prefix, sep=""), all.files)
if(length(file.id) == 0) {
	if(verbal)
		print(paste("No files starting with '", prefix, "' were found in directory: ", dir, sep=""))
	return(def.return)
}
all.files <- all.files[file.id]

# 2. Extract files with desired ending
file.id <- grep(paste(ending, "$", sep=""), all.files)
if (length(file.id) == 0) {

	if(verbal) {
		if(prefix != "")
			print(paste("No files ending with '", ending, "' (and starting with '", prefix, "') were found in directory: ", dir, sep=""))
		else
			print(paste("No files ending with '", ending, "' were found in directory: ", dir, sep=""))
	}

	return(def.return)
}
all.files <- all.files[file.id]

# 3. Extract files with desired key
file.id <- grep(key, all.files)
if (length(file.id) == 0) {

	if(verbal) {
		if(prefix != "" && ending != "") {
			print(paste("No files containg '", key, "' (and starting with '", prefix, "', and ending with '", ending, "') were found in directory: ", dir, sep=""))
		} else if (prefix != "") {
			print(paste("No files containg '", key, "' (and starting with '", prefix, "') were found in directory: ", dir, sep=""))

		} else if (ending != "") {
			print(paste("No files containg '", key, "' (and ending with '", ending, "') were found in directory: ", dir, sep=""))
		} else {
			print(paste("No files containg '", key, "' were found in directory: ", dir, sep=""))
		}
	}
        return(def.return)
}
all.files <- all.files[file.id]

if(onename) { # return just the last one.
	return(all.files[length(all.files)]) 
}

return(all.files)


}
