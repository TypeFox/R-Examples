XML4R2list <- function(file){

	# CHECK IF FILE EXISTS
	if(!file.exists(file)) stop(paste0("'", file, "' does not exist."))

	# READ IN FILE
	read_lines <- readLines(file)
	
	# IF FILE IS NOT IN XML FORMAT, RETURN RESULT OF READ LINES
	if(sum(grepl('(>)', read_lines)) == 0) return(read_lines)
	
	# COLLAPSE AT LINE BREAKS IN FILE
	read_lines <- paste(read_lines, collapse="\n")

	# ADD LINE BREAKS WHERE NOT IN ORIGINAL FILE
	read_lines <- gsub('(>)([[:print:]])', '\\1\n\\2', read_lines)
	read_lines <- gsub('([[:print:]])(</)', '\\1\n\\2', read_lines)

	# SPLIT AT LINE BREAKS
	lines <- strsplit(x=read_lines, split="\n")[[1]]
	
	# REMOVE TABS
	lines <- gsub("^[\t]*", "", lines)

	# REMOVE EMPTY LINES
	lines <- lines[lines != ""]

	# DEFAULT INITIAL OBJECT TYPE
	object_type <- 'vector'

	# CHECK IF A SINGLE NODE/TAG ENCLOSES ALL OTHER CONTENTS
	enclosing_tag_added <- FALSE
	if(sum(grepl('(<)', lines[1:2])) < 2){
		
		# ADD ENCLOSING TAG
		lines <- c('<enclose_all type=list >', lines)
		lines <- c(lines, '</enclose_all>')
		
		# ENCLOSING TAG ADDED
		enclosing_tag_added <- TRUE
	}
	
	# SET RETURN LIST
	read_xml_lines <- XML4R2listLines(lines)
	
	if(enclosing_tag_added){
		rlist <- read_xml_lines$rlist
	}else{
		rlist <- list(read_xml_lines$rlist)
		names(rlist)[1] = read_xml_lines$obj.name
	}
	
	rlist

}