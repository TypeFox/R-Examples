readXML4R <- function(file){

	# CHECK IF FILE EXISTS
	if(!file.exists(file)) stop(paste0("'", file, "' does not exist."))

	# READ IN FILE
	read_lines <- readLines(file)
	
	# COLLAPSE AT LINE BREAKS IN FILE
	read_lines <- paste(read_lines, collapse="\n")

	# ADD LINE BREAKS WHERE NOT IN ORIGINAL FILE
	read_lines <- gsub('(>)([[:print:]])', '\\1\n\\2', read_lines)
	read_lines <- gsub('([[:print:]])(</)', '\\1\n\\2', read_lines)

	# SPLIT AT LINE BREAKS
	lines <- strsplit(x=read_lines, split="\n")[[1]]
	
	lines <- gsub("^[\t]*", "", lines)

	# DEFAULT INITIAL OBJECT TYPE
	object_type <- 'vector'

	# SET RETURN LIST
	read_xml_lines <- readXMLLines(lines)
	
	rlist <- list(read_xml_lines$rlist)

	names(rlist)[1] = read_xml_lines$obj.name
	
	rlist
}