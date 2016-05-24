#scratch code for development of inlinedocs handling tabs and indentation

tmp.f <- function(){
	library(debug)
	setwd("inlinedocs/pkg/inlinedocs/tests")
	source("../R/package.skeleton.dx.R")
	prefix <- "^[ \t]*###[ \t]"			#changed the pattern to handle tabs at the beginning
	#mtrace(extract.docs.chunk)
	res <- suppressWarnings(test.pkg("indent"))
}

tmp.f2 <- function(){
	line="\t,third ,\t##<< an argument with comma before and online comment"
	arg.pat <- paste("^[^=#]*?([\\w\\.]+)\\s*([=,].*|\\)\\s*)?",
		"<<\\s*(\\S.*?)\\s*$",
		sep="##") # paste avoids embedded trigger fooling the system
	grep(arg.pat,line,perl=TRUE)
	gsub(arg.pat,"\\3",line,perl=TRUE);	#comment
	gsub(arg.pat,"\\\\item\\{\\1\\}",line,perl=TRUE) #arg
	
	(arg.pat <- paste("^[^=#]*?([\\w\\.]+)", sep="##"))
	arg.pat <- "[^[^=,#]*"
	arg.pat <- "[^[^=,#]*?"
}
