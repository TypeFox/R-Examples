helloWorld <- function
### function for which an example should be added
( ){
  cat("Hello World!\n")
  invisible(NULL)
  ### None (invisible NULL)
}

attr(helloWorld,"ex") <- function(){
	# all text including comments inside this block 
	# should go to the examples section of 
	# function helloWorld
	helloWorld()	# prints Hello World
}

.tmp.f <- function(
	### some non-exported interactive debugging code.
){
	# first source test.R from inlinedocs
	fname <- file.path("inst","testfiles","attrExamples.R")
	#library(debug); mtrace(extract.file.parse)
	#save.test.result(fname)
	#test.file(fname)
}

.result <- 
	list(helloWorld = list(definition = "helloWorld <- function\n### function for which an example should be added\n( ){\n  cat(\"Hello World!\\n\")\n  invisible(NULL)\n  ### None (invisible NULL)\n}",  
			description = "function for which an example should be added",  
			value = "None (invisible NULL)", examples = "\n# all text including comments inside this block \n# should go to the examples section of \n# function helloWorld\nhelloWorld()\t# prints Hello World\n",
               format="",title="helloWorld"))
