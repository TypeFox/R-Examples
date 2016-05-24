.onLoad <- function(libname, pkgname) {
	if (!.jniInitialized) {
	  .jinit(parameters=c("-Xrs", "-Xss1m",
	                      paste0("-Djava.io.tmpdir=", tempdir())))
	}
	.jpackage(pkgname)
	.jpackage("JavaGD")
	
	jars <- c("afcommons",
			"commons-logging", "commons-lang", "jgoodies-common", "forms",  
			"iText", "jhlir.jar", "jxlayer", 
			"log4j", "swing-worker")
	
	loadJars(jars)
	
  # The following few lines are based on the code of the rJava .jpackage function
  classes <- system.file("jri", package = "rJava", lib.loc = NULL)
	if (nchar(classes)) {
		.jaddClassPath(classes)
		jars <- grep(".*\\.jar", list.files(classes, full.names = TRUE), TRUE, value = TRUE)
		if (length(jars)) { 
			.jaddClassPath(jars)
		}		
	}
  
	## we supply our own JavaGD class
	Sys.setenv("JAVAGD_CLASS_NAME"="org/mutoss/gui/JavaGD")  
		
	rJavaVersion <- utils::sessionInfo()$otherPkgs$rJava$Version
	
	# If we have a rJava version 0.9-4 load JRIEngine.jar and REngine.jar    
	if (!is.null(rJavaVersion) && rJavaVersion == "0.9-4") {
	  classes <- system.file("JRI", package = "CommonJavaJars", lib.loc = NULL)
	  if (nzchar(classes)) {
	    .jaddClassPath(classes) # Necessary?!
	    jars <- grep(".*\\.jar", list.files(classes, full.names = TRUE), TRUE, value = TRUE)
	    if (length(jars)) { 
	      .jaddClassPath(jars)
	    }		
	  }
	}
	
	assign(".summary_table",  buildSummaryTable(), envir=Crossover.env)
}  

.onAttach <- function(libname, pkgname) {
    #packageStartupMessage("************************************************************\n",
    #                      "* THIS PACKAGE IS AN ALPHA VERSION!                        *\n",
    #                      "* It may give you good designs or it may make mistakes.    *\n",
    #                      "* Please check for every design whether it fits your needs.*\n",
    #                      "************************************************************")
}