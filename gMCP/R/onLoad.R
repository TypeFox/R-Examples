.onLoad <- function(libname, pkgname) {
	if (!.jniInitialized) {
		.jinit(parameters=c("-Xrs", "-Xss1m",
		                    paste0("-Djava.io.tmpdir=", tempdir())))
    # Remark - rJava 0.9-5: detect support for -Xrs and enable it by default (this prevents
    # Java from killing R process on interrupt)
	}
	.jpackage(pkgname)	
	.jpackage("JavaGD")
	
	jars <- c("afcommons", "commons-collections", "commons-lang", 
			"commons-logging", "commons-validator", "jgoodies-common", "forms", 
			"iText", "javax.json", "jhlir.jar", "jlatexmath", "jxlayer", 
			"log4j", "swing-worker")
	
	loadJars(jars)
	
	# Some of following lines are based on the code of the rJava .jpackage function
	if (length(grep("64-bit", sessionInfo()$platform))>0) { # Necessary?!
	  classes <- system.file("jri/x64", package = "rJava", lib.loc = NULL)
	} else {
	  classes <- system.file("jri/i386", package = "rJava", lib.loc = NULL)
	}
	if (nzchar(classes)) {
	  .jaddClassPath(classes)
	}
	classes <- system.file("jri", package = "rJava", lib.loc = NULL)
	if (nzchar(classes)) {
		#.jaddClassPath(classes)
		jars <- grep(".*\\.jar", list.files(classes, full.names = TRUE), TRUE, value = TRUE)
		if (length(jars)) { 
			.jaddClassPath(jars)
		}		
	}
  # Adding poi jar files to Classpath:
	classes <- system.file("java", package = "xlsxjars", lib.loc = NULL)
	if (nzchar(classes)) {	  
	  jars <- grep(".*\\.jar", list.files(classes, full.names = TRUE), TRUE, value = TRUE)
	  if (length(jars)) { 
	    .jaddClassPath(jars)
	  }		
	}
	
	rJavaVersion <- utils::sessionInfo()$otherPkgs$rJava$Version
	
	# If we have a rJava version < 0.8-3 load JRIEngine.jar and REngine.jar
  if (!is.null(rJavaVersion)) {
  	if (rJavaVersion < "0.8-3") {
			classes <- system.file("R28", package = "CommonJavaJars", lib.loc = NULL)
			if (nzchar(classes)) {
				.jaddClassPath(classes) # Necessary?!
				jars <- grep(".*\\.jar", list.files(classes, full.names = TRUE), TRUE, value = TRUE)
				if (length(jars)) { 
					.jaddClassPath(jars)
				}		
			}
		}	
		# If we have a rJava version 0.9-4 load JRIEngine.jar and REngine.jar    
		if (rJavaVersion == "0.9-4") {
		  classes <- system.file("JRI", package = "CommonJavaJars", lib.loc = NULL)
		  if (nzchar(classes)) {
		    .jaddClassPath(classes) # Necessary?!
		    jars <- grep(".*\\.jar", list.files(classes, full.names = TRUE), TRUE, value = TRUE)
		    if (length(jars)) { 
		      .jaddClassPath(jars)
		    }		
		  }
		}
  }
	
	# TODO Should we add "&& interactive()" or move it to some other place?
	if (!"jri.jar" %in% tolower(sapply(.jclassPath(), function(x) {substring(x, first=nchar(x)-6)}))) {
	  warning(paste(c("JRI.jar seems to be missing from the classpath.",
	                  "The graphical user interface will most likely not be available.",
	                  "Compile R with shared library enabled (--enable-R-shlib option)",
	                  "and reinstall rJava to use JRI functionality."), sep="\n"))
	}
  
  java.info <- getJavaInfo(FALSE, FALSE, TRUE)
	if (length(grep("-Xss1m", java.info))==0) {
	  warning(paste("JVM was already initialized with unknown memory settings:",strsplit(java.info, split="Input Arguments:")[1][[1]][2]))
	}
	
	## We supply our own JavaGD class
	Sys.setenv("JAVAGD_CLASS_NAME"="org/mutoss/gui/JavaGD")  
	
	# Optional Deducer integration:
	#deducer.addMenuItem <- if (exists("deducer")) get("deducer") else function(...) {}
	#jgr.addMenuSeparator <- if (exists("jgr.addMenuSeparator")) get("jgr.addMenuSeparator") else function(...) {}
	#jgr.addMenuItem <- if (exists("jgr.addMenuItem")) get("jgr.addMenuItem") else function(...) {}
	#deducer.addMenuItem("Multiple Test Graph", NULL, "graphGUI()", "Analysis")
	#jgr.addMenuSeparator("Analysis")
	#jgr.addMenuItem("Analysis", "Multiple Test Graph", "graphGUI()")	
}  

.onUnload <- function(libpath) {
	# TODO Unload jars?
}
