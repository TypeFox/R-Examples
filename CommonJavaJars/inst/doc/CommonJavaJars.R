### R code from vignette source 'CommonJavaJars.Rnw'

###################################################
### code chunk number 1: CommonJavaJars.Rnw:31-36
###################################################

options(width=140)
options(prompt=" ", continue=" ")
options(digits=4)



###################################################
### code chunk number 2: CommonJavaJars.Rnw:47-50 (eval = FALSE)
###################################################
## 
## .jpackage("CommonJavaJars", jars=c("forms-1.2.0.jar", "iText-2.1.4.jar"))
## 


###################################################
### code chunk number 3: CommonJavaJars.Rnw:58-61 (eval = FALSE)
###################################################
## 
## loadJars(c("forms", "iText"))
## 


###################################################
### code chunk number 4: CommonJavaJars.Rnw:90-102 (eval = FALSE)
###################################################
## # The following few lines are based on the code of the rJava .jpackage function
## if (!is.null(sessionInfo()$otherPkgs$rJava$Version) && sessionInfo()$otherPkgs$rJava$Version < "0.8-3") {
## 	classes <- system.file("R28", package = "CommonJavaJars", lib.loc = NULL)
## 	if (nchar(classes)) {
## 		.jaddClassPath(classes)
## 		jars <- grep(".*\\.jar", list.files(classes, full.names = TRUE), TRUE, value = TRUE)
## 		if (length(jars)) { 
## 			.jaddClassPath(jars)
## 		}		
## 	}
## }
## # Otherwise load from rJava. 


