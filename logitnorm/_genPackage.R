#------------ generating package tw.DEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")

.tmp.f <- function(){
	library(twMisc)
	source("R/logitnorm.R")
	
	twUtestF()
}

.tmp.f <- function(){
	pkg <- "logitnorm"
	
	library(inlinedocs)
	unlink( file.path("man","*.Rd") )	
	package.skeleton.dx(".")
	
	# generate the HTML  files
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --html ",pkg, sep="") )
	setwd(prevWd)
	
	# show in Browser
	htmlRoot <- file.path( system.file(package = pkg), "html" )
	html_viewer <- function(path) {
		browser <- getOption("browser")
		if(is.null(browser) && .Platform$OS.type == "windows")
			shell.exec(chartr("/", "\\", path))
		else browseURL(paste("file://", URLencode(path), sep=""))
	}
	html_viewer(file.path(htmlRoot,"00Index.html"))
	
	updateVersionAndDate()
}

#R CMD check --no-vignettes --no-latex --no-install logitnorm
#R CMD check --no-vignettes --no-latex --no-codoc logitnorm
#R CMD INSTALL --html logitnorm

.tmp.f <- function(){
	library(sos)
	fres1 <- findFn("logitnormal") 
} 
