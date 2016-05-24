



get.pimgs <- function(rootdir = NULL, f2load = "p_ImageList.Rdata") {
	## basic collation of p_ImageLists
	##   - we pick up all the p_ImageLists, and then make a plot


	if (is.null(rootdir)) stop("rootdir must be set, with full path")
	## we set that as our working directory
	setwd(rootdir)  ## just like setting rootdir in /File/Change dir/

	## list every file in there
	fs <- list.files()

	## keep only the directories, we will visit each of these
	fs <- fs[file.info(fs)$isdir]

	## empty list to collect all the pImages
	l.p <- list()

	## name of the file we expect in each directory (could be "pFat_ImageList.Rdata", for that case)
	f2load <- "p_ImageList.Rdata"

	## now loop over each directory, check for the file we want, load and collect it if it is there
	for (f in fs)  {
	   ## make sure we are in the right place
	    setwd(rootdir)
	   setwd(f)
	   ## only load if the file we want is actually here
	   if (file.exists(f2load)) {
	       load(f2load)
	       p <- get("p")
	       l.p[[f]] <- p  ## add a named element to our list, containing this pImage (same name as the directory)
   		}
	}
	setwd(rootdir)  ## reset
	l.p
}
