SweaveAll <- function(SweaveFiles, make=1, PostSweaveHook=NULL, weave=utils::Sweave, ...) {
    i <- 0
    result <- character()
    while (i < length(SweaveFiles)) {
        i <- i+1
        suppressWarnings(remove(".SweaveFiles", ".TexRoot", ".PostSweaveHook", ".SweaveMake",
                                envir=globalenv()))
        thisfile <- weave(SweaveFiles[i], ...)
    	result <- c(result, thisfile)
    	.PostSweaveHook <- PostSweaveHook
    	if (exists(".PostSweaveHook", envir=globalenv())) 	
    	    .PostSweaveHook <- get(".PostSweaveHook", envir=globalenv())
    	if (exists(".SweaveMake", envir=globalenv()))
    	    make <- get(".SweaveMake", envir=globalenv())
    	if (!is.null(.PostSweaveHook)) {
    	    .PostSweaveHook <- match.fun(.PostSweaveHook)
    	    .PostSweaveHook(thisfile)
    	}    
    	if (make && exists(".SweaveFiles", envir=globalenv())) {
    	    newfiles <- setdiff(get(".SweaveFiles", globalenv()), SweaveFiles)
            if (length(newfiles)) {
            	if (make == 1) {
            	    tex <- paste(tools::file_path_sans_ext(newfiles), ".tex", sep="")
            	    SweaveFiles <- c(SweaveFiles, newfiles[!file_test("-f", tex) | file_test("-nt", newfiles, tex)])
            	} else 
            	    SweaveFiles <- c(SweaveFiles, newfiles)
            }
        }
        if (exists(".TexRoot", envir=globalenv())) {
            TexRoot <- get(".TexRoot", globalenv())
            result <- c(TexRoot, result[result != TexRoot])
        }
    }
    result
}
