merlin = function(x, markers=seq_len(x$nMark), model=TRUE, theta=NULL, options="", verbose=FALSE, generate.files=TRUE, cleanup=generate.files, logfile="") {
	
	clean = function(cleanup, verbose, files) if (cleanup) {unlink(files); if(verbose) cat("Files successfully removed\n")}
	
	if (x$nMark == 0) stop("No markers exist for this linkdat object.")
	if(model && is.null(x$model)) stop("No model is set for this object")
    x = removeMarkers(x, seq_len(x$nMark)[-markers])
	map = .getMap(x, na.action=1, verbose=F)
	mNames = map$MARKER
	
	extensions = c("ped", "dat", "map", "freq", if(model) "model")
	
	if (generate.files) {
		files = write.linkdat(x, prefix="_merlin", what=extensions, merlin=TRUE)
		if (verbose) {cat("Files successfully generated:\n");  print(files)}
	}	
	
	options = paste(options, "--markerNames --quiet ")
	
	if (nonz_theta <- any(theta > 0)) {
		if (length(markers) > 1) {
			clean(cleanup, verbose, files)
			stop("Nonzero 'theta' values are possible only with a single marker.")
		}
        if (any(theta > 0.5)) {
			clean(cleanup, verbose, files)
			stop("Recombination fractions cannot exceed 0.5.")
		}
		pos = as.numeric(map[1, 3]) - 50 * log(1 - 2 * theta) #Haldane's map: Converting rec.fractions to cM positions.
		options = paste(options, " --positions:", paste(pos, collapse=","), sep="")
	}
	
	program = if(identical(x$model$chrom, "X")) "minx" else "merlin"
    command = paste(program, " -p _merlin.ped -d _merlin.dat -m _merlin.map -f _merlin.freq ", 
					if(model) "--model _merlin.model --tabulate ", options, sep="")
					
	if (verbose) cat("\nExecuting the following command:\n", command, "\n\n", sep="")

	merlinout = suppressWarnings(system(command, intern=T))
	clean(cleanup, verbose, files)
	if (nzchar(logfile))		write(merlinout, logfile)
	if (any(substr(merlinout,1,11) == "FATAL ERROR")) {
        cat("\n====================================\n", paste(merlinout[-(2:10)], collapse="\n"), 
        '====================================\n\n'); return(invisible())}
	if(verbose) {cat("Merlin run completed\n"); print(merlinout)}
	if (!is.na(skipped <- which(substr(merlinout,3,9) == "SKIPPED")[1])) stop(paste(merlinout[c(skipped-1, skipped)], collapse="\n"))
	
	if(!model) return(merlinout)

	##Extract LOD scores
    res = read.table("merlin-parametric.tbl", sep="\t", header=T, colClasses=c('numeric','numeric','character','NULL','numeric','NULL','NULL')) # chrom, pos, marker names and LOD.
    if (cleanup) unlink("merlin-parametric.tbl")
    
    if(nonz_theta) {
		mlodsdim = c(length(theta), 1)
		dimnam = list(theta, mNames)
	}
	else {
        markernames = res$LABEL
        if(!all(markernames %in% mNames)) {
            markernames = paste(res$CHR, res$POS*100, sep="_") # create new map with markernames formed as "chr_pos". Merlin weirdness: POS is in Morgans??
            map = data.frame(CHR=res$CHR, MARKER = markernames, POS = res$POS*100)
        }
        mlodsdim = c(1, nrow(res))
        dimnam = list(0, markernames)
    }
    
    lodres = structure(res$LOD, dim=mlodsdim, dimnames=dimnam, analysis="mlink", map=map, class="linkres")
    return(lodres)
}
