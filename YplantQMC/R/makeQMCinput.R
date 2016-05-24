
makeQMCinput <- function(plant, inputfile=NULL, outputfile=NULL, writefile=TRUE){
  
	if(is.character(plant$pfile)){
		proot <- gsub("\\.p","",plant$pfile,ignore.case=TRUE)
    proot <- basename(proot)
	} else {
		proot <- paste0("Plant",plant$nleaves)
	}
	if(is.null(inputfile))inputfile <- paste0(proot,".dat")
	if(is.null(outputfile))outputfile <- paste0(proot,".out")

	unlink(inputfile)
	unlink(outputfile)
	
	# Make sure there are no spaces in the file names (replace with '_').
	inputfile <- gsub(" ","_",inputfile)
	outputfile <- gsub(" ","_",outputfile)
	
	outr <- c()
  for(i in 1:plant$nleaves){  
    outr <- c(outr, paste(i, "E(0)"))
    outr <- c(outr, "Dbegin")
    outr <- c(outr, "polygon")
    
    xyz <- plant$leaves[[i]]$XYZ
    		
    # check if last point is the same as the first,
    # and then delete it (QuasiMC does not like it).
    if(sqrt(sum(xyz[1,] - xyz[nrow(xyz),])^2) < 1E-09)
    	xyz <- xyz[-nrow(xyz),]
    		
    # Note: QuasiMC assumes order X,Z,Y.
    xyz <- xyz[,c(1,3,2)]
  	
    outr <- c(outr, apply(xyz, 1, paste, collapse=" "))
    outr <- c(outr, "Dend")
  }
  outr <- c(outr, "Control: 8")
  outr <- c(outr, "")  # not sure if needed...

	if(!writefile){
		unlink(inputfile)
		return(list(qmcinput=outr, inputfile=inputfile, outputfile=outputfile))
	} else {
    writeLines(outr, inputfile)
		return(c(inputfile,outputfile))
	}
}

