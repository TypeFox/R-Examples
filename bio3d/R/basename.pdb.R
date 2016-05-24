basename.pdb <- function(x, mk4=FALSE) {
	##
	##- Extract PDB basename/identifier from filenames
	##   like "basename()" for PDB files
	##  E.g.:
	##       basename.pdb("/somedir/somewhere/1bg2_myfile.pdb")
	##  Will give: 1bg2_myfile
	##       basename.pdb("/somedir/somewhere/1bg2_myfile.pdb", TRUE)
        ##  Will give: 1bg2

  	y <- sub("\\.pdb$","", basename(x))
  	if(mk4) { y <- substr(y,1,4) } 
	names(y) <- x
	return(y)
}

