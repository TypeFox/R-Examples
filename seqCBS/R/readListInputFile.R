readListInputFile <-
function(inputFilename, sep="\t") {
	inputFiles = read.table(inputFilename, sep=sep, header=TRUE, as.is=TRUE)
	colInput = colnames(inputFiles)
	if("File" %in% colInput && "Type" %in% colInput) {
		if(sum(inputFiles$Type=="Normal")==0 || sum(inputFiles$Type=="Tumor")==0) {
			print("Need Both Normal and Tumor Samples to work")
		}
		else {
			normalFiles = inputFiles$File[inputFiles$Type=="Normal"]
			tumorFiles = inputFiles$File[inputFiles$Type=="Tumor"]
			return(list(normalFiles=normalFiles, tumorFiles=tumorFiles))
		}
	}
	else {
		print("Input List incomplete or Column Name error")
	}
}

