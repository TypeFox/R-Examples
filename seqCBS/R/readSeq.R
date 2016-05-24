readSeq <-
function(filename, formatName) {
	if(formatName=="ELANDPaired") {
		return(readSeqELANDPaired(filename))
	}
	else if(formatName=="Chiang") {
		return(readSeqChiang(filename))
	}
	else {
		print("Input Format Not Supported")
	}
}

