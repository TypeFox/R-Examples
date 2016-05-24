fasta.read <-
function(file, seqlength = 648, pos1 = 1, pos2 = 3){
	fastaFile <- readLines(file)

	lines <- matrix(NA, nrow = length(fastaFile), ncol = 2)
	n <- 0
	for(i in 1:length(fastaFile)){
		if(substr(fastaFile[i], 1, 1) == ">"){
			n <- n + 1
			lines[n,1] <- substr(fastaFile[i],2,nchar(fastaFile[i]))
		}
		else{
			if(is.na(lines[n,2])){
				lines[n,2] = fastaFile[i]
			}
			else{
				lines[n,2] <- paste(lines[n,2], fastaFile[i], sep = "")
			}
		}
	}
	lines <- lines[-(n+1:length(fastaFile)),]

	sequence <- matrix(NA, nrow = nrow(lines), ncol = 2+seqlength)
	for(r in 1:nrow(lines)){
		sequence[r, 1:2] <- strsplit(lines[r,1], split = "\\|")[[1]][c(pos1,pos2)]
		sequence[r, 3:(2+seqlength)] <- strsplit(lines[r,2], split = "")[[1]]
	}
	return(sequence)
}
