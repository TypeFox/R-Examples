read.BOLD <- function(IDs){
	allSeqs <- list()
	
	for(i in 1:length(IDs)){
		URL <- paste("http://v3.boldsystems.org/index.php/Public_RecordView?processid=", IDs[i], sep = "")
		res <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
		#DNA sequence
		dnaLineNo <- grep("Locus", res)
		dnaBlock <- paste(res[dnaLineNo: (dnaLineNo + 19)], collapse="")
		dnaSeq <- strsplit(dnaBlock, split = "<pre>|</pre>")[[1]][2]
		gene <- strsplit(strsplit(dnaBlock, split = "td")[[1]][4], split = "<|>")[[1]][2]
		
		#Taxonomy---Subfamily, genus, species, BIN number
		taxLineNo <- grep("TAXONOMY", res)
		taxBlock <- paste(res[taxLineNo: (taxLineNo + 30)], collapse="")
		taxFields <- strsplit(gsub("\\t", "", taxBlock), split = "<tr>|</tr>")
		sciname <- strsplit(taxFields[[1]][8], split = ">|<")[[1]][29]
		BIN <- strsplit(taxFields[[1]][10], split = ">|<")[[1]][31]
		
		seqs <- as.DNAbin(list(strsplit(tolower(dnaSeq), split="")[[1]]))
		nam <- paste(IDs[i], sciname, BIN, sep = "|")
		names(seqs) <- nam
		attr(seqs, "species") <- sciname
		attr(seqs, "BIN") <- BIN
		attr(seqs, "accession_num") <- IDs[i]
		attr(seqs, "gene") <- gene
		allSeqs[[i]] <- seqs
	}
	collSeqs <- do.call(c, allSeqs)
	attr(collSeqs, "species") <- unlist(lapply(allSeqs, function(x) attr(x, "species")))
	attr(collSeqs, "accession_num") <- unlist(lapply(allSeqs, function(x) attr(x, "accession_num")))
	attr(collSeqs, "BIN") <- unlist(lapply(allSeqs, function(x) attr(x, "BIN")))
	attr(collSeqs, "gene") <- unlist(lapply(allSeqs, function(x) attr(x, "gene")))
	collSeqs
}