ConvFas <-
function (fil = NULL, type = c("fas", "nxs", "phy")) 
{
    #if(is.null(fil))
	#{stop("Please specify the input data.")}
	match.arg(type)
	#dna = readLines(fil)
	dna = fil
    if (type == "fas") {
        seqNamPos = grep("^>", dna)
        pos = c(seqNamPos, length(dna) + 1)
        seqNam = dna[seqNamPos]
    }
    if (type == "nxs") {
         dna = dna[(grep("matrix", dna, ignore.case = TRUE) + 1):length(dna)]
         dna = dna[-which(dna == "" | dna == "end;" | dna == ";")]
         seqNam = unique(substr(dna, 1, regexpr(" ", dna) - 1))
        }
    if (type == "phy") {
         dna = dna[regexpr("[ATGC-]", dna) > 0]
         seqNam = substr(dna, 1, regexpr(" ", dna) - 1)
         seqNam = seqNam[-which(seqNam == "")]
        }
	   
	
    nSeq = length(seqNam)
    for (i in 1:nSeq) {
        if (type == "fas") {
            st = pos[i] + 1
            ed = pos[i + 1] - 1
            stri = gsub(" ", "", paste(dna[st:ed], collapse = ""))
        }
        if (type == "nxs" | type == "phy") {
            nBlock = length(dna)/length(seqNam)
            rNam = ((1:nBlock) - 1) * nSeq + i
            stri = gsub("[ -]", "", gsub(seqNam[i], "", paste(dna[rNam], 
                   collapse = "")))
		    stri = toupper(stri)
        }
        Nam = paste(">", seqNam[i], sep = "")
        if (i == 1){ 
            DNA = c(Nam, stri)
	     }
        if (i > 1){ 
            DNA = c(DNA, Nam, stri)
	     }          
	    }
    DNA = gsub(">>", ">", DNA)
	class(DNA) <- "fasta"
    return(DNA)
}
