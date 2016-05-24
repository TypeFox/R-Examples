revComp <- 
function (dna = NULL) 
{
    if(is.null(dna))
	{stop("Which sequence do you want to convert?")}
	dna = toupper(dna)
    stri = gsub("V1", "B", gsub("B", "V", gsub("V", "V1", gsub("D1", 
        "H", gsub("H", "D", gsub("D", "D1", gsub("M1", "K", gsub("K", 
            "M", gsub("M", "M1", gsub("R1", "Y", gsub("Y", "R", 
                gsub("R", "R1", gsub("G1", "C", gsub("C", "G", 
                  gsub("G", "G1", gsub("A1", "T", gsub("T", "A", 
                    gsub("A", "A1", dna))))))))))))))))))
    bb = strsplit(stri, "")[[1]]
    rev.stri = bb[order(seq_along(bb), decreasing = TRUE)]
    DNA = paste(rev.stri, sep = "", collapse = "")
    return(DNA)
}
