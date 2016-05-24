ConvertPeptide <- function(sequence, output = "elements", IAA = TRUE) {

	peptideVector <- strsplit(sequence, split = "")[[1]]

	if(output == "elements") {
	
		FindElement <- function(residue) {
		
			if(residue == "A") element <- c(C = 3, H = 5, N = 1, O = 1, S = 0)
			if(residue == "R") element <- c(C = 6, H = 12, N = 4, O = 1, S = 0)
			if(residue == "N") element <- c(C = 4, H = 6, N = 2, O = 2, S = 0)
			if(residue == "D") element <- c(C = 4, H = 5, N = 1, O = 3, S = 0)
			if(residue == "E") element <- c(C = 5, H = 7, N = 1, O = 3, S = 0)
			if(residue == "Q") element <- c(C = 5, H = 8, N = 2, O = 2, S = 0)
			if(residue == "G") element <- c(C = 2, H = 3, N = 1, O = 1, S = 0)
			if(residue == "H") element <- c(C = 6, H = 7, N = 3, O = 1, S = 0)
			if(residue == "I") element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
			if(residue == "L") element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
			if(residue == "K") element <- c(C = 6, H = 12, N = 2, O = 1, S = 0)
			if(residue == "M") element <- c(C = 5, H = 9, N = 1, O = 1, S = 1)
			if(residue == "F") element <- c(C = 9, H = 9, N = 1, O = 1, S = 0)
			if(residue == "P") element <- c(C = 5, H = 7, N = 1, O = 1, S = 0)
			if(residue == "S") element <- c(C = 3, H = 5, N = 1, O = 2, S = 0)
			if(residue == "T") element <- c(C = 4, H = 7, N = 1, O = 2, S = 0)
			if(residue == "W") element <- c(C = 11, H = 10, N = 2, O = 1, S = 0)
			if(residue == "Y") element <- c(C = 9, H = 9, N = 1, O = 2, S = 0)
			if(residue == "V") element <- c(C = 5, H = 9, N = 1, O = 1, S = 0)
			
			if(residue == "C" & IAA == FALSE) element <- c(C = 3, H = 5, N = 1, O = 1, S = 1)
			if(residue == "C" & IAA == TRUE) element <- c(C = 5, H = 8, N = 2, O = 2, S = 1)
			
			return(element)
			
		}
		
		resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
		for(i in 1:length(peptideVector)) { resultsVector <- FindElement(peptideVector[i]) + resultsVector }
    
    resultsVector <- resultsVector + c(C = 0, H = 2, N = 0, O = 1, S = 0)   # add water
		
		return(as.list(resultsVector))
	}

	if(output == "3letter") {
		
		FindCode <- function(residue) {
		
			if(residue == "A") let <- "Ala"
			if(residue == "R") let <- "Arg"
			if(residue == "N") let <- "Asn"
			if(residue == "D") let <- "Asp"
			if(residue == "C") let <- "Cys"
			if(residue == "E") let <- "Glu"
			if(residue == "Q") let <- "Gln"
			if(residue == "G") let <- "Gly"	
			if(residue == "H") let <- "His"
			if(residue == "I") let <- "Ile"
			if(residue == "L") let <- "Leu"
			if(residue == "K") let <- "Lys"
			if(residue == "M") let <- "Met"
			if(residue == "F") let <- "Phe"
			if(residue == "P") let <- "Pro"
			if(residue == "S") let <- "Ser"
			if(residue == "T") let <- "Thr"
			if(residue == "W") let <- "Trp"
			if(residue == "Y") let <- "Tyr"
			if(residue == "V") let <- "Val"
			
			return(let)
				
		} 

		codes <- sapply(peptideVector, FindCode)
		return(paste(codes, collapse = ""))
	
	}
	
}
