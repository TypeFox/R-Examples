read.seerstat = function(fileName, keep.missing = TRUE) {
	dicFile = paste(fileName, "dic", sep = ".");
	lines = readLines(dicFile);
	header = NULL;
	
	for (i in 1:length(lines)) {
		if (length(grep("Page Variables", lines[i], fixed = TRUE)) == 1) {
			k = 1;
			while(lines[i + k] != "") {
				words = unlist(strsplit(lines[i + k], "Name="));
				if (length(words) == 2) {
					header = c(header, words[2]);
				}
				
				k = k + 1;
			}
		}
	}
	print(header);
	header2 = gsub("[,()<>={}!@#$%^&*+-]", "", header);
	header3 = gsub(" ", "_", header2);
	
	txtFile = paste(fileName, "txt", sep = ".");
	M = scan(txtFile, na.strings = ".", sep = "\t");
	M1 = matrix(M, length(header3), length(M) / length(header3));
	M2 = t(M1);
	M2 = read.table(txtFile, header = FALSE, na.strings = ".");
	names(M2) = header3;
	
	if (!keep.missing) {
		tempNa = apply(M2, 1, mean);
		bNa = sapply(tempNa, is.na);
		M2 = M2[!bNa, ];
	}
	return(M2);
}

