toNexus = function(Pop, fileName = 'output.nex'){
    if(class(Pop) != "population")
        stop("Pop must be an object of class 'population'")

    nLoci = Pop$nLoci
    ns = Pop$nSubpops

    f1 = file(fileName, 'w')

    if(!isOpen(f1)){
        stop(paste("Could not open", fileName, "for writing\n"))
    }

    writeLines(c("#nexus", "", "begin gdadata;"), f1)

    line = paste("\tdimensions nloci=", nLoci, "npops=", ns, ";")
    writeLines(line, f1)

    writeLines(c("\tformat tokens labels missing=? datapoint=standard;",
                 "\tlocusallelelabels"), f1);

    locusLines = paste('\t\t', 1:nLoci, ' ',
                       "'", paste('Locus_', Pop$Freqs$loci, sep = ''),
                       "',", sep = "")
    locusLines[nLoci] = substr(locusLines[nLoci], 0,
                               nchar(locusLines[nLoci]) - 1)

    writeLines(locusLines, f1)
    writeLines(c(';', "MATRIX"), f1)

    pop = matrix(Pop$profiles, ncol = 2 * nLoci, byrow = T)
    Ns = nrow(pop) / ns

    toProf = function(x){
        n = length(x) / 2;
        i1 = 2 * (1:n) - 1; i2 = i1 + 1;
        paste(paste(x[i1], '/', x[i2], sep = ""), collapse = " ")
    }

    for(s in 1:ns){
        popLabel = paste('Pop_', s, ':', sep = '')
        writeLines(popLabel, f1)

        i1 = (s - 1)*Ns + 1
        i2 = s*Ns

        writeLines(paste('\t\t', 1:Ns, apply(pop[i1:i2,], 1, toProf)), f1)

        if(s <= (ns - 1))
            writeLines(',', f1)
    }
    writeLines(c("\t;", "END;"), f1)
    close(f1)

}
