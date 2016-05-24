readFreqs = function(strPath, FSIGenFormat = TRUE, delim = ','){
    Freqs = list(loci = NULL, freqs = NULL, counts = NULL)

    if(FSIGenFormat){
        ## expect the information to start with Alleles, and end with n
        f1 = file(strPath, 'r')

        if(isOpen(f1)){
            Lines = readLines(f1)
            close(f1)

            alleleStart = -1

            ## look for a line starting with alleles

            if(!any(grepl('^[Aa]llele', Lines))){
                stop("The header line should start with 'Allele'")
            }else{
                alleleStart = grep('^[Aa]llele', Lines)
                if(length(alleleStart) > 1){
                    stop("There are two or more lines starting with Allele")
                }
            }


            ## look for line starting with n
            countStart = -1

            if(!any(grepl('^[Nn]', Lines))){
                stop("The frequency file should have locus allele counts")
            }else{
                countStart = grep('^[Nn]', Lines)
                if(length(countStart) > 1){
                    warning("There should only one locus allele count line. Using first ")
                    countStart = countStart[1]
                }
            }

            if(countStart == -1 || alleleStart == -1){
                stop("Header and footer lines not found")
            }else{
                n = length(Lines)
                Lines = Lines[alleleStart:countStart]

                if((n - length(Lines)) > 0){
                    cat(paste("Dropped", n - length(Lines), "lines\n"))
                }

                locusLine = Lines[1]
                Tokens = unlist(strsplit(locusLine, delim))
                Loci = toupper(Tokens[-1])
                nLoci = length(Loci)
                cat(paste("Detected", nLoci, "loci\n"))

                countLine = Lines[length(Lines)]
                Tokens = unlist(strsplit(countLine, delim))
                counts = as.numeric(Tokens[-1])

                Lines = Lines[-c(1,length(Lines))]

                freqs = vector(length = nLoci, mode = "list")
                names(freqs) = Loci

                for(line in Lines){
                    Tokens = unlist(strsplit(line, delim))

                    a = Tokens[1]
                    fx = as.numeric(Tokens[-1])


                    for(loc in 1:nLoci){
                        if(!is.na(fx[loc]) && fx[loc] > 0){
                            f = fx[loc]

                            if(f > 1){
                                stop("Error: frequency greater than 1")
                            }

                            freqs[[loc]] = c(freqs[[loc]], f)
                            names(freqs[[loc]])[length(freqs[[loc]])] = a
                        }
                    }
                }

                cat(paste("Successfully processed", length(Lines), "lines\n"))

                Freqs$freqs = freqs
                Freqs$loci = Loci
                Freqs$counts = counts
            }
        }else{
            stop(paste("Couldn't open", strPath, "for reading"))
        }
    }else{ ## Curran format
        f1 = file(strPath, "r")

        if(isOpen(f1)){
            Lines = readLines(f1)
            nLines = length(Lines)
            ## cat(paste("Read", nLines, "lines from", strPath,"\n"))
            close(f1)
        }else{
            stop(paste("Couldn't open :", strPath, "for reading"))
        }

        Freqs = NULL

        nLoci = as.numeric(Lines[1])
        ## cat(paste(nLoci, "\n"))

        Loci = rep("", nLoci)
        freqs = vector(nLoci, mode = "list")
        currLine = 2

        for(nLoc in 1:nLoci){
            Tokens = unlist(strsplit(Lines[currLine], ","))
            currLine = currLine + 1

            Loci[nLoc] = Tokens[1]
            nAlleles = as.numeric(Tokens[2])
            Alleles = rep("", nAlleles)

            freqs[[nLoc]] = rep(0, nAlleles)

            for(nA in 1:nAlleles){
                Tokens = unlist(strsplit(Lines[currLine], ","))
                currLine = currLine + 1
                Alleles[nA] = Tokens[2]
                freqs[[nLoc]][nA] = as.numeric(Tokens[3])
            }

            names(freqs[[nLoc]]) = Alleles
        }

        names(freqs) = Loci
        Freqs$loci = Loci
        Freqs$freqs = freqs
    }
    return(Freqs)
}
