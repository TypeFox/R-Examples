fetchBMdata = function(){
    f1 = file('http://www.fbi.gov/about-us/lab/forensic-science-communications/fsc/july1999/dnaloci.txt', 'r')

    if(!isOpen(f1))
        stop("Can't read from URL: http://www.fbi.gov/about-us/lab/forensic-science-communications/fsc/july1999/dnaloci.txt. You need to be connected to the internet for this to work")

    Lines = readLines(f1)
    close(f1)

#     ## Data is on a single line, number 432, in this file
# 
#     Lines = Lines[432]
#     Lines = unlist(strsplit(Lines, "<br />"))
    ## 2014-03-12 The format of the data has changed again. I will search for the start of the data 'Table 1' and the end '3409'
    Lines = unlist(strsplit(Lines, "<br />"))
    i1 = grep("Table 1.", Lines)
    i2 = grep("^3409", Lines)
    Lines = Lines[i1:i2]
    Lines = gsub("[&][#]13;", "", Lines)
    Lines = gsub("<[/]*p>","",Lines)

    ## All rares are encoded 108.1
    Lines = gsub("&(gt|lt);[0-9]+","108.1",Lines)

    ## drop the empty lines

    Lines = Lines[nchar(Lines)>0]

    ## drop the first line

    Lines = Lines[-1]

    ## The lines that have length < 50 are the population names

    popLocations = which(nchar(Lines) < 50)
    popNames = Lines[popLocations]

    ## Locus names are the same in all six populations
    ## However, last 3 pops only have CODIS loci
    ## First tag on locus line is ID# so drop that

    Loci = unlist(strsplit(Lines[2], "[\t ]+"))[-1]
    nLoci = length(Loci)

    ## get the alleles possible at each locus

    Alleles = vector(length = nLoci, mode = "list")

    allPops = strsplit(Lines[-c(popLocations,popLocations+1)],"[\t ]+")
    allPopsM = matrix(0, ncol = 2*nLoci, nrow = length(allPops))

    i1 = which(sapply(allPops, length) == 2*nLoci + 1)

    for(i in i1){
        allPopsM[i,] = allPops[[i]][-1]
    }

    i1 = (1:nrow(allPopsM))[-i1]

    for(i in i1){
        allPopsM[i,1:26] = allPops[[i]][-1]
    }

    allPops = allPopsM

    for(loc in 1:nLoci){
        i1 = 2*(loc - 1)  + 1
        i2 = 2*loc

        tbl = table(allPops[,i1:i2])
        A = names(tbl)
        A = A[A!="0"]

        if(!any(grepl("[^0-9.]", A))){
            A = sort(as.numeric(A))
        }

        Alleles[[loc]] = A
    }

    names(Alleles) = Loci

    db = vector(length = 6, mode = "list")
    names(db) = popNames

    i1 = popLocations + 2
    i2 = c(popLocations[-1] - 1, length(Lines))

    for(pop in 1:6){
        if(pop > 3){
            nLoci = 13 ## the last three pops only have 13 loci
            Loci = Loci[1:nLoci]
        }

        db[[pop]] = list(loci = Loci,
                         freqs = vector(mode = "list", length = nLoci))

        names(db[[pop]]$freqs) = Loci
        #browser()
        
        for(loc in 1:nLoci){
            db[[pop]]$freqs[[loc]] = rep(0, length(Alleles[[loc]]))
            names(db[[pop]]$freqs[[loc]]) = Alleles[[loc]]
        }

        popLines = strsplit(Lines[i1[pop]:i2[pop]],"[\t ]+")
        popM = matrix(0, nrow = length(popLines), ncol = 2*nLoci)


        j1 = which(sapply(popLines, length) == 2*nLoci + 1)

        for(j in j1){
            popM[j,] = popLines[[j]][-1]
        }

        j1 = (1:nrow(popM))[-j1]

        for(i in j1){
            popM[j,1:26] = popLines[[j]][-1]
        }

        for(loc in 1:nLoci){
            j1 = 2*(loc - 1)  + 1
            j2 = 2*loc

            tbl = NULL

            ## 1st 14 loci have numeric alleles
            if(loc <=14){
                tbl = table(as.numeric(popM[,j1:j2]))
                A = as.numeric(names(tbl))
                tbl = tbl[A!=0]
                A = A[A!=0]

                tbl = tbl/sum(tbl)

                pos = match(A, Alleles[[loc]])

                db[[pop]]$freqs[[loc]][pos] = tbl
            }else{
                tbl = table(popM[,j1:j2])
                A = names(tbl)
                tbl = tbl[A!="0"]
                A = A[A!="0"]

                tbl = tbl/sum(tbl)

                pos = match(A, Alleles[[loc]])

                db[[pop]]$freqs[[loc]][pos] = tbl
            }
        }
    }

    return(db)
}
