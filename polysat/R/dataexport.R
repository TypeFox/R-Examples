write.Structure <- function(object, ploidy, file="",
                            samples=Samples(object), loci=Loci(object),
                            writepopinfo=TRUE, extracols=NULL,
                            missingout=-9){
    if(!all(!is.na(Ploidies(object, samples=samples, loci=loci))))
        stop("Ploidies needed.")

    if(missing(ploidy))
    stop("Ploidy of file must be specified seperately from ploidies of dataset.")

    structdata <- data.frame(rowlabel=c("missing",rep(samples, each=ploidy)))
    thiscol = 1 # which column we are working with

    # Fill in PopInfo if that option is selected
    if(writepopinfo){
        if(!all(!is.na(PopInfo(object)[samples])))
            stop("PopInfo needed, or set writepopinfo to FALSE")
        thiscol <- thiscol+1
        structdata[[thiscol]]<-c(0,rep(PopInfo(object)[samples], each=ploidy))
        names(structdata)[thiscol] <- "PopInfo"
    }

    # Fill in extra columns.  Ignored if extracols=NULL.
    for(excol in dimnames(extracols)[[2]]){
        thiscol <- thiscol+1
        structdata[[thiscol]]<-c(0,rep(extracols[samples, excol], each=ploidy))
        names(structdata)[thiscol]<-excol
    }

    # Fill in the genotypes
    for(L in loci){
        thiscol<-thiscol+1
        alleles<-missingout #set up the vector of alleles
        for(s in samples){
            # Is genotype missing?
            if(isMissing(object,s,L)){
                thesealleles <- rep(missingout, ploidy)
            } else {
            #If missing data must be inserted to reflect a lower ploidy level:
            if(Ploidies(object, s, L) < ploidy){
                if(length(Genotype(object,s,L)) == Ploidies(object,s,L)){
                    #fully heterozygous genotype
                    thesealleles <- c(Genotype(object,s,L),
                                    rep(missingout,
                                        times=ploidy-Ploidies(object,s,L)))
                } else {
                  if(length(Genotype(object,s,L)) < Ploidies(object,s,L)){
                      #duplicate the first allele to get to the right ploidy
                      thesealleles <- c(Genotype(object,s,L),
                                        rep(Genotype(object,s,L)[1],
                                            times=Ploidies(object,s,L)-
                                            length(Genotype(object,s,L))),
                                        rep(missingout,
                                            times=ploidy-Ploidies(object,s,L)))
                  } else {
                      #randomly choose alleles to use if there are too many
                      thesealleles <- c(sample(Genotype(object,s,L)
                                       ,Ploidies(object,s,L),replace=FALSE),
                                        rep(missingout,
                                            times=ploidy-Ploidies(object,s,L)))
                      cat(paste("Randomly removed alleles:",s,L), sep="\n")
                  }
              }
            #If the individual has equal or greater ploidy to that
                #being used in the file:
            } else {
                if(length(Genotype(object,s,L)) == ploidy){
                    #genotype fills available spaces
                    thesealleles <- Genotype(object,s,L)
                } else {
                    if(length(Genotype(object,s,L)) < ploidy){
                        #duplicate the first allele to get to the right ploidy
                        thesealleles<-c(Genotype(object,s,L),
                                        rep(Genotype(object,s,L)[1],
                                            times= ploidy -
                                            length(Genotype(object,s,L))))
                    } else {
                        #randomly choose alleles to use if there are too many
                        thesealleles<-sample(Genotype(object,s,L),ploidy,
                                             replace=FALSE)
                        cat(paste("Randomly removed alleles:",s,L), sep="\n")
                    }
                }
            }}
            alleles<-c(alleles,thesealleles)
        }
        structdata[[thiscol]]<-alleles
        names(structdata)[thiscol]<-L
    }
   write.table(structdata, file=file, sep="\t", row.names=FALSE,
               quote=FALSE, col.names=TRUE)
}

write.GenoDive<-function(object,
                         digits=2, file="", samples=Samples(object),
                         loci=Loci(object)){
    if(!all(!is.na(Ploidies(object, samples=samples, loci=loci))))
        stop("Ploidies needed.
              Use estimatePloidy to get max number of alleles.")
    if(!all(!is.na(PopInfo(object)[samples])))
            stop("PopInfo needed.")
    if(!all(!is.na(Usatnts(object)[loci])))
        stop("Usatnts needed.")
    # get some info for the second line of the file
    numsam<-length(samples)
    numloc<-length(loci)
    numpop<-length(unique(PopInfo(object)[samples]))
    maxploidy<-max(Ploidies(object,samples,loci))

    # start a character vector to contain all the lines of the file
    lines<-c(Description(object),
             paste(numsam,numpop,numloc,maxploidy,digits, sep="\t"))
    lines[3:(numpop+2)]<-
        PopNames(object)[sort(unique(PopInfo(object)[samples]))]
    lines[numpop+3]<-paste("Population\tIndividual",
                           paste(loci,sep="",collapse="\t"),sep="\t")

    # enter population and sample names
    for(s in 1:numsam){
        lines[numpop+3+s]<-paste(PopInfo(object)[samples[s]],
                                 samples[s], sep="\t")
    }

    # process alleles and write them to lines
    for(L in loci){
        # replace missing data with zeros
        repmiss<-function(value){
            if(value[1]==Missing(object)){ 0 } else {value}
        }
        convertedalleles<-mapply(repmiss, Genotypes(object,samples,L),
                                 SIMPLIFY=FALSE)
        # convert alleles to repeat number
        divallele<-function(value){floor(value/Usatnts(object)[L])}
        convertedalleles<-mapply(divallele,convertedalleles,SIMPLIFY=FALSE)
        # subtract if necessary to get them to the right number of digits
        suballele<-function(value){if(value[1]!=0){
            value-(10^(digits-1))} else{
            0}}
        while(max(mapply(max,convertedalleles)) >= 10^digits){
            convertedalleles<-mapply(suballele, convertedalleles,
                                     SIMPLIFY=FALSE)
        }
        # For each sample, concatenate into strings and add to lines
        for(s in 1:numsam){
            # Convert alleles to character and set up string for concatenation
            charalleles<-as.character(convertedalleles[[s]])
            allelestring<-""
            # For each allele:
            for(ca in charalleles){
                allele<-ca
                # Add zeros if necessary
                while(nchar(allele) < digits){
                    allele<-paste("0",allele,sep="")
                }
                # add this allele to the string
                allelestring<-paste(allelestring,allele,sep="")
            }
            # Add the allele string to the line for that sample.
            lines[numpop+3+s]<-paste(lines[numpop+3+s],allelestring,sep="\t")
        }
    }

    # write the file
    cat(lines, file=file, sep="\n")
}

write.SPAGeDi<-function(object,samples=Samples(object),
                        loci=Loci(object),
                        allelesep="/", digits=2, file="",
                        spatcoord=data.frame(X=rep(1,length(samples)),
                                             Y=rep(1,length(samples)),
                                             row.names=samples)
                        ){
    if(!all(!is.na(Ploidies(object,samples,loci))))
        stop("Ploidies needed.")
    if(!all(!is.na(PopInfo(object)[samples])))
            stop("PopInfo needed.")
    if(!all(!is.na(Usatnts(object)[loci])))
        stop("Usatnts needed.")
    # name the rows of the spatial coordinates frame if not already done
    if(identical(row.names(spatcoord), as.character(1:dim(spatcoord)[1]))){
        row.names(spatcoord)<-samples}

    # set up data frame to contain genotypes
    gentable<-data.frame(Ind=samples,
                         Cat=PopNames(object)[PopInfo(object)[samples]],
                         spatcoord[samples,])

    # get genotype data by locus (column)
    for(L in loci){
        genotypesL<-Genotypes(object,samples,L)
        # replace missing data with zeros
        genotypesL[isMissing(object, samples, L)]<-0
        # divide and subtract to convert to repeats
        divallele<-function(value){floor(value/Usatnts(object)[L])}
        genotypesL<-mapply(divallele,genotypesL,SIMPLIFY=FALSE)
        suballele<-function(value){if(value[1]!=0){
            value-(10^(digits-1))} else{
            0}}
        while(max(mapply(max,genotypesL)) >= 10^digits){
            genotypesL<-mapply(suballele, genotypesL, SIMPLIFY=FALSE)
        }
        names(genotypesL) <- samples
        # add zeros up to correct ploidy, delete alleles if necessary
        zerostoadd<-Ploidies(object,samples,L) - mapply(length,genotypesL)
        names(zerostoadd)<-samples
        for(s in samples){
            if(length(genotypesL[[s]])==1){
                # replicate allele if totally homozygous
                genotypesL[[s]]<-rep(genotypesL[[s]],Ploidies(object,s,L))
            } else {
                if(zerostoadd[s] < 0){
                    # randomly remove alleles if there are too many
                    genotypesL[[s]]<-sample(genotypesL[[s]],
                                              Ploidies(object,s,L),
                                            replace=FALSE)
                    cat("Alleles randomly removed to get to ploidy:",L,s,"\n")
                } else {
                    # add zeros for partial heterozygotes
                    genotypesL[[s]]<-c(genotypesL[[s]],rep(0,zerostoadd[s]))
                }
            }
            # also make each allele the right number of digits if necessary
            if(allelesep==""){
                genotypesL[[s]]<-as.character(genotypesL[[s]])
                for(a in 1:length(genotypesL[[s]])){
                    while(nchar(genotypesL[[s]][a]) < digits){
                        genotypesL[[s]][a]<-paste(0,genotypesL[[s]][a],
                                                    sep="")
                    }
                }
            }
        }
        # concatenate into strings
        genvect<-mapply(paste,genotypesL,collapse=allelesep)
        # add the vector to the data frame
        gentable<-data.frame(gentable,genvect)
        names(gentable)[dim(gentable)[2]]<-L
    }
    # write file
    write.table(gentable,file="SpagTemp.txt",sep="\t",row.names=FALSE,
                col.names=TRUE,quote=FALSE)
    cat(paste(length(samples),length(unique(PopInfo(object)[samples])),
              dim(spatcoord)[2],
                length(loci),digits,max(Ploidies(object,samples,loci)),
              sep="\t"),
        "0",
        readLines("SpagTemp.txt"),"END",sep="\n",file=file)
}

write.GeneMapper<-function(object,file="",samples=Samples(object),
                           loci=Loci(object)){
    if(!all(!is.na(Ploidies(object,samples,loci))))
        stop("Ploidies needed to determine number of columns.")
    # figure out how many allele columns are needed
    numallelecol<-max(Ploidies(object,samples,loci))
    # figure out how many rows are needed
    numrows<-length(samples)*length(loci)
    # set up data frame
    gentable<-data.frame(Sample.Name=rep(samples, times=length(loci)),
                         Marker=rep(loci, each=length(samples)))
    # put empty allele columns into data frame
    alcollabels<-paste("Allele.",1:numallelecol,sep="")
    for(ac in 1:numallelecol){
        gentable[[ac+2]]<-rep("",times=numrows)
        names(gentable)[ac+2]<-alcollabels[ac]
    }

    # put alleles into their respective cells
    currentrow<-1
    for(L in loci){
        for(s in samples){
            thesealleles<-as.character(Genotype(object, s, L))
            while((length(thesealleles) + 2) > dim(gentable)[2])
                gentable <- data.frame(gentable, rep("", numrows),
                                       stringsAsFactors=FALSE)
            gentable[currentrow,3:(length(thesealleles)+2)]<-thesealleles
            currentrow<-currentrow+1
        }
    }

    # write the table to file
    write.table(gentable, file=file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

write.Tetrasat<-function(object, samples=Samples(object), loci=Loci(object),
                         file=""){
    if(!all(!is.na(PopInfo(object)[samples])))
            stop("PopInfo needed.")
    if(!all(!is.na(Usatnts(object)[loci])))
        stop("Usatnts needed.")
    if(!all(Ploidies(object,samples,loci) == 4))
        stop("Ploidy must be 4.")

    # a vector of populations to cycle through
    allpops<-unique(PopInfo(object)[samples])
    # set up vector of lines to write to file
    lines<-c(Description(object)[1],loci)
    # where does population/genotype data begin?
    datastart<-length(lines) + 1

    # make lines with "Pop" and sample names, and an index of where samples are
    sampleindex<-as.integer(c())
    currentline<-datastart
    for(pop in allpops){
        lines[currentline]<-"Pop"
        currentline<-currentline+1
        # make a vector of all samples that belong to this population
        thesesamples<-samples[PopInfo(object)[samples]==pop]
        # make an index of the lines that they will go on
        indextoadd<-currentline:(currentline+length(thesesamples)-1)
        names(indextoadd)<-thesesamples
        # add this to the sample index
        sampleindex<-c(sampleindex,indextoadd)
        # write each sample name (plus spaces up to 20 characters) to the
        # appropriate line
        for(s in thesesamples){
            samname<-s
            while(nchar(samname)<20){
                samname<-paste(samname," ", sep="")
            }
            lines[currentline]<-samname
            currentline<-currentline+1
        }
    }

    # now go through the loci, convert the genotypes, and fill them in
    for(L in loci){
        # convert alleles to repeat number
        divallele<-function(value){
            if(value[1]==Missing(object)){
                value
            } else {
                floor(value/Usatnts(object)[L])
            }
        }
        convertedalleles<-mapply(divallele, Genotypes(object,samples,L),
                                 SIMPLIFY=FALSE)
        names(convertedalleles) <- samples
        # make sure alleles have two digits or fewer
        suballele<-function(value){
            if(value[1]==Missing(object)){
                value
            } else {
                value-10
            }
        }
        while(max(mapply(max,convertedalleles)) >= 100){
            convertedalleles<-mapply(suballele, convertedalleles,
                                     SIMPLIFY=FALSE)
        }
        # go through the genotypes by sample
        for(s in samples){
            # convert missing data to blank spaces,
            # convert to alleles to characters
            if(convertedalleles[[s]][1]==Missing(object)){
                charalleles<-"         "
            } else {
                charalleles<-as.character(convertedalleles[[s]])
                # duplicate allele if fully homozygous
                if(length(charalleles)==1){
                    charalleles<-rep(charalleles,4)
                }
                # randomly remove alleles if there are more than 4
                if(length(charalleles)>4){
                    charalleles<-sample(charalleles,4,replace=FALSE)
                    cat(c("Alleles randomly removed:",s,L,"\n"),sep=" ")
                }
            }
            # concatenate all alleles into one string
            allelestring<-""
            for(ca in charalleles){
                allele<-ca
                # Add zeros if necessary
                while(nchar(allele) < 2){
                    allele<-paste("0",allele,sep="")
                }
                # add this allele to the string
                allelestring<-paste(allelestring,allele,sep="")
            }
            # add spaces up to nine characters
            while(nchar(allelestring) < 9){
                allelestring<-paste(allelestring," ", sep="")
            }
            # add the concatenated genotype to the appropriate line
            lines[sampleindex[s]]<-paste(lines[sampleindex[s]],allelestring,
                                         sep="")
        }
    }

    cat(lines, sep="\n", file=file)
}

write.ATetra<-function(object, samples=Samples(object), loci=Loci(object),
                       file=""){
    if(!all(!is.na(PopInfo(object)[samples])))
            stop("PopInfo needed.")
    if(!all(Ploidies(object,samples,loci) == 4))
        stop("Ploidy must be 4.")

    # set up a character vector to hold the lines for the file
    lines<-paste("TIT",Description(object)[1], sep=",")
    currentline<-2 # a variable to say which line to write to next
    # make numbers to go with loci, pops, and samples
    locnums<-1:length(loci)
    names(locnums)<-loci
    samnums<-1:length(samples)
    names(samnums)<-samples
    popnames <- PopNames(object)[sort(unique(PopInfo(object)[samples]))]
    popnums<-sort(unique(PopInfo(object)[samples]))
    names(popnums)<-popnames

    # fill in data using a loop
    for(L in loci){
        # write a line describing the locus
        lines[currentline]<-paste("LOC",locnums[L],L,sep=",")
        currentline<-currentline+1
        for(p in popnames){
            # write a line describing the population
            lines[currentline]<-paste("POP",locnums[L],popnums[p],p, sep=",")
            currentline<-currentline+1
            # get a vector of individuals in this population
            thesesamples<-Samples(object, populations=p)
            thesesamples<-thesesamples[thesesamples %in% samples]
            # write lines for individuals
            for(s in thesesamples){
                # first put sample info into the line
                lines[currentline]<-paste("IND",locnums[L],popnums[p],
                                          samnums[s],s,sep=",")
                # get the alleles and print a warning if there is missing data
                thesealleles<-Genotype(object, s, L)
                if(isMissing(object, s, L)){
                    thesealleles<-""
                    cat("Missing data:",s,L,"\n",sep=" ")
                }
                # take a random sample if there are more than 4
                if(length(thesealleles)>4){
                    thesealleles<-sample(thesealleles,4,replace=FALSE)
                    cat("More than 4 alleles:",s,L,"\n",sep=" ")
                }
                # add alleles to the line
                for(a in 1:4){
                    if(length(thesealleles)>=a ){
                      lines[currentline]<-paste(lines[currentline],
                                                thesealleles[a],sep=",")
                  } else {
                      lines[currentline]<-paste(lines[currentline],"",sep=",")
                  }
                }
                currentline<-currentline+1
            }
        }
    }

    # write the "END-record"
    lines[currentline]<-"END"

    # write the file
    cat(lines, sep="\n", file=file)
}

write.POPDIST <- function(object, samples=Samples(object),
                          loci=Loci(object), file=""){
    # error messages
    if(!all(!is.na(PopInfo(object)[samples])))
        stop("PopInfo needed")
    if(!all(!is.na(Ploidies(object,samples,loci))))
        stop("Ploidies needed")
    # start a vector of lines to write to the file
    Lines <- c(Description(object)[1], loci)
    # make a vector to index all the lines containing samples
    samindex <- integer(0)

    # add one population at a time
    for(p in unique(PopInfo(object)[samples])){
        Lines <- c(Lines, "Pop")
        # add each individual
        for(s in samples[samples %in% Samples(object, populations=p)]){
            Lines <- c(Lines, paste(PopNames(object)[p],"\t,",sep=""))
            samindex <- c(samindex, length(Lines))
            names(samindex)[length(samindex)] <- s
        }
        # warning message for mixed ploidy
        if(length(unique(
             Ploidies(object,samples[samples %in% Samples(object,
                                                           populations=p)],
                      loci)))>1)
            warning(paste("Mixed ploidy in population",p,
                          "; POPDIST may reject file."))
    }

    # Convert alleles to two digits then add to lines
    for(L in loci){
        thesegen <- Genotypes(object, samples, L)
        if(max(mapply(max, thesegen)) > 99){
            if(!is.na(Usatnts(object)[L]) && Usatnts(object)[L] > 1){
                thesegen <- array(mapply(function(value){
                        if(value[1]==Missing(object)){
                            value
                        } else {
                            floor(value/Usatnts(object)[L])
                        }},
                                   thesegen, SIMPLIFY=FALSE),
                                  dim=c(length(samples),1),
                                  dimnames=list(samples, L))
            }

            while(max(mapply(max, thesegen)) > 99){
                thesegen <- array(mapply(function(value){
                  if(value[1]==Missing(object)){
                      value
                  } else {
                      value - 10
                  }},
                                         thesegen, SIMPLIFY=FALSE),
                                  dim=c(length(samples),1),
                                  dimnames=list(samples, L))
            }
        }

        for(s in samples){
            if(thesegen[[s,L]][1]==Missing(object)){
                alleles <- "00"
            } else {
                alleles <- as.character(thesegen[[s,L]])
            }
            if(length(alleles) > Ploidies(object,s,L)){
                alleles <- sample(alleles, Ploidies(object,s,L))
                warning(paste("Alleles randomly removed:",s,L))
            }
            Lines[samindex[s]] <- paste(Lines[samindex[s]], " ", sep="")
            for(a in 1:Ploidies(object,s,L)){
                if(a > length(alleles)){
                    Lines[samindex[s]] <- paste(Lines[samindex[s]],"00",sep="")
                } else {
                    if(nchar(alleles[a])==1)
                        alleles[a] <- paste("0",alleles[a],sep="")
                    Lines[samindex[s]] <- paste(Lines[samindex[s]],
                                                alleles[a],sep="")
                }
            }
        }
    }

    cat(Lines, file=file, sep="\n")
}

gendata.to.genind <- function(object, samples=Samples(object),
                              loci=Loci(object)){
  # Errors
  if(!is(object, "gendata")) stop("genambig or genbinary object needed.")
  ploidy <- unique(Ploidies(object, samples, loci))
  if(length(ploidy)>1) stop("Single ploidy needed for genind.")
  if(is.na(ploidy)) stop("Please specify ploidy.")

  # Get genbinary if necessary
  if(is(object, "genambig"))
    object <- genambig.to.genbinary(object)

  # export genotypes table
  tab <- Genotypes(object, samples, loci)
  # convert missing data to NA
  tab[tab == Missing(object)] <- NA
  # make column headings usable
  dimnames(tab)[[2]] <- gsub(".","-", dimnames(tab)[[2]], fixed=TRUE)

  # create genind object
  x <- adegenet::genind(tab, pop=PopNames(object)[PopInfo(object)[samples]],
              ploidy=ploidy, type="PA")
  return(x)
}
