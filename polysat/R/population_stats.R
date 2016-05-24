# Non-iterative function for calculating allele frequencies
# under polysomic inheritance.  For allele copy number ambiguity,
# assumes that all alleles have an equal chance of being present in more
# than one copy.
simpleFreq <- function(object, samples=Samples(object), loci=Loci(object)){
  # subset the object to avoid a lot of indexing later in the function
  object <- object[samples,loci]
  
    # we need PopInfo and Ploidies; check for these
    if(!all(!is.na(PopInfo(object)))){
        cat("PopInfo needed for samples:",
            samples[is.na(PopInfo(object))], sep="\n")
        stop("PopInfo needed for this function.")
    }
    if(!all(!is.na(Ploidies(object)))){
        stop("Ploidies needed for this function.")
    }

    # we need a genbinary object for this function; make one if necessary
    if(is(object, "genambig")){
        object <- genambig.to.genbinary(object)
    }
    Present(object) <- as.integer(1)
    Absent(object) <- as.integer(0)

    # get the populations that will be evaluated
    pops <- PopNames(object)[unique(PopInfo(object))]

    # Get ploidies into a format where you can total over the populations
    # (and don't index by locus if you don't need to)
    if(is(object@Ploidies, "ploidymatrix")){
      object <- reformatPloidies(object, output="collapse", erase=FALSE)
    }
    if(!is(object@Ploidies, "ploidysample") &&
       length(unique(Ploidies(object, samples, loci)))==1){
      object <- reformatPloidies(object, output="sample", erase=FALSE)
    }
    if(is(object@Ploidies, "ploidylocus")){
      object <- reformatPloidies(object, output="matrix", erase=FALSE)
    }

    # Get total genomes per population if loci are uniform ploidy
    if(is(object@Ploidies, "ploidysample")){
      # get the total number of genomes per population
      totgenomes <- rep(0, length(pops))
      names(totgenomes) <- pops
      for(p in pops){
          totgenomes[p] <- sum(Ploidies(object)[Samples(object,populations=p)])
      }

          # set up data frame to contain allele frequencies
      freqtable <- data.frame(Genomes=totgenomes,row.names=pops)
    } else { # or if loci have different ploidies
      freqtable <- data.frame(row.names=pops)
    }
    
    # loop to get frequency data
    for(L in loci){
        # get all samples without missing data at this locus
        xsamples <- samples[!isMissing(object, samples, L)]
        # get the total number of genomes per population
        totgenomes <- rep(0, length(pops))
        names(totgenomes) <- pops
        for(p in pops){
            totgenomes[p] <- sum(Ploidies(object, xsamples[xsamples %in%
                                                           Samples(object,
                                                           populations=p)], L))
        }
        # write genomes to the table if necessary
        if(is(object@Ploidies, "ploidymatrix")){
          totg <- data.frame(totgenomes)
          names(totg) <- paste(L, "Genomes", sep=".")
          freqtable <- cbind(freqtable, totg)
        }
        # make a conversion factor to weight allele presence based on ploidy
        # of each individual and number of alleles at this locus
        numalleles <- rep(0, length(xsamples))
        names(numalleles) <- xsamples
        for(s in xsamples){
            numalleles[s] <- sum(Genotype(object, s , L))
        }
        convf <- as.vector(Ploidies(object,xsamples,L))/numalleles
        # make table of weighted allele presence
        loctable <- Genotypes(object, xsamples, L) * convf
        # loop through all alleles at this locus
        for(al in dimnames(loctable)[[2]]){
            theseallelefreqs <- rep(0, length(pops))
            names(theseallelefreqs) <- pops
            # loop through populations
            for(p in pops){
                theseallelefreqs[p]<-sum(loctable[xsamples[xsamples %in%
                                                Samples(object, populations=p)],
                                                  al])/totgenomes[p]
            }
            freqtable <- cbind(freqtable,theseallelefreqs)
            names(freqtable)[length(freqtable)] <- al
        }
    }

    # return data frame
    return(freqtable)
}

# Internal functions from De Silva et al. 2005, for use in deSilvaFreq
# and meandistance.matrix2.

    # G function from de Silva et al
    .G <- function(q, n){
        return(factorial(n+q)/(factorial(q+1)*factorial(n-1)))
    }

    # INDEXG function from de Silva et al
    .indexg <- function(ag1, na1, m2){
        x <- 1 + ag1[m2] - ag1[m2-1]
        for(q in 1:(m2-2)){
            x <- x + .G(q,na1+1-ag1[m2-q-1]) - .G(q,na1+1-ag1[m2-q])
        }
        x <- x + .G(m2-1,na1) - .G(m2-1,na1+1-ag1[1])
        return(x)
    }

    # GENLIST function from de Silva et al
    .genlist <- function(ng, na1, m2){
        # set up temporary genotype vector and genotype array
        ag1 <- rep(1, m2)
        ag <- array(0, dim=c(ng, m2))

        # fill genotype array with all possible combinations
        g <- 1
        ag[g,] <- ag1
        a <- m2
        while(a>0){
            if(ag1[a]==na1){
                a <- a-1
            } else {
                if(a > 0){
                    ag1[a] <- ag1[a] + 1
                    if(a < m2){
                        for(a1 in (a+1):m2){
                            ag1[a1] <- ag1[a]
                        }
                    }
                    g <- g+1
                    ag[g,] <- ag1
                    a <- m2
                }
            }
        }
        return(ag)
    }

    # RANMUL function
    .ranmul <- function(ng, na1, ag, m2){
        # RANMUL subroutine
        # rmul is multiplier to get genotype freq under random mating
        # arep shows how many copies of each allele each genotype has
        rmul <- rep(1, ng)
        arep <- matrix(0, nrow=ng, ncol=na1)
        for(g in 1:ng){
            ag1 <- ag[g,]
            arep[g,ag1[1]] <- 1
            for(j in 2:m2){
                rmul[g] <- rmul[g]*j
                if(ag1[j] == ag1[j-1]){
                    arep[g,ag1[j]] <- arep[g,ag1[j]] + 1
                    rmul[g] <- rmul[g]/arep[g,ag1[j]]
                } else {
                    arep[g,ag1[j]] <- 1
                }
            }
        }
        return(list(rmul, arep))
    }

    # SELFMAT function
    .selfmat <- function(ng, na1, ag, m2){
        # SELFMAT subroutine
        # smat is selfing matrix, and smatdiv is divisor for selfing matrix
        m <- m2/2
        smat <- matrix(0, nrow=ng, ncol=ng)
        al1 <- rep(0, m)
        al2 <- rep(0, m)
        ag1 <- rep(0, m2)
        for(g in 1:ng){
            al1[1] <- 1
            if(m > 1){
                for(j in 2:m) al1[j] <- al1[j-1] + 1
            }
            al1[m] <- al1[m] - 1
            a1 <- m
            while(a1 > 0){
                if(al1[a1] == (m+a1)){
                    a1 <- a1 - 1
                } else {
                    if(a1 > 0){
                        al1[a1] <- al1[a1] + 1
                        if(a1 < m){
                            for(a3 in (a1+1):m) al1[a3] <- al1[a3-1]+1
                        }
                        al2[1] <- 1
                        if(m > 1){
                            for(j in 2:m) al2[j] <- al2[j-1] + 1
                        }
                        al2[m] <- al2[m]-1
                        a2 <- m
                        while(a2 > 0){
                            if(al2[a2] == (m+a2)){
                                a2 <- a2 - 1
                            } else {
                                if(a2 > 0){
                                    al2[a2] <- al2[a2] + 1
                                    if(a2 < m){
                                        for(a3 in (a2+1):m) al2[a3]<-al2[a3-1]+1
                                    }
                                    #UPDATESMAT subroutine
                                    j1 <- 1
                                    j2 <- 1
                                    k1 <- al1[j1]
                                    k2 <- al2[j2]
                                    for(i in 1:m2){
                                        if(k1 < k2){
                                            ag1[i] <- ag[g,k1]
                                            if(j1 == m){
                                                k1 <- 999
                                            } else {
                                                j1 <- j1 + 1
                                                k1 <- al1[j1]
                                            }
                                        } else {
                                            ag1[i] <- ag[g,k2]
                                            if(j2 == m){
                                                k2 <- 999
                                            } else {
                                                j2 <- j2 + 1
                                                k2 <- al2[j2]
                                            }
                                        }
                                    }
                                    g2 <- .indexg(ag1, na1, m2)
                                    smat[g,g2] <- smat[g,g2] + 1
                                    # end UPDATESMAT subroutine
                                    a2 <- m
                                }
                            }
                        }
                        a1 <- m
                    }
                }
            }
        }
        return(smat)
    }


# Iterative estimation of allele frequencies under polysomic inheritance,
# with a uniform, even-numbered ploidy and a known selfing rate.
# Much of this code is translated directly from:
# De Silva HN, Hall AJ, Rikkerink E, McNeilage MA, and LG Fraser (2005).
# Estimation of allele frequencies in polyploids under certain patterns
# of inheritance.  Heredity 95:327-334.
deSilvaFreq <- function(object, self,
                        samples=Samples(object),
                        loci=Loci(object), initNull = 0.15,
                        initFreq=simpleFreq(object[samples,loci]),
                        tol = 0.00000001){
    # make sure self argument has been given
    if(missing(self)) stop("Selfing rate required.")

    # make sure initFreq is in right format
    if(!"Genomes" %in% names(initFreq))
      stop("initFreq must have single Genomes column.")

    # convert object to genambig if necessary
    if("genbinary" %in% class(object)) object <- genbinary.to.genambig(object,
                                                                       samples,
                                                                       loci)
    
    # get the ploidy (m2), and check that there is only one and that it is even
    m2 <- unique(Ploidies(object,samples,loci))
    if(length(m2) != 1)
      stop("Only one ploidy allowed.  Try running subsets of data one ploidy at a time.")
    if(is.na(m2)) stop("Function requires information in Ploidies slot.")
    if(m2 %% 2 != 0) stop("Ploidy must be even.")

    # check and set up initNull
    if(!length(initNull) %in% c(1,length(loci)))
        stop("Need single value for initNull or one value per locus.")
    if(length(initNull) == 1)
        initNull <- rep(initNull, times=length(loci))
    if(is.null(names(initNull)))
        names(initNull) <- loci

    # calculate the number of alleles from one gamete
    # (from de Silva et al's INIT subroutine)
    m <- m2/2

    # get the populations that will be used
    pops <- PopNames(object)[unique(PopInfo(object)[samples])]
    if(!identical(pops, row.names(initFreq)))
        stop("Population names in initFreq don't match those in object.")

    # set up the data frame where final allele frequencies will be stored
    # has columns from initFreq, plus columns for nulls
    finalfreq <- data.frame(row.names=pops, Genomes=initFreq$Genomes,
                            matrix(0, nrow=length(pops),
                                   ncol=dim(initFreq)[2]-1+length(loci),
                                   dimnames=list(NULL, sort(c(
                                   paste(loci, ".null", sep=""),
                                   names(initFreq)[2:dim(initFreq)[2]])))))

    # get the divisor for the selfing matrices
    smatdiv <- (.G(m-1,m+1))^2

    # INDEXF function from de Silva et al
    # af1 = vector of alleles in phenotype
    # m = number of observable alleles in phenotype
    # na = number of alleles, calculated for each locus
    indexf <- function(m, af1, na){
        x <- 1
        if(m == 1){
            x <- x + af1[1]
        } else {
            if ( m >1 ){
                for(q in 1:(m-1)){
                    x <- x + .G(q-1,na-q+1)
                }
                x <- x + .G(m-1,na+1-m) - .G(m-1,na+2-m-af1[1])
                if(m > 2){
                    for(q in 1:(m-2)){
                        x <- x + .G(q,na-q-af1[m-q-1]) -
                            .G(q,na+1-q-af1[m-q])
                    }
                }
                x <- x + af1[m] - af1[m-1]
            }
        }
        return(x)
    }

    # FENLIST function
    fenlist <- function(na){
        # set up temporary holding vector for phenotpes (FENLIST)
        af1 <- rep(0, m2)
        # set up array to contain all phenotypes (FENLIST)
        af <- array(0, dim=c(1, m2))
        # set up vector for number of alleles in each phenotype (FENLIST)
        naf <- integer(0)

        # fill in af and naf (FENLIST)
        # This is done with rbind rather than setting up whole array,
        # because the method for calculating the number of phenotypes
        # doesn't seem to work for certain numbers of alleles.
        f <- 1
        naf <- 0
        for(m in 1:min(m2, na)){
            af1[1] <- 1
            if(m > 1){
                for(j in 2:m){
                    af1[j] <- af1[j-1] + 1
                }
            }
            f <- f + 1
            naf[f] <-  m
            af <- rbind(af, af1)
            a <- m
            while(a > 0){
                if(af1[a] == (na+a-m)){
                    a <- a - 1
                } else {
                    if(a > 0){
                        af1[a] <- af1[a] + 1
                        if(a < m){
                            for(a1 in (a+1):m) af1[a1] <- af1[a1-1] + 1
                        }
                        f <- f + 1
                        naf[f] <- m
                        af <- rbind(af, af1)
                        a <- m
                    }
                }
            }
        }

        return(list(af, naf))
    }

    # CONVMAT function
    convmat <- function(ng, nf, na1, ag){
        na <- na1-1
        af1 <- rep(0, m2)
        # CONVMAT subroutine - make a matrix for conversion from
        # genotypes to phenotypes
        cmat <- matrix(0, nrow=nf, ncol=ng)
        for(g in 1:ng){
            ag1 <- ag[g,]
            # find the phenotype (af1) that matches this genotype (ag1)
            if(ag1[1] == na1){
                naf1 <- 0 # this is the homozygous null genotype
            } else {
                naf1 <- 1
                af1[naf1] <- ag1[1]
                for(a in 2:m2){
                    if(ag1[a] == na1) break # exit loop if null allele
                    if(ag1[a] > af1[naf1]){
                        naf1 <- naf1 + 1
                        af1[naf1] <- ag1[a]
                    }
                }
            }
            # fill in the extra alleles with zeros
            if(naf1 < m2){
                for(j in (naf1 + 1):m2){
                    af1[j] <- 0
                }
            }
            f <- indexf(naf1, af1, na) # This is the one phenotype to match
                                       # this genotype.
            cmat[f,g] <- 1
        }
        return(cmat)
    }

    # Loop to do calculations one locus and population at a time
    for(L in loci){
        cat(paste("Starting", L), sep="\n")

        ## Begin looping through populations
        for(pop in pops){
            cat(paste("Starting", L, pop), sep="\n")
            # get a list of samples in this pop being analyzed
            psamples <- Samples(object, populations=pop)[!isMissing(object,
                                        Samples(object, populations=pop),
                                        L)]
            psamples <- psamples[psamples %in% samples]

            # extract the initial allele frequencies and add a null
            subInitFreq <- initFreq[pop, grep(paste(L,".",sep=""),
                                              names(initFreq), fixed=TRUE)]

            # set up matrices if this has not already been done
                templist <- names(subInitFreq)[subInitFreq !=0]
                templist <- strsplit(templist, split=".", fixed=TRUE)
                alleles <- rep(0, length(templist))
                for(i in 1:length(alleles)){
                    alleles[i] <- as.integer(templist[[i]][2])
                }
            alleles <- sort(alleles)
                na <- length(alleles)

                # get the number of alleles with a null
                na1 <- na + 1
                # get the number of genotypes (from the GENLIST function)
                ng <- na1
                for(j in 2:m2){
                    ng <- ng*(na1+j-1)/j
                }

            ag <- .genlist(ng, na1, m2)
            temp <- fenlist(na)
            af <- temp[[1]]
            naf <- temp[[2]]
            nf <- length(naf)
            temp <- .ranmul(ng, na1, ag, m2)
            rmul <- temp[[1]]
            arep <- temp[[2]]
            smatt <- .selfmat(ng, na1, ag, m2)/smatdiv
            cmat <- convmat(ng, nf, na1, ag)
            rm(temp)

            # calculate pp, the frequency of each phenotype in this population
            pp <- rep(0, nf)
            for(s in psamples){
                phenotype <- sort(unique(Genotype(object, s, L)))
                phenotype <- match(phenotype, alleles)
                f <- indexf(length(phenotype), phenotype, na)
                pp[f] <- pp[f] + 1
            }
            pp <- pp/sum(pp)

            # Get initial allele frequencies
            p1 <- rep(0, na1) # vector to hold frequencies
            p1[na1] <- initNull[L] # add null freq to the last position
            # make sure everything will sum to 1
            subInitFreq <- subInitFreq * (1-initNull[L])/sum(subInitFreq)
            # get each allele frequency
            for(a in alleles){
                p1[match(a, alleles)] <- subInitFreq[1,
                                                     grep(a, names(subInitFreq),
                                                          fixed=TRUE)]
            }

            ## Begin the EM algorithm
            converge <- 0
            niter <- 1
            oneg <- rep(1, ng)
            while(converge == 0){
                # Expectation step
                # GPROBS subroutine
                pa <- rep(0, na1)
                pa[na1] <- 1
                for(j in 1:na){
                    pa[j] <- p1[j]
                    pa[na1] <- pa[na1]-p1[j]
                    # it seems like this is the same as pa <- p1
                    # unless p1 is not supposed to have the null allele
                    # (conflicting information in orignal code comments)
                }
                rvec <- rep(0, ng)
                for(g in 1:ng){
                    rvec[g] <- rmul[g]
                    for(j in 1:m2){
                        rvec[g] <- rvec[g]*pa[ag[g,j]]
                    }
                    # This gives prob of genotype by multiplying
                    # probs of alleles and coefficients for a multinomial
                    # distribution.
                }

                # make an identity matrix
                id <- diag(nrow=ng)
                # calculate gprob using selfing and outcrossing matrices
                s3 <- id - self * smatt
                s3inv <- solve(s3)
                gprob <- (1-self) * s3inv %*% rvec
                # end GPROBS

                # equation (12) from the paper
                xx1 <- matrix(0, nrow=nf, ncol=ng)
                for(i in 1:nf){
                    xx1[i,] <- cmat[i,] * gprob
                }
                xx2 <- xx1 %*% oneg
                xx3 <- matrix(0, nrow=nf, ncol=ng)
                for(i in 1:ng){
                    xx3[,i] <- xx1[,i] * (xx2^(-1))
                }
                EP <- t(xx3) %*% pp

                # Maximization step
                p2 <- t(arep) %*% EP/m2

                # check for convergence
                pB <- p1 + p2
                pT <- p1 - p2
                pT <- pT[pB != 0]
                pB <- pB[pB != 0]

                if(length(pB) == 0){
                    converge <- 1
                } else {
                    if(sum(abs(pT)/pB) <= tol){
                        converge <- 1
                    }
                }

                niter <- niter + 1
                p1 <- p2
            }

            # write frequencies in p2 to finalfreqs
            for(a in alleles){
                finalfreq[pop, match(paste(L, a, sep="."), names(finalfreq))]<-
                    p2[match(a, alleles)]
            }
            finalfreq[pop, match(paste(L, "null", sep="."), names(finalfreq))]<-
                p2[na1]
            # print the number of reps
            cat(paste(niter-1, "repetitions for", L, pop),
                 sep="\n")
        }
    }
    return(finalfreq)
}

calcFst<-function(freqs, pops=row.names(freqs),
                  loci=unique(as.matrix(as.data.frame(strsplit(names(freqs),
                  split=".",
                                                      fixed=TRUE),
                                             stringsAsFactors=FALSE))[1,])){
    # Clean up loci
    loci<-loci[loci!="Genomes"]
    # Set up matrix for Fst values
    fsts<-matrix(0,nrow=length(pops),ncol=length(pops),dimnames=list(pops,pops))
    # Get genome number from the table
    if("Genomes" %in% names(freqs)){
      genomes<-freqs$Genomes
      names(genomes)<-pops
      GbyL <- FALSE
    } else {
      GbyL <- TRUE
    }
    for(m in 1:length(pops)){
        for(n in m:length(pops)){
            # set up array for HT and HS values
            hets<-array(0,dim=c(length(loci),2),
                        dimnames=list(loci,c("HT","HS")))
            if(!GbyL){
              genomesM <- genomes[pops[m]]
              genomesN <- genomes[pops[n]]
            }
            for(L in loci){
              if(GbyL){
                genomesM <- freqs[pops[m],paste(L, "Genomes", sep=".")]
                genomesN <- freqs[pops[n],paste(L, "Genomes", sep=".")]
              }
              # get just the frequencies for these pops and this locus
              thesefreqs<-freqs[c(pops[m],pops[n]),
                                grep(paste(L,".",sep=""),
                                     names(freqs),fixed=TRUE)]
              thesefreqs <-
                thesefreqs[,names(thesefreqs)!=paste(L,"Genomes",sep=".")]
              # get average allele frequencies weighted by genomes/pop
              avgfreq<-(thesefreqs[1,]*genomesM +
                        thesefreqs[2,]*genomesN)/
                  (genomesM + genomesN)
              # estimate H by 1 - sum of squared allele frequencies
              # put the heterozygositites in the array
              hets[L,"HT"]<-1-sum(avgfreq^2)
              hets[L,"HS"]<-((1-sum(thesefreqs[1,]^2))*genomesM +
                             (1-sum(thesefreqs[2,]^2))*genomesN)/
                                 (genomesM + genomesN)
            }
            HT<-mean(hets[,"HT"])
            HS<-mean(hets[,"HS"])
            fsts[m,n]<-(HT-HS)/HT
            fsts[n,m]<-(HT-HS)/HT
        }
    }
    # return matrix of Fst values
    return(fsts)
}

# put allele frequencies into a format for SPAGeDi for it to use in estimates
# of kinship and relatedness coefficients
write.freq.SPAGeDi <- function(freqs, usatnts, file="", digits=2,
                               pops=row.names(freqs),
                  loci=unique(as.matrix(as.data.frame(strsplit(names(freqs),
                  split=".", fixed=TRUE), stringsAsFactors=FALSE))[1,])){
    if(is.null(names(usatnts))) names(usatnts) <- loci
    loci <- loci[loci != "Genomes"]
    # subset the data frame
    freqs <- freqs[pops,]
    # make a list to contain the columns to write
    datalist <- list(0)
    length(datalist) <- ( length(loci) * 2 )
    item <- 1
    genbyloc <- !"Genomes" %in% names(freqs)
    if(!genbyloc) genomes <- freqs$Genomes

    for(L in loci){
        # find the alleles
        alleles <- as.matrix(as.data.frame(strsplit(names(freqs)[grep(paste(L,".",sep=""),
                                                            names(freqs),
                                                            fixed=TRUE)],
                                          split=".", fixed=TRUE))[2,])
        alleles <- as.vector(alleles)
        if(genbyloc){
          alleles <- alleles[alleles!="Genomes"]
        }
        # convert alleles to numbers used in SPAGeDi file
        alleles[alleles == "null"] <- 0
        alleles <- as.integer(alleles)
        alleles <- floor(alleles/usatnts[L])
        suballele<-function(value){if(value[1]!=0){
            value-(10^(digits-1))} else{
            0}}
        while(max(alleles) >= 10^digits){
            alleles<-mapply(suballele, alleles)
        }
        # add alleles to datalist
        datalist[[item]] <- alleles
        names(datalist)[item] <- L
        item <- item + 1
        # make a weighted average of allele frequencies across populations
        if(genbyloc){
          genomes <- freqs[[paste(L, "Genomes", sep=".")]]
          locindex <- grep(paste(L,".",sep=""), names(freqs),fixed=TRUE)
          locindex <- locindex[!locindex %in% grep("Genomes", names(freqs))]
        } else {
          locindex <- grep(paste(L,".",sep=""), names(freqs),fixed=TRUE)
        }
        avgfreq <- (genomes %*% as.matrix(freqs[,locindex])) / sum(genomes)
        # add frequencies to list
        datalist[[item]] <- as.vector(avgfreq)
        names(datalist)[item] <- length(avgfreq)
        item <- item + 1
    }

    # find the maximum number of alleles
    maxal <- max(mapply(length, datalist))
    # set up data frame to write
    fr <- data.frame(row.names=1:maxal)
    # put list elements into the data frame
    for(i in datalist){
        fr <- data.frame(fr, c(i, rep("", maxal-length(i))))
    }
    names(fr) <- names(datalist)

    # write the file
    write.table(fr, file=file, sep="\t", col.names=TRUE, row.names=FALSE,
                quote=FALSE)
}

freq.to.genpop <- function(freqs, pops=row.names(freqs),
            loci=unique(as.matrix(as.data.frame(strsplit(names(freqs),
                  split=".",fixed=TRUE),stringsAsFactors=FALSE))[1,])){
    # clean up loci
    loci <- loci[loci!="Genomes"]

    # errors
    if(!"Genomes" %in% names(freqs))
      stop("Only one Genomes column allowed.")

    # get population sizes
    popsizes <- freqs[pops, "Genomes"]

    # get an index of columns containing these loci
    loccol <- integer(0)
    for(L in loci){
        loccol <- c(loccol, grep(paste(L,".",sep=""), names(freqs), fixed=TRUE))
    }
    # get a table just for these populations and loci
    subfreq <- freqs[pops, loccol]

    # convert allele frequencies to allele counts
    for(i in 1:length(subfreq)){
        subfreq[[i]] <- round(popsizes * subfreq[[i]])
    }

    return(subfreq)
}

## Genotype diversity functions
# p is a vector of genotype counts
# base = exp(1) for natural log, or base = 2 for log base 2
Shannon <- function(p, base=exp(1)){
    N <- sum(p)
    return(-sum(p/N * log(p/N, base=base)))
}

Simpson <- function(p){
    N <- sum(p)
    return(sum(p*(p-1)/(N*(N-1))))
}

Simpson.var <- function(p){
    N <- sum(p)
    p <- p/N
    v <- (4*N*(N-1)*(N-2)*sum(p^3) + 2*N*(N-1)*sum(p^2) -
            2*N*(N-1)*(2*N-3)*(sum(p^2))^2)/((N*(N-1))^2)
    return(v)
}

# ... is passed to index. (So you can adjust base of Shannon index.)
# threshold is highest distance between two individuals that can be considered
# to be one clone.
genotypeDiversity <- function(genobject,
                   samples = Samples(genobject), loci = Loci(genobject),
                    d = meandistance.matrix(genobject,samples,loci,
                    all.distances=TRUE,distmetric=Lynch.distance),
                    threshold = 0, index = Shannon,
                              ...){
    # get populations
    if(all(is.na(PopInfo(genobject)[samples]))){
        PopInfo(genobject)[samples] <- rep(1, length(samples))
    }
    if(!all(!is.na(PopInfo(genobject)[samples]))){
        stop("PopInfo needed.")
    }
    pops <- PopNames(genobject)[unique(PopInfo(genobject)[samples])]

    # set up results matrix
    results <- matrix(NA, ncol=length(loci)+1, nrow=length(pops),
                      dimnames=list(pops, c(loci,"overall")))

    # statistics for individual loci
    for(L in loci){
        for(p in pops){
            # get all samples for this pop with data for this locus
            psamples <- samples[samples %in% Samples(genobject, populations=p)]
            psamples <- psamples[!isMissing(genobject, psamples, L)]

            # assign to genotype groups
            dsub <- d[[1]][L, psamples, psamples]
            if(is.vector(dsub)){ # fix for if there was only one sample
                dsub <- array(dsub, dim=c(1,1),
                              dimnames=list(psamples, psamples))
            }
            clones <- assignClones(dsub,
                                   threshold=threshold)
            # get genotype frequencies
#            n <- length(clones) # number of individuals
            cl <- length(unique(clones)) # number of groups
            counts <- rep(NA, cl) # vector to hold counts of individuals
            for(i in 1:cl){
                counts[i] <- length(clones[clones == i])
            }
            # get diversity index
            results[p,L] <- index(counts, ...)
        }
    }

    # get a mean array for only the loci being used here
    if(length(loci) > dim(d[[1]])[1]){
        d2 <- meandist.from.array(d[[1]], samples=samples, loci=loci)
    } else {
        d2 <- d[[2]]
    }
    # statistics for multilocus genotypes
    for(p in pops){
        psamples <- samples[samples %in% Samples(genobject, populations=p)]
        clones <- assignClones(d2, psamples, threshold)
#        n <- length(clones) # number of individuals
        cl <- length(unique(clones)) # number of groups
        counts <- rep(NA, cl) # vector to hold counts of individuals
        for(i in 1:cl){
            counts[i] <- length(clones[clones == i])
        }
        # get diversity index
        results[p,"overall"] <- index(counts, ...)
    }

    return(results)
}

# function to get unique alleles and counts
alleleDiversity <- function(genobject, samples=Samples(genobject),
                            loci=Loci(genobject),
                            alleles=TRUE, counts=TRUE){
  if(!alleles && !counts)
    stop("At least one of alleles and counts must be set to TRUE")
  # get populations
    if(all(is.na(PopInfo(genobject)[samples]))){
        PopInfo(genobject)[samples] <- rep(1, length(samples))
    }
    if(!all(!is.na(PopInfo(genobject)[samples]))){
        stop("PopInfo needed.")
    }
    pops <- PopNames(genobject)[unique(PopInfo(genobject)[samples])]

    # set up results tables
    rcounts <- matrix(NA, nrow=length(pops)+1, ncol=length(loci),
                     dimnames=list(c(pops,"overall"),loci))
    ralleles <- array(list(NA), dim=c(length(pops)+1,length(loci)),
                      dimnames=dimnames(rcounts))
    # get unique alleles
    for(L in loci){
      for(p in pops){
        psamples <- samples[samples %in% Samples(genobject, populations=p)]
        xalleles <- .unal1loc(genobject, psamples, L)
        ralleles[[p,L]] <- xalleles
        rcounts[p,L] <- length(xalleles)
      }
      xalleles <- .unal1loc(genobject, samples, L)
      ralleles[["overall",L]] <- xalleles
      rcounts["overall",L] <- length(xalleles)
    }

    if(alleles && counts)
      return(list(alleles=ralleles, counts=rcounts))
  if(alleles && !counts)
    return(ralleles)
  if(!alleles && counts)
    return(rcounts)
}
