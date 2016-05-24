Bruvo.distance <- function(genotype1, genotype2, maxl=9, usatnt=2, missing=-9) {
    
if(is.na(usatnt)) stop("Bruvo.distance needs info from Usatnts slot.")

if(identical(genotype1, genotype2)&&genotype1[1]!=missing) {dist <- 0} else {
    if((length(genotype1)>maxl & length(genotype2)>maxl)|
       genotype1[1]==missing|genotype2[1]==missing){
        dist <- NA} else {
        if(length(genotype1) >= length(genotype2)) {
            genotypeL <- genotype1/usatnt; genotypeS <- genotype2/usatnt} else {
                genotypeL <- genotype2/usatnt; genotypeS <- genotype1/usatnt}

# if genotypes are identical, just return a distance of zero without doing the rest
        # of the calculation
# if genotypes are both longer than 9, skip this calculation because it will take hours
# whichever genotype has more alleles, make this genotypeL (long) and the other
        # genotypeS (short)
# convert alleles into repeat counts by dividing by usatnt

kl <- length(genotypeL) # sets the ploidy level for this genotype comparison
ks <- length(genotypeS) # number of alleles in the shorter genotype

allele.distances <- array(0 , c(kl,ks))
# Create an empty matrix to contain the raw distances between alleles

for(n in 1:kl) { for(m in 1:ks) {allele.distances[n,m] <- genotypeL[n] - genotypeS[m]}}
# fills the array with the differences in allele repeat count

geometric.distances <- array(1 - 2^-abs(allele.distances) , c(kl,ks))
# geometric transformation based on mutation probabilities

#Next, find the minimum distance sum among all permutations

column <- 1:ks # an index of all columns (genotypeS alleles)
row <- 1:kl # an index of all rows (genotypeL alleles)

combinations <- combn(row, ks, FUN = NULL, simplify=FALSE)
# all combinations of alleles in genotypeL that can be matched to non-virtual
        # alleles in genotypeS

permutations <- combinat::permn(ks)
# all possible orders that alleles within these combinations can go in

mindist <- Inf # this variable will store the minimum sum encountered so far.

for(i in 1:length(combinations)) {
# the loop to go through every possible sum of compatible allele comparisons

rowcomb <- combinations[[i]] # choose one combination of rows for this round

for(l in 1:length(permutations)){ # go through all orders of this combinations of rows

sum <- 0 # this is si, the sum of allele comparisons

for(j in 1:ks){
    sum <- sum + geometric.distances[rowcomb[permutations[[l]][j]],column[j]]}
# the loop to calculate the sum for this permutation

if(sum < mindist) {mindist <- sum} # is this the minimum sum found so far?

}}

dist <- (mindist+kl-ks)/kl
# add 1 for each infinite virtual allele, then divide by the ploidy

}}
return(dist)
}

Lynch.distance<-function(genotype1,genotype2,usatnt=NA,missing=-9){
    if(genotype1[1]==missing || genotype2[1]==missing){
        # return NA if there is any missing data
        distance<-NA
    } else {
        # make sure each allele is only listed once
        genotype1 <- unique(genotype1)
        genotype2 <- unique(genotype2)
        # get the average number of bands for the two genotypes
        meanbands<-(length(genotype1)+length(genotype2))/2
        # find how many bands the genotypes have in common
        commonbands<- sum(genotype1 %in% genotype2)
        # calculate the distance
        distance<- 1-(commonbands/meanbands)
    }

    # return distance
    return(distance)
}


meandistance.matrix <- function(object, samples=Samples(object),
                                 loci=Loci(object), all.distances=FALSE,
                                 distmetric=Bruvo.distance,
                                progress=TRUE, ...){
    # error if not a genambig object
    if(!is(object, "genambig")) stop("\"genambig\" object needed.")
    # subset the object so that samples can be numbered
    object <- object[samples, loci]

# create an array containing all distances by locus and sample
    loci.matrices<-array(dim=c(length(loci),length(samples),length(samples)),
                         dimnames=list(loci,samples,samples))
    for(L in loci){
       for(m in 1:length(samples)){
           for(n in m:length(samples)){
               thisdistance <- distmetric(Genotype(object, m, L),
                                          Genotype(object, n, L),
                                          usatnt = Usatnts(object)[L],
                                          missing = Missing(object),
                                          ...)
               loci.matrices[L,m,n] <- thisdistance
               loci.matrices[L,n,m] <- thisdistance
               if(progress) print(c(L, samples[m], samples[n]))
           }
       }
    }

    # calculate the mean matrix across all loci
    mean.matrix <- meandist.from.array(loci.matrices)

    # return either the mean matrix and possibly the array as well
    if(all.distances){
        return(list(DistByLoc=loci.matrices, MeanMatrix=mean.matrix))
    } else {
        return(mean.matrix)
    }
}

meandist.from.array<-function(distarray, samples=dimnames(distarray)[[2]],
                              loci=dimnames(distarray)[[1]]){
    # get the array to be averaged
    subarray<-distarray[loci,samples,samples]
    # calculate means
    mean.matrix <- colMeans(subarray, na.rm=TRUE)
    return(mean.matrix)
}

find.na.dist<-function(distarray, samples=dimnames(distarray)[[2]],
                              loci=dimnames(distarray)[[1]]){
    # number of missing distances
    nmissing <- sum(is.na(distarray[loci, samples, samples]))
    # set up vectors for data frame to contain info on where missing data is
    Locus <- character(nmissing)
    Sample1 <- character(nmissing)
    Sample2 <- character(nmissing)
    # current row in the data frame
    currrow<-1
    # go through the array, find NA, and put the index into the data frame
    for(L in loci){
        for(s1 in samples){
            for(s2 in samples){
                if(is.na(distarray[L,s1,s2])){
                    Locus[currrow]<-L
                    Sample1[currrow]<-s1
                    Sample2[currrow]<-s2
                    currrow<-currrow+1
                }
            }
        }
    }
    #return data frame
    return(data.frame(Locus=Locus, Sample1=Sample1, Sample2=Sample2, stringsAsFactors=FALSE))
}

find.missing.gen<-function(object, samples=Samples(object),
                           loci=Loci(object)){
    # set up vectors to contain the indices to put into the data frame
    Locus <- c("")
    Sample <- c("")
    # current row in the data frame
    currrow<-1
    # find which data are missing
    for(L in loci){
        for(s in samples){
            if(isMissing(object, s, L)){
                Locus[currrow]<-L
                Sample[currrow]<-s
                currrow<-currrow+1
            }
        }
    }
    # return the data frame
    return(data.frame(Locus=Locus,Sample=Sample, stringsAsFactors=FALSE))
}

# find NA distances that aren't the result of missing data
find.na.dist.not.missing<-function(object, distarray,
                                   samples=dimnames(distarray)[[2]],
                                   loci=dimnames(distarray)[[1]]){
    # get the data frames of NA distances and missing genotypes
    na.dist<-find.na.dist(distarray,samples=samples,loci=loci)
    missing.gen<-find.missing.gen(object,samples=samples,loci=loci)
    # set up vectors for data frame to contain results
    Locus<-""
    Sample1<-""
    Sample2<-""
    currrow<-1
    # for each row of na.dist, look for that locus and samples in missing.gen
    for(i in 1:length(na.dist$Locus)){
        L<-na.dist$Locus[i]
        s1<-na.dist$Sample1[i]
        s2<-na.dist$Sample2[i]
        missthislocus<-missing.gen[missing.gen$Locus==L,]
        if(identical(c(s1,s2) %in% missthislocus$Sample, c(FALSE,FALSE))){
            Locus[currrow]<-L
            Sample1[currrow]<-s1
            Sample2[currrow]<-s2
            currrow<-currrow+1
        }
    }
    # return data frame
    return(data.frame(Locus=Locus, Sample1=Sample1, Sample2=Sample2,
                      stringsAsFactors=FALSE))
}

genotypeProbs <- function(object, sample, locus, freq=NULL, gprob=NULL,
                          alleles=NULL){
    pl <- Ploidies(object,sample,locus) # get ploidy
    gen <- Genotype(object, sample, locus) # get ambiguous genotype

    # determine which method we are using
    u <- c(is.null(freq),is.null(gprob),is.null(alleles))
    ok <- FALSE
    if(identical(u, c(FALSE,TRUE,TRUE))){
        selfing <- FALSE
        ok <- TRUE
    }
    if(identical(u, c(TRUE, FALSE,FALSE))){
        selfing <- TRUE
        ok <- TRUE
    }
    if(!ok) stop("Supply freq, or gprob and alleles.")

    # Errors
    if(is.na(pl)) stop("Ploidies required.")
    if(length(gen) >  pl)
        stop(paste("There are too many alleles for ploidy:",sample,locus))

    # What to do with unambiguous genotypes
    if(length(gen)==1){
        if(pl==0 && isMissing(object, sample, locus)) pl <- 1
               # 8/30/14 - above edit is in anticipation of recodeAllopoly
        results <- list(probs=1, genotypes=matrix(rep(gen, pl), nrow=1,
                                 ncol=pl))
    }
    if(length(gen)==pl){
        results <- list(probs=1, genotypes=matrix(gen, nrow=1, ncol=pl))
    }

    # What to do with ambiguous genotypes
    if(!length(gen) %in% c(1,pl)){
        if(!selfing){
          pop <- PopNames(object)[PopInfo(object)[sample]]
          f <- freq[pop, paste(locus, gen, sep=".")] # get allele frequencies
          f <- as.vector(f/sum(f)) # normalize frequencies
          names(f) <- gen
        }
        a1 <- length(gen) # number of alleles
        a2 <- pl-a1 # number of slots to fill

        # set up function for recursive building of genotypes
        makegen <- function(gmat){
            gmat2 <- matrix(nrow=0, ncol=dim(gmat)[2]+1)
            for(i in 1:dim(gmat)[1]){
                alsofar <- gmat[i,]
                lastallele <- alsofar[length(alsofar)]
                for(j in gen[match(lastallele,gen):length(gen)]){
                    gmat2 <- rbind(gmat2, matrix(c(alsofar,j),nrow=1,
                                                 ncol=dim(gmat2)[2]))
                }
            }
            if(dim(gmat2)[2]==a2){
                return(gmat2)
            } else {
                return(makegen(gmat2))
            }
        }
        # alleles in ambiguous slots
        allelefill <- matrix(gen, nrow=a1, ncol=1)
        if(a2 > 1) allelefill <- makegen(allelefill)
        ng <- dim(allelefill)[1] # number of genotypes
        results <- list(probs=c(), genotypes=matrix())
        # genotypes (unambiguous + ambiguous alleles)
        results[["genotypes"]] <- cbind(matrix(rep(gen,each=ng),
                                               nrow=ng,ncol=a1),
                                        allelefill)
        # Probability of each genotype & allele sorting
        for(i in 1:ng){
            # sort the alleles in this genotype
            results$genotypes[i,] <- sort(results$genotypes[i,])
            if(!selfing){  ## random mating method
              # multiply allele probabilities
              fmult <- prod(f[as.character(allelefill[i,])])
              # get the polynomial coefficient
              allelecopies <- c()
              for(a in unique(allelefill[i,])){
                  allelecopies <- c(allelecopies,
                                  length(allelefill[i,][allelefill[i,]==a]))
              }
              pc <- factorial(a2)/prod(mapply(factorial, allelecopies))
              # get the probability of this genotype
              results[["probs"]][i] <- pc * fmult
            } else {  ## partial selfing method
                # get unambig genotype in terms of genlist numbers
                ugen <- match(as.integer(results$genotypes[i,]),alleles)
                results$probs[i] <- gprob[.indexg(ugen, length(alleles), pl)]
            }
        }
        if(selfing){
            # with selfing method, need to normalize probabilities.
            results$probs <- results$probs/sum(results$probs)
        }
    }

    # Return a list containing unambiguous genotypes and their probabilities
    return(results)
}

meandistance.matrix2 <- function(object, samples=Samples(object),
                                 loci=Loci(object),
                            freq=simpleFreq(object, samples, loci), self=0,
                            all.distances=FALSE, distmetric = Bruvo.distance,
                            progress=TRUE, ...){
    # Errors
    if(!all(!is.na(Ploidies(object,samples,loci))))
        stop("Ploidies needed.")
    if(!all(!is.na(PopInfo(object)[samples])))
        stop("PopInfo needed.")
    if(!all(PopNames(object)[unique(PopInfo(object)[samples])] %in%
            row.names(freq))) stop("PopNames must match row.names of freq.")

    # calculate probabilities of every unambiguous genotype from all genotypes
    gprobs <- array(list(NA), dim=c(length(samples), length(loci)),
                    dimnames=list(samples, loci))
    if(self == 0){  # calculate genotype probabilities under random mating
        for(s in samples){
            for(L in loci){
                gprobs[[s,L]] <- genotypeProbs(object, s, L, freq=freq)
            }
        }
    } else {  # calculate genotype probabilities under partial selfing
        pops <- PopNames(object)[unique(PopInfo(object)[samples])]
        for(p in pops){
            psamples <- samples[samples %in% Samples(object, populations=p)]
            
            for(L in loci){
                m2 <- unique(Ploidies(object, psamples, L))
                if(length(m2) != 1)
                  stop("Only one ploidy allowed per pop*locus when self > 0.")
                if(is.na(m2))
                    stop("Function requires information in Ploidies slot.")
                if(m2 %% 2 != 0) stop("Ploidy must be even.")
                
                cat("Setting up genotype probabilities...",sep="\n")
                subfreq <- freq[p, grep(L, names(freq), fixed=TRUE)]
                templist <- names(subfreq)[subfreq !=0]
                templist <- strsplit(templist, split=".", fixed=TRUE)
                alleles <- rep(0, length(templist))
                for(i in 1:length(alleles)){
                    alleles[i] <- templist[[i]][2]
                }
                # let's just not do null allelles for now
                # (Hey, I'm trying to graduate.)
                alleles <- sort(as.integer(alleles[alleles != "null"]))
                na1 <- length(alleles)
                # setting stuff up from De Silva method
                ng <- na1 # number of genotypes
                for(j in 2:m2){
                    ng <- ng*(na1+j-1)/j
                }
                ag <- .genlist(ng, na1, m2)
                temp <- .ranmul(ng, na1, ag, m2)
                rmul <- temp[[1]]
                arep <- temp[[2]]
                rm(temp)
                smat <- .selfmat(ng, na1, ag, m2)
                m <- m2/2
                smatdiv <- (.G(m-1,m+1))^2
                p1 <- rep(0, na1) # vector to hold frequencies
                for(a in alleles){
                    p1[match(a, alleles)] <- subfreq[1,
                                                     grep(a, names(subfreq),
                                                          fixed=TRUE)]
                }
                p1 <- p1/sum(p1) # normalize to sum to 1
                rvec <- rep(0,ng)
                for(g in 1:ng){
                    rvec[g] <- rmul[g]
                    for(j in 1:m2){
                        rvec[g] <- rvec[g]*p1[ag[g,j]]
                    }
                }
                id <- diag(nrow=ng)
                smatt <- smat/smatdiv
                s3 <- id - self * smatt
                s3inv <- solve(s3)
                gprob <- (1-self) * s3inv %*% rvec
                for(s in psamples){
                    gprobs[[s,L]] <- genotypeProbs(object, s, L, gprob=gprob,
                                                   alleles=alleles)
                }
            }
        }
    }

    # set up a 3D array to hold distances by sample x sample x locus
    loci.matrices<-array(dim=c(length(loci),length(samples),length(samples)),
                         dimnames=list(loci,samples,samples))

    # cycle through calculations
    for(L in loci){
        u <- Usatnts(object)[L]
        for(m in samples){
            for(n in samples[match(m,samples):length(samples)]){
                totdist <- 0 # total distance so far
                # cycle through all unambiguous genotypes for these two samples
                for(i in 1:length(gprobs[[m,L]][["probs"]])){
                    for(j in 1:length(gprobs[[n,L]][["probs"]])){
                        # raw distance
                        d <- distmetric(
                                 as.vector(gprobs[[m,L]][["genotypes"]][i,]),
                                 as.vector(gprobs[[n,L]][["genotypes"]][j,]),
                                        usatnt=u, missing=Missing(object), ...)
                        # probability of this combination
                        p <- gprobs[[m,L]][["probs"]][i] *
                             gprobs[[n,L]][["probs"]][j]
                        totdist <- totdist + (d*p)
                    }
                }
                loci.matrices[L,m,n] <- totdist
                loci.matrices[L,n,m] <- totdist
                if(progress) print(c(L, m, n))
            }
        }
    }

    # calculate the mean matrix across all loci
    mean.matrix <- meandist.from.array(loci.matrices)

    # return either the mean matrix and possibly the array as well
    if(all.distances){
        return(list(DistByLoc=loci.matrices, MeanMatrix=mean.matrix))
    } else {
        return(mean.matrix)
    }
}
# expect ambiguous genotypes to be a non-zero distance from themselves!

assignClones <- function(d, samples=dimnames(d)[[1]], threshold=0){
    # set up results vector
    results <- rep(NA, length(samples))
    names(results) <- samples

    results[1] <- 1 # the first individual is clone 1
    numclones <- 1 # there is one clone so far

     # for just one sample, we're done
    if(length(samples)>1){

    # assign the rest
    for(i in 2:length(samples)){
        s <- samples[i] # i is the numerical index, s is the character index
        sam <- samples[1:(i-1)] # samples that have already been assigned
        if(all(is.na(d[s,sam]))) next # skip missing data

        if(min(d[s, sam], na.rm=TRUE) > threshold){ # new clone
            numclones <- numclones + 1
            results[i] <- numclones
        } else { # match to existing clone(s)
            msam <- sam[d[s,sam] <= threshold] # samples that match
            cl <- unique(results[msam]) # the matching clone(s)
            cl <- cl[!is.na(cl)]
            if(length(cl) == 1){ # these are all one clone
                results[i] <- cl
            } else { # multiple clones are being merged due to this individual
                results[i] <- cl[1]
                results[results %in% cl[2:length(cl)]] <- cl[1]
            }
        }
    }

    # cleanup; just give sequential numbers
    clones <- unique(results)
    clones <- clones[!is.na(clones)]
    results <- match(results, clones)
    names(results) <- samples
  }

    return(results)
}

# Function for the genotype loss/addition model of the Bruvo measure.
# Figures out all the virtual alleles, then passes equal-length genotypes
# to Bruvo.distance.
# See Fig. 1b-1c and Eqn. 6 of Bruvo et al. 2004.
Bruvo2.distance <- function(genotype1, genotype2, maxl=7, usatnt=2, missing=-9,
                            add=TRUE, loss=TRUE){
    if(length(genotype1)==length(genotype2) ||
       genotype1[1] == missing ||
       genotype2[1] == missing ||
       (!add && !loss)){
        # Let normal Bruvo.distance return NA for missing genotypes
        # or pass directly to Bruvo.distance if there are no virtual alleles,
        # or if user isn't using genome loss or genome addition model.
        d <- Bruvo.distance(genotype1, genotype2, maxl=maxl, usatnt=usatnt,
                            missing=missing)
    } else {
        if(length(genotype1) > maxl || length(genotype2) > maxl){
            # if any has two many alleles, don't bother with enumeration
            d <- NA
        } else {
            # Assign the long and short genotypes
            if(length(genotype1) >= length(genotype2)) {
                genotypeL <- genotype1
                genotypeS <- genotype2
            } else {
                genotypeL <- genotype2
                genotypeS <- genotype1
            }
            # get the difference in length
            diff <- length(genotypeL) - length(genotypeS)
            # get allele combinations to try
            alcomb <- function(genotype){
                lg <- length(genotype)
                mat <- matrix(nrow=0, ncol=lg^diff)
                for(i in 1:diff){
                    mat <- rbind(mat, rep(genotype, times=lg^(i-1),
                                          each=lg^(diff-i)))
                }
                return(mat)
            }

            if(add){
                alleleadd <- alcomb(genotypeS)
            }
            if(loss){
                alleleloss <- alcomb(genotypeL)
            }
            # set up vectors to hold results
            distadd <- rep(NA, times=length(genotypeS)^diff)
            distloss <- rep(NA, times=length(genotypeL)^diff)
            # loop through calculations
            if(add){
            for(i in 1:length(distadd)){
                if(length(unique(alleleadd[,i])) > 1){
                    # check to see if this calculation has already been done
                    for(j in 1:(i-1)){
                        if(identical(sort(alleleadd[,i]),sort(alleleadd[,j]))){
                            distadd[i] <- distadd[j]
                            break
                        }
                    }
                }
                if(is.na(distadd[i])){
                    # do the calculation is this allele combo is new
                distadd[i] <- Bruvo.distance(genotypeL,
                                             c(genotypeS, alleleadd[,i]),
                                             maxl=maxl, usatnt=usatnt,
                                             missing=missing)}
            }}
            if(loss){
            for(i in 1:length(distloss)){
                if(length(unique(alleleloss[,i])) > 1){
                    for(j in 1:(i-1)){
                        if(identical(sort(alleleloss[,i]),
                                     sort(alleleloss[,j]))){
                            distloss[i] <- distloss[j]
                            break
                        }
                    }
                }
                if(is.na(distloss[i])){
                distloss[i] <- Bruvo.distance(genotypeL,
                                             c(genotypeS, alleleloss[,i]),
                                             maxl=maxl, usatnt=usatnt,
                                             missing=missing)}
            }}
            # get average
            d <- mean(c(mean(distadd), mean(distloss)), na.rm=TRUE)
        }
    }
    return(d)
}
