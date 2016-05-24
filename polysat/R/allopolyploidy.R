# "Resolving microsatellite genotype ambiguity in populations of
# allopolyploid and diploidized autopolyploid organisms using negative
# correlations between alleles"

# function to create simulated data
simAllopoly <- function(ploidy=c(2,2),
                         n.alleles=c(4,4), n.homoplasy=0,
                         n.null.alleles=rep(0, length(ploidy)),
                         alleles=NULL,
                         freq=NULL, meiotic.error.rate=0,
                         nSam=100, locname="L1"){
    # number of different genomes in the allopolyploid (1 for autopolyploid)
    nGenomes <- length(ploidy)
    # no meiotic error if there is only one genome
    if(nGenomes==1 && meiotic.error.rate != 0)
      stop("No meiotic error possible if there is only one subgenome.")

    # make allele names if needed
    if(is.null(alleles)){
        if(length(n.alleles)!=nGenomes)
              stop("'n.alleles' must be a vector of the same length as 'ploidy'.")
          if(nGenomes>6) stop("More than six subgenomes not supported at this time.")
          GenomeNames <- c('A','B','C','D','E','F')[1:nGenomes]
          alleles <- list()
          length(alleles) <- nGenomes
          for(i in 1:nGenomes){
              alleles[[i]] <- paste(GenomeNames[i], 1:n.alleles[i], sep="-")
          }
          if(length(n.homoplasy)!=1) stop("Only one value allowed for n.homoplasy.")
          if(n.homoplasy>0){ # if there are more than 0 homoplasious alleles
              HAlleles <- paste('H', 1:n.homoplasy, sep="-")
              for(i in 1:nGenomes){ # put homoplasious alleles at end of vector
                  alleles[[i]][(n.alleles[i]-n.homoplasy+1):n.alleles[i]] <- HAlleles
              }
          }
          # if there are null alleles, put at beginning of vector
          if(length(n.null.alleles)!=nGenomes)
              stop("'n.null.alleles' must be a vector the same length as 'ploidy'.")
          if(!all(n.null.alleles<=n.alleles))
              stop("'n.null.alleles' values must be smaller than 'n.alleles' values.")
          for(i in 1:nGenomes){
              if(n.null.alleles[i]>0){
                  alleles[[i]][1:n.null.alleles[i]] <- "N"
              }
          }
    } else {
        # if alleles are provided, get n.alleles to match
        n.alleles <- sapply(alleles, length)
    }

    # randomly generate allele frequencies if needed
    if(is.null(freq)){
        freq <- list()
        length(freq) <- nGenomes
        for(i in 1:nGenomes){
            temp <- sample(100, n.alleles[i], replace=TRUE)
            freq[[i]] <- temp/sum(temp)
        }
    }
    
    # errors
    if(!is.list(alleles) || !is.list(freq))
        stop("alleles and freq must each be a list of vectors")
    if(!identical(sapply(freq, length), sapply(alleles, length)))
        stop("Need same number of values for alleles and frequencies.")
    if(nGenomes != length(alleles) || nGenomes != length(freq))
      stop("Number of genomes nees to match between ploidy, alleles, and freq")
    if(length(meiotic.error.rate)!=1 || meiotic.error.rate<0 ||
       meiotic.error.rate>0.5)
        stop("meiotic.error.rate should be a single value between 0 and 0.5.")
    if(meiotic.error.rate>0 && !all(ploidy %% 2 == 0))
        stop("Even (not odd) ploidy is required if simulating meiotic error.")

    # create genambig object to return at end
    gen <- new("genambig", samples=1:nSam, loci=locname)

    # fixing the sample function (Thanks, help file on "sample"!)
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    # create genotypes
    for(s in 1:nSam){
        # determine if there was a meiotic error
        gameteploidy <- ploidy/2
        if(sample(c(TRUE,FALSE),size=1,
                  prob=c(meiotic.error.rate, 1-meiotic.error.rate))){
            isolocusLost <- sample(nGenomes,1)
            isolocusGained <- resample((1:nGenomes)[-isolocusLost],1)
            gamete1ploidy <- gameteploidy
            gamete1ploidy[isolocusLost] <-  gamete1ploidy[isolocusLost]-1
            gamete1ploidy[isolocusGained] <- gamete1ploidy[isolocusGained]+1
        } else {
            gamete1ploidy <- gameteploidy
        }
        if(sample(c(TRUE,FALSE),size=1,
                  prob=c(meiotic.error.rate, 1-meiotic.error.rate))){
            isolocusLost <- sample(nGenomes,1)
            isolocusGained <- resample((1:nGenomes)[-isolocusLost],1)
            gamete2ploidy <- gameteploidy
            gamete2ploidy[isolocusLost] <-  gamete2ploidy[isolocusLost]-1
            gamete2ploidy[isolocusGained] <- gamete2ploidy[isolocusGained]+1
        } else {
            gamete2ploidy <- gameteploidy
        }
        tempploidy <- gamete1ploidy + gamete2ploidy
        # get genotype
        thisgen <- unique(unlist(mapply(resample, alleles,
                                                       replace= TRUE,
                                      size = tempploidy, prob = freq,
                                        SIMPLIFY=FALSE)))
        thisgen <- thisgen[thisgen!=0 & thisgen!="N"] # remove null alleles
        if(length(thisgen)==0){
            Genotype(gen,s,1) <- Missing(gen)
        } else {
            Genotype(gen,s,1) <- thisgen
        }
    }

    return(gen)
}

# n.start is passed to the kmeans function
alleleCorrelations <- function(object, samples=Samples(object), locus=1,
                         alpha=0.05, n.subgen=2, n.start=50){
    # subset by samples if needed
    object <- object[samples,]
    # check missing data
    mymissing <- isMissing(object, loci=locus)
    if(all(mymissing)) stop("All data at locus are missing.")
    # only do one locus at a time
    if(length(locus) > 1) stop("More than one locus")
    # subset dataset object by this locus
    object <- object[!mymissing,locus]
    # perform class conversion if necessary
    if(is(object, "genambig")) object <- genambig.to.genbinary(object)
    # error if not the right class
    if(!is(object, "genbinary")) stop("genambig or genbinary object needed.")
    # set Absent to 0 and Present to 1
    Absent(object) <- as.integer(0)
    Present(object) <- as.integer(1)

    gentable <- Genotypes(object)  # get the 1/0 table of genotypes
    # remove alleles that aren't in dataset
    gentable <- gentable[,apply(gentable, 2,
                                function(x) !all(x==Absent(object)))]
    gentable.out <- gentable # to be output
    # find any non-variable alleles (present in all genotypes)
    nonvar <- apply(gentable, 2, function(x) all(x==Present(object)))
    # are there too many non-variable alleles?
    mingen <- sum(nonvar)
    if(mingen > n.subgen){
        stop(paste("Locus",locus,": too many non-variable alleles"))
    }
    nAl <- dim(gentable)[2]        # number of alleles
    # get allele names for output
    mysplit <- strsplit(dimnames(gentable)[[2]], split=".", fixed=TRUE)
    alleles <- sapply(mysplit, function(x) x[2])
    # set up matrix for clusters
    clustmat1 <- matrix(0, nrow=n.subgen, ncol=nAl,
                        dimnames=list(NULL,alleles))

    # if there are non-variable alleles, assign each to its own genome
    if(sum(nonvar) > 0){
        for(i in 1:sum(nonvar)){
            clustmat1[i,alleles[nonvar][i]] <- 1
        }
        sigmatNeg <- sigmatPos <- oddsRatio <- NULL
        pValuesNeg <- pValuesPos <- mydist <- NULL
        mytotss <- mybss <- NA
        clustmethod <- "fixed alleles" # gets overwritten later if necessary
        gentable <- gentable[,!nonvar]
        nAl <- nAl - sum(nonvar) # adjust number of alleles to only variable ones.
        alleles <- alleles[!nonvar]

        # if each isolocus has a non-variable allele, but there are
        # variable alleles left, make them all homoplasious
        if(sum(nonvar)==n.subgen){
            clustmat1[,alleles] <- 1
        }
    }
    if(sum(rowSums(clustmat1)==0)==1){
        # if everything else can be assigned to one genome, do that
        clustmat1[rowSums(clustmat1)==0,colSums(clustmat1)==0] <- 1
        sigmatNeg <- sigmatPos <- oddsRatio <- NULL
        pValuesNeg <- pValuesPos <- mydist <- NULL
        mytotss <- mybss <- NA
        clustmethod <- "fixed alleles"
    }
    # make a second copy of assignment matrix for UPGMA
    clustmat2 <- clustmat1
    
    ### how the function will normally work (all alleles variable): ###
    if(sum(rowSums(clustmat1)==0)>1){

        clustmethod <- "K-means and UPGMA"

    # set up matrices for results of Fisher's exact test
    pValuesNeg <- matrix(NA, nrow=nAl, ncol=nAl) # significance of negative corr.
    pValuesPos <- matrix(NA, nrow=nAl, ncol=nAl) #significance of positive corr.
    oddsRatio <- matrix(NA, nrow=nAl, ncol=nAl) # (00*11)/(01*10)

    # loop to run Fisher's exact test on each pair of alleles
    for(a in 1:(nAl-1)){
        for(b in (a+1):nAl){
            X <- fisher.test(gentable[,a], gentable[,b], alternative="less",
                             conf.int=FALSE)
            pValuesNeg[a,b] <- X$p.value; pValuesNeg[b,a] <- X$p.value
            Y <- fisher.test(gentable[,a], gentable[,b], alternative="greater",
                             conf.int=FALSE)
            pValuesPos[a,b] <- Y$p.value; pValuesPos[b,a] <- Y$p.value
            oddsRatio[a,b] <- oddsRatio[b,a] <- X$estimate
        }
    }

    # Holm-Bonferroni correction to get p-value threshold
    # Negative correlations (for assigning alleles)
    pVect <- sort(pValuesNeg[lower.tri(pValuesNeg)]) # smallest to largest pvalues
    m <- length(pVect) # number of hypothesis tests
    kTest <- pVect > alpha/(m + 1 - (1:m))
    if(all(!kTest)){ # reject all hypotheses
        sigmatNeg <- matrix(TRUE, nrow=nAl, ncol=nAl)
    } else {
        k <- min((1:m)[kTest])
        maxP <- pVect[k-1] # maximum p value for significance
        if(k > 1){
            sigmatNeg <- pValuesNeg <= maxP
        } else {
            sigmatNeg <- matrix(FALSE, nrow=nAl, ncol=nAl)
        }
    }
    # Positive correlations (for issuing warnings)
    pVect <- sort(pValuesPos[lower.tri(pValuesPos)]) # smallest to largest pvalues
    kTest <- pVect > alpha/(m + 1 - (1:m))
    if(all(!kTest)){ # reject all hypotheses
        sigmatPos <- matrix(TRUE, nrow=nAl, ncol=nAl)
    } else {
        k <- min((1:m)[kTest])
        maxP <- pVect[k-1] # maximum p value for significance
        if(k > 1){
            sigmatPos <- pValuesPos <= maxP
        } else {
            sigmatPos <- matrix(FALSE, nrow=nAl, ncol=nAl)
        }
    }
    if(!all(!sigmatPos)){
        cat(paste("Warning: Significant positive correlations between alleles at locus", Loci(object), "; population structure or scoring error may bias results."), sep="\n")
    }

    dimnames(sigmatNeg) <- list(alleles,alleles)
    dimnames(sigmatPos) <- list(alleles,alleles)
    dimnames(pValuesNeg) <- list(alleles, alleles)
    dimnames(pValuesPos) <- list(alleles, alleles)
    dimnames(oddsRatio) <- list(alleles, alleles)

    # replace NA with 0 in the p-value matrix
    mydist <- pValuesNeg
    mydist[is.na(mydist)] <- 0
    # adjust number of clusters if there were non-variable alleles
    n.subgen.adj = n.subgen - sum(rowSums(clustmat1)>0)
    # kmeans clustering of alleles based on pvalues
    kres <- kmeans(mydist, centers=n.subgen.adj, nstart=n.start)
    myclust1 <- kres$cluster
        mytotss <- kres$totss
        mybss <- kres$betweenss
        
    # UPGMA clustering    
        myclust2 <- cutree(hclust(as.dist(mydist), method="average"),
                           k = n.subgen.adj)
    
    # put clusters into matrix
    for(i in (n.subgen - n.subgen.adj + 1):n.subgen){
        clustmat1[i,alleles[myclust1==(i - n.subgen + n.subgen.adj)]] <- 1
        clustmat2[i,alleles[myclust2==(i - n.subgen + n.subgen.adj)]] <- 1
    } } # end of if clause for if there aren't non-variable alleles

    return(list(locus=Loci(object), clustering.method=clustmethod,
        significant.neg=sigmatNeg, significant.pos=sigmatPos,
                p.values.neg=pValuesNeg, p.values.pos=pValuesPos,
                odds.ratio=oddsRatio,
                Kmeans.groups=clustmat1, UPGMA.groups=clustmat2,
                heatmap.dist=mydist,
                totss=mytotss, betweenss=mybss, gentable=gentable.out))
}

# function to test consistency of genotypes with allele groups
# will also make adjustment to the groups
# genobject is genambig or genbinary
# fisherResults is output of AlleleCorrelations
# SGploidy is the ploidy of each subgenome, e.g. 2 for an allotetraploid
# null.weight should be 0 if there are definitely null alleles, or 1 if
# there are definitely not.
# tolerance is the proportion of genotypes that can be allowed to be inconsistent
# with the allele assignments.
# rare.al.check is a threshold for the proportion of genotypes in which an allele
# should be present to be considered rare.  Rare alleles are checked to see
# if they should be switched to a different isolocus.
testAlGroups <- function(object, fisherResults, SGploidy=2,
                         samples=Samples(object),
                         null.weight=0.5, tolerance=0.05,
                         rare.al.check=0.2){
    # make sure genotype data is right class
    if(!is(object, "genambig") && !is(object, "genbinary"))
        stop("genambig or genbinary object needed.")
    object <- object[samples,]
    # locus
    L <- fisherResults$locus
    # remove missing data
    mymissing <- isMissing(object, loci=L)
    # number of individuals
    nind <- sum(!mymissing)
    # subset genotypes
    genobject <- object[!mymissing,L]
    # change class if necessary
    if(is(genobject, "genbinary")) genobject <- genbinary.to.genambig(genobject)
    # alleles
    alleles <- dimnames(fisherResults$Kmeans.groups)[[2]]

    # number of subgenomes
    n.subgen <- dim(fisherResults$Kmeans.groups)[1]

    # get the number of alleles in each genotype
    numAl <- sapply(Genotypes(genobject), length)
    # check to see that none have too many alleles
    if(!all(numAl <= n.subgen*SGploidy))
        stop(paste("One or more genotypes have more than",
                   n.subgen*SGploidy, "alleles."))
    # check to see that alleles match in genobject and fisherResults
    if(!all(as.character(.unal1loc(genobject, samples=Samples(genobject),
                        locus=L)) %in% alleles))
        stop("Alleles found in genobject that aren't in fisherResults.")

    # Compare K-means and UPGMA results
    # Sort them so groups are the same
    Kmat <- fisherResults$Kmeans.groups
    Umat <- fisherResults$UPGMA.groups
    KsortV <- UsortV <- integer(length(alleles))
    Uswap <- Kswap <- integer(n.subgen) # index = new number, number=row in matrix
    for(i in 1:n.subgen){
        KfirstLeft <- match(0, KsortV) # first unassigned allele
        UfirstLeft <- match(0, UsortV)
        # match rows to new subgenomes
        Kswap[i] <- match(1, Kmat[,KfirstLeft])
        Uswap[i] <- match(1, Umat[,UfirstLeft])
        # assign alleles to new subgenomes
        KsortV[Kmat[Kswap[i],] == 1] <- i
        UsortV[Umat[Uswap[i],] == 1] <- i
    }
    # were the results identical?
    if(identical(KsortV, UsortV)){
        G <- Kmat # use K-means assignments (same as UPGMA assignments)
    } else {
        # if K-means and UPGMA gave different results, see which is more consistent
        # with genotypes.
        Kok <- 0 # count consistent with K-means assignments
        Uok <- 0 # count consistent with UPGMA assignments
        for(s in Samples(genobject)){
            # get genotype and assignment matrices for this individual
            gen <- as.character(Genotype(genobject, s, L))
            Ksubg <- Kmat[,gen, drop=FALSE]
            Usubg <- Umat[,gen, drop=FALSE]
            # check whether assignments are okay
            if(all(rowSums(Ksubg)>0) &
               all(rowSums(Ksubg[,colSums(Ksubg)==1, drop=FALSE]) <= SGploidy)){
                Kok <- Kok + 1
            }
            if(all(rowSums(Usubg)>0) &
               all(rowSums(Usubg[,colSums(Usubg)==1, drop=FALSE]) <= SGploidy)){
                Uok <- Uok + 1
            }
        }
        if(Uok > Kok){
            G <- Umat # UPGMA was better, use this matrix
        } else {
            G <- Kmat # Use K-means if it was better or if there was a tie.
        }
    }

    # First, check if rare alleles should be swapped to a different isolocus.
    # determine which alleles are rare
    AlOcc <- colSums(fisherResults$gentable)/dim(fisherResults$gentable)[1]
    AlToCheck <- alleles[AlOcc <= rare.al.check]
    continueChecking <- TRUE
    while(continueChecking){
        continueChecking <- FALSE
    for(a in AlToCheck){
        # set up alternative assignment matrices
        G.alts <- list()
        length(G.alts) <- n.subgen
        for(i in 1:n.subgen){
            G.temp <- G
            G.temp[i,a] <- 1
            G.temp[-i,a] <- 0
            G.alts[[i]] <- G.temp
        }
        # vector to count genotypes consistent with alternative matrices
        ok.count <- rep(0, n.subgen)
        # check genotypes against alternative matrices
        for(s in Samples(genobject)){
            gen <- as.character(Genotype(genobject, s, L)) # the genotype
            if(!a %in% gen) next # skip if it doesn't have this allele
            for(i in 1:n.subgen){
                subg <- G.alts[[i]][,gen, drop=FALSE]
                if(all(rowSums(subg)>0) &
                   all(rowSums(subg[,colSums(subg)==1, drop=FALSE]) <= SGploidy)){
                    ok.count[i] <- ok.count[i] + 1
                }
            }
        }
        # only change G if there is not a tie
        if(sum(ok.count==max(ok.count))==1 &&
           !identical(G, G.alts[[which(ok.count==max(ok.count))]])){
            G <- G.alts[[which(ok.count==max(ok.count))]]
            continueChecking <- TRUE # go through all remaining alleles again
            AlToCheck <- AlToCheck[AlToCheck!=a]
        }
    }}

    # set up internal function to check consistency for adding homoplasy
    makeA <- function(){
        A <- matrix(0, nrow=n.subgen, ncol=length(alleles),
                dimnames=list(NULL,alleles))
        num.error <- 0 # number of genotypes with scoring or meiotic error
        for(s in Samples(genobject)){
            gen <- as.character(Genotype(genobject, s, L)) # the genotype
 
            # get 1/0 matrix for presence of these alleles in subgenomes
            subG <- G[,gen, drop=FALSE]
            # Start checking the genotype.
            # It is okay if no genome has more alleles than is possible AND
            # no genome has zero alleles.
            # Otherwise, figure out which alleles might need to be copied:
            if(!all(rowSums(subG[,colSums(subG) == 1, drop=FALSE])
                    <= SGploidy)){
                # If one subgenome has more alleles than are possible
                donor <-
                (1:n.subgen)[rowSums(subG[,colSums(subG) == 1, drop=FALSE]) >
                             SGploidy]
                # they might be copied to a subgenome(s) that has fewer than the
                # maximum.
                recipient <- (1:n.subgen)[rowSums(subG) < SGploidy]
                if(length(recipient)==0){
                    recipient <- (1:n.subgen)[rowSums(subG[,colSums(subG) == 1,
                                                   drop=FALSE]) < SGploidy]
                }
                gendonor <- gen[colSums(subG[recipient,,drop=FALSE]==0)>0 &
                                colSums(subG[donor,,drop=FALSE]==1)>0]
                A[recipient,gendonor] <- A[recipient,gendonor] + 1

                num.error <- num.error + 1
            } else {
                # If one subgenome has no alleles
                if(!all(rowSums(subG) > 0) && null.weight > 0){
                    recipient <- (1:n.subgen)[rowSums(subG)==0]
                    A[recipient,gen] <- A[recipient,gen] + null.weight
                    num.error <- num.error + 1
                }
            }
        } # end of Samples loop
        
        # calculate the proportion inconsistent with allele assignments
        prop.error <- num.error/nind

        # Make sure not to copy alleles to where they have already been
        # assigned (problem in hexaploids).
        A[G==1] <- 0

        # weight the A matrix by allele occurance
        # omitted; is only helpful for rare homoplasious alleles.
 #       A <- sweep(A, 2, AlOcc, "/")

        return(list(A, prop.error))
    } # end of makeA function

    # check consistency once before starting assingnment and re-checking loop
    temp <- makeA()
    A <- temp[[1]]
    prop.error <- temp[[2]]

    # loop to copy alleles to other genomes, then re-check
    numloops <- 0
    while(prop.error > tolerance && numloops < 1000){
        numloops <- numloops + 1 # precaution to keep it from getting stuck

        A <- A[,AlOcc != 1] # don't copy fixed alleles

        # choose best allele to make homoplasious for this round
        tocopy <- which(A == max(A), arr.ind=TRUE)
        if(dim(tocopy)[1]==1){ # if one is clearly the best
            thisallele <- dimnames(A)[[2]][tocopy[,"col"]]
            thissubgen <- tocopy[,"row"]
        } else { # choose among several with the highest scores
            if(!is.null(fisherResults$p.values.neg)){ # use p values
                meanp <- numeric(dim(tocopy)[1])
                for(i in 1:dim(tocopy)[1]){
                    # alleles already in this genome
                    thisgenal <- alleles[G[tocopy[i,"row"],] == 1]
                    # remove any fixed alleles
                    thisgenal <- thisgenal[!thisgenal %in% alleles[AlOcc == 1]]
                    # get mean p values between this allele and alleles already in genome
                    thisal <- dimnames(A)[[2]][tocopy[i,"col"]]
                    meanp <- mean(fisherResults$p.values.neg[thisgenal,thisal])
                }
                bestp <- which(meanp == min(meanp))
                if(length(bestp) > 1){ # if there's a tie, choose at random
                    bestp <- sample(bestp, size=1)
                }
                thisallele <- dimnames(A)[[2]][tocopy[bestp,"col"]]
                thissubgen <- tocopy[bestp,"row"]
            } else { # pick at random if there are no p-values to check
                myrand <- sample(dim(tocopy)[1], size=1)
                thisallele <- dimnames(A)[[2]][tocopy[myrand,"col"]]
                thissubgen <- tocopy[myrand,"row"]
            }
        }
        G[thissubgen, thisallele] <- 1

        # check genotypes against new G matrix
        temp <- makeA()
        A <- temp[[1]]
        prop.error <- temp[[2]]
    } # end of while loop
    return(list(locus=L, SGploidy=SGploidy, assignments=G))
}

# function to assign alleles to subgenomes based on the procedure of
# Catalan et al. 2006 (Genetics 172:1939-1953).
# If any genotypes have fewer alleles than n.subgen, there must be homoplasy,
# so exclude locus. (In their paper, <2 alleles for allotetraploid.)
# Any genotypes with as many alleles as n.subgen (2 for allotetraploid) must be
# homozygous at each isolocus, so these are used first to assign alleles to
# isoloci.
# Check that this assignment is consistent with the genotypes that have a
# number of alleles > n.subgen, or if there is homoplasy.
# This is not stated explicitly in paper, but genotypes with n.subgen*SGploidy
# alleles could be used to further resolve ambiguity if genotypes with
# n.subgen alleles are not sufficient.
catalanAlleles <- function(object, samples=Samples(object), locus=1,
                           n.subgen=2, SGploidy=2, verbose=FALSE){
    # subset object if needed
    object <- object[samples,]
    # check missing data
    mymissing <- isMissing(object, loci=locus)
    if(all(mymissing)) stop("All data at locus are missing.")
    # only do one locus at a time
    if(length(locus) > 1) stop("More than one locus")
    # subset dataset object by this locus
    object <- object[!mymissing,locus]
    # perform class conversion if necessary
    if(is(object, "genbinary")) object <- genbinary.to.genambig(object)
    # error if not the right class
    if(!is(object, "genambig")) stop("genambig or genbinary object needed.")
    # error if just one subgenome
    if(n.subgen<2) stop("Need at least two subgenomes.")

    # get the number of alleles in each genotype
    numAl <- sapply(Genotypes(object), length)
    # check to see that none have too many alleles
    if(!all(numAl <= n.subgen*SGploidy))
        stop(paste("One or more genotypes have more than",
                   n.subgen*SGploidy, "alleles."))

    # check to see if any genotypes have too few alleles, implying homoplasy
    if(!all(numAl >= n.subgen)){
        result <- list(locus=Loci(object), SGploidy=SGploidy,
    assignments="Homoplasy or null alleles: some genotypes have too few alleles")
    } else {
        # set up results matrix 
        alleles <- .unal1loc(object, Samples(object), locus)
        G <- matrix(0, nrow=n.subgen, ncol=length(alleles),
                    dimnames=list(NULL, as.character(alleles)))
        # start with genotypes that have n.subgen alleles
        X <- Genotypes(object)[numAl==n.subgen,]
        if(length(X)==0){  # quit if there aren't any
            result <- list(locus=Loci(object), SGploidy=SGploidy,
                assignments="Unresolvable: no genotypes with n.subgen alleles.")
        } else {
            X <- unique(X) # get unique genotypes w/ n. subgen alleles
            # vector to rank genotypes by allele frequencies
            AlFreqInX <- sapply(alleles,
                                function(a) sum(sapply(X,
                                                       function(x) a %in% x)))
            names(AlFreqInX) <- as.character(alleles)
            Xrank <- sapply(X, function(x) sum(AlFreqInX[as.character(x)]))
            # sort genotypes, putting the most interconnected first
            X <- X[order(Xrank, decreasing=TRUE)]
            
            # use the first one
            gen <- as.character(X[[1]])
            for(i in 1:n.subgen){
                G[i,gen[i]] <- 1
            }
            ToAssign <- alleles[colSums(G)==0] # alleles to assign
            Xa <- X[sapply(X, function(x) sum(ToAssign %in% x)==1)]
            
            # try using the rest of the genotypes with n.subgen alleles
             while(length(Xa)>0){
                 gen <- as.character(Xa[[1]])
                    g <- G[,gen]
                    thissubgen <- rowSums(g)==0
                    thisallele <- gen[colSums(g)==0]
                    G[thissubgen,thisallele] <- 1
                    ToAssign <- alleles[colSums(G)==0] # alleles to assign
                    Xa <- X[sapply(X, function(x) sum(ToAssign %in% x)==1)]
             }

            # if needed, use heterozygous genotypes to further resolve.
            if(!all(colSums(G)!=0)){
                # get genotypes that will be informative
                X <- Genotypes(object)[,1]
                X <- unique(X)
                     # one allele to assign, one subgenome has no alleles yet
                Xa <- X[sapply(X, function(x) (sum(ToAssign %in% x)==1 &&
                        sum(rowSums(G[,as.character(x)])==0)==1) ||
                     # OR, all subgenomes but one have max number of alleles
                        (sum(ToAssign %in% x) > 0 &&
                    sum(rowSums(G[,as.character(x)])==SGploidy)==n.subgen-1))]
                while(length(Xa)>0){
                    # use the first genotype in the list to assign allele(s)
                    g <- G[,as.character(Xa[[1]])]
                    thissubgen <- rowSums(g) == min(rowSums(g))
                    thisallele <- dimnames(g)[[2]][colSums(g)==0]
                    G[thissubgen,thisallele] <- 1
                    
                    # update list of unassigned alleles and useful genotypes
                    ToAssign <- alleles[colSums(G)==0]
                    Xa <- X[sapply(X, function(x) (sum(ToAssign %in% x)==1 &&
                        sum(rowSums(G[,as.character(x)])==0)==1) ||
                        (sum(ToAssign %in% x) > 0 &&
                    sum(rowSums(G[,as.character(x)])==SGploidy)==n.subgen-1))]
                }
            }

            # check results
            if(!all(colSums(G)!=0)){
                result <- list(locus=Loci(object), SGploidy=SGploidy,
                                   assignments="Unresolvable")
                if(verbose) print(G)
            } else {
                # loop through each genotype to check consistency with results
                badgen <- list()
                bgIndex <- 1
                for(s in Samples(object)){
                    g <- G[,as.character(Genotype(object, s, locus))]
                    gr <- rowSums(g)
                    if(!all(gr > 0 & gr <= SGploidy)){
                        badgen[[bgIndex]] <- Genotype(object,s,locus)
                        bgIndex <- bgIndex+1
                    }
                }
                if(length(badgen)==0){
                    result <- list(locus=Loci(object), SGploidy=SGploidy,
                                   assignments=G)
                } else {
                    if(verbose){
                        cat("Allele assignments:", sep="\n")
                        print(G)
                        cat("Inconsistent genotypes:", sep="\n")
                        print(badgen)
                    }
                    result <- list(locus=Loci(object), SGploidy=SGploidy,
                                   assignments="Homoplasy or null alleles")
                }
            }
        }
    }

    if(verbose) print(result)
    return(result)
}

# function to merge results from the same locus from multiple populations.
# x is a list of lists; each element of x is a result such as that
# produced by catalanAlleles or testAlGroups.
mergeAlleleAssignments <- function(x){
    loci <- sapply(x, function(z) z$locus)   # extract locus names
    SGp <- sapply(x, function(z) z$SGploidy) # extract ploidies
    lociU <- unique(loci)
    if(is.list(loci) || is.list(SGp) || !all(!is.na(lociU)) || !all(!is.na(SGp)))
        stop("Locus and SGploidy value required for each list element.")
    # set up results object
    results <- list()
    length(results) <- length(lociU)

    # function to help with merging
    onezero <- function(y,z){
        res <- y+z
        res[res==2] <- 1
        return(res)
    }

    # loop through each locus
    for(L in 1:length(lociU)){
        locus <- lociU[L]
        thisSGp <- unique(SGp[loci==locus]) # get ploidy for this locus
        if(length(thisSGp)>1) stop("More than one SGploidy value per locus.")
        # get allele assignments
        A <- lapply(x[loci==locus], function(z) z$assignments)
        # throw out any that are not a matrix
        A <- A[sapply(A, function(z) is.matrix(z))]
        if(length(A)==0){
            G <- "No assignment"
        } else {
        # merge assignment matrices, or just keep the first if there is only one
            # sort A so that matrices with the most alleles are used first
            A <- A[order(-sapply(A, function(x) dim(x)[2]),1:length(A))]
            # get unique alleles
            alleles <- unique(unlist(lapply(A, function(z) dimnames(z)[[2]])))
            # set up matrix
            G <- matrix(0, nrow=dim(A[[1]])[1], ncol=length(alleles),
                        dimnames=list(NULL, alleles))
            # fill in from first matrix
            G[,dimnames(A[[1]])[[2]]] <- A[[1]]
            A <- A[-1]
            # remaining matrices
            stuckcount <- 0
            while(length(A)>0){
                if(stuckcount==length(A)){
                    G <- "No assignment"
                    break
                }
                stuckcount <- 0  # to prevent endless looping
                toremove <- integer(0)
                for(i in 1:length(A)){
                    a <- A[[i]]
                    assigned <- alleles[colSums(G)>0]
                    # find alleles present in current matrix that are already assigned
                    overlap <- assigned[assigned %in% dimnames(a)[[2]][colSums(a)>0]]
                    if(length(overlap)==0){
                        stuckcount <- stuckcount+1
                        next
                    }
                    # try to match up rows
                    groups <- kmeans(rbind(G[,overlap],a[,overlap]),dim(G)[1])
                    # check if matching makes sense
                    if(!all(1:dim(G)[1] %in% groups$cluster[1:dim(G)[1]]) ||
                       groups$betweenss/groups$totss <= 0.5){
                        stuckcount <- stuckcount + 1
                        next
                    }
                    # combine matrices
                    G[groups$cluster[1:dim(G)[1]],dimnames(a)[[2]]] <-
                        onezero(G[groups$cluster[1:dim(G)[1]],dimnames(a)[[2]]],
                                a[groups$cluster[-(1:dim(G)[1])],])
                    toremove <- c(toremove,i) # done with this matrix
                }
                # get rid of matrices that have already been used
                if(length(toremove)>0) A <- A[-toremove] 
            }
        }
        results[[L]] <- list(locus=locus, SGploidy=thisSGp,
                             assignments=G)
    }
    return(results)
}

# function to rewrite genambig object using allele assignments
# object = genambig or genbinary object
# x = list of results (e.g. list of catalanAlleles outputs)
# allowAneuploidy = let individuals be aneuploid vs. inserting missing data
recodeAllopoly <- function(object, x, allowAneuploidy=TRUE,
                           samples=Samples(object), loci=Loci(object)){
    # data conversion and error checking
    if(is(object, "genbinary")) object <- genbinary.to.genambig(object)
    if(!is(object, "genambig"))
        stop("'genbinary' or 'genambig' object required.")
    object <- object[samples, loci]
    lociA <- sapply(x, function(z) z$locus)   # extract locus names
    SGp <- sapply(x, function(z) z$SGploidy) # extract ploidies
    if(is.list(lociA) || is.list(SGp) || !all(!is.na(lociA)) ||
       !all(!is.na(SGp)))
        stop("Locus and SGploidy value required for each list element in x.")
    if(!all(lociA %in% loci)){
        warning(paste("Loci", paste(lociA[!lociA %in% loci], collapse=TRUE),
                      "found in 'x' but not 'loci'."))
        x <- x[lociA %in% loci]
        if(length(x)==0)
            stop("No recoding possible because locus names do not match")
        lociA <- sapply(x, function(z) z$locus)   # extract locus names
        SGp <- sapply(x, function(z) z$SGploidy) # extract ploidies
    }
    # merge assignments if any loci are represented more than once
    if(!all(table(lociA)==1)){
        x <- mergeAlleleAssignments(x)
        lociA <- sapply(x, function(z) z$locus)   # extract locus names
        SGp <- sapply(x, function(z) z$SGploidy) # extract ploidies
    }
    # finish getting info by locus set up.
    # lociA = loci to recode; SGp = ploidy of isoloci, by locus;
    # nSG = number of isoloci, by locus; A = allele assignments, by locus.
    names(SGp) <- lociA
    nSG <- sapply(x, function(z) ifelse(is.matrix(z$assignments),
                                        dim(z$assignments)[1], NA))
    names(nSG) <- lociA
    A <- lapply(x, function(z) z$assignments)
    names(A) <- lociA
    # get rid of loci with no assignments
    lociA <- lociA[!is.na(nSG)]
    SGp <- SGp[lociA]
    nSG <- nSG[lociA]
    A <- A[lociA]

    # set up new object
    # duplicate the loci that should be duplicated
    newloci <- paste(rep(lociA, times=nSG),unlist(lapply(nSG, function(y) 1:y)),
                     sep="-")
    extraloci <- loci[!loci %in% lociA]
    newloci <- c(newloci, extraloci)
    result <- new("genambig", samples, newloci)
    PopNames(result) <- PopNames(object)
    PopInfo(result) <- PopInfo(object)[samples]
    Description(result) <- Description(object)
    Missing(result) <- Missing(object)
    newploidies <- matrix(rep(SGp, times=nSG), nrow=length(samples),
                             ncol=sum(!newloci %in% extraloci),
                             byrow=TRUE)
    # transfer ploidy data for loci that aren't being recoded
    if(length(extraloci) > 0){ 
        if(is(object@Ploidies, "ploidymatrix")){
            newploidies <- cbind(newploidies,
                                 Ploidies(object, samples=samples,
                                          loci=extraloci))
        }
        if(is(object@Ploidies, "ploidylocus")){
            newploidies <- cbind(newploidies,
                                 matrix(Ploidies(object, loci=extraloci),
                                        nrow=length(samples), ncol=length(extraloci),
                                        byrow=TRUE))
        }
        if(is(object@Ploidies, "ploidysample") || is(object@Ploidies, "ploidyone")){
            newploidies <- cbind(newploidies,
                                 matrix(Ploidies(object, loci=extraloci),
                                        nrow=length(samples), ncol=length(extraloci),
                                        byrow=FALSE))
        }
    }
    Ploidies(result) <- newploidies

    # collapse ploidy if we don't need it as a matrix
    if(!allowAneuploidy && plCollapse(result@Ploidies, na.rm=FALSE, returnvalue=FALSE)){
      result <- reformatPloidies(result, output="collapse")
    }

    # start assignments
    for(L in loci){
        # get positions of isoloci in new object
        Lindex <- which(sapply(strsplit(newloci, "-"), function(y) y[1]) == L)
        # copy over microsat repeat lengths
        Usatnts(result)[Lindex] <- Usatnts(object)[L]
        if(length(Lindex)==1){ # for loci with no assignments:
            Genotypes(result, loci=Lindex) <- Genotypes(object, loci=L)
            Ploidies(result)[,Lindex] <- Ploidies(object, Samples(object), L)
        } else {               # for loci with assignments:
            # test if this locus is numeric or integers, to convert back later
            xsamples <- samples[!isMissing(object, samples, L)]
            LisInt <- all(sapply(Genotypes(object, xsamples, L),
                                 is.integer, USE.NAMES=FALSE))
            LisNum <- all(sapply(Genotypes(object, xsamples, L),
                                 is.numeric, USE.NAMES=FALSE))
            # If there are any alleles that are unassigned, make them
            # fully homoplasious.
            AL <- A[[L]]  
            allelesL <- as.character(.unal1loc(object, xsamples, L))
            allelesUA <- allelesL[!allelesL %in% dimnames(AL)[[2]]]
            if(length(allelesUA)>0){
                toadd <- matrix(1, nrow=length(Lindex), ncol=length(allelesUA),
                                dimnames=list(NULL, allelesUA))
                AL <- cbind(AL, toadd)
            }
        for(s in xsamples){
            thisgen <- as.character(Genotype(object, s, L))
            thisA <- AL[,thisgen, drop=FALSE]
            # too many alleles overall?
            if(length(thisgen) > nSG[L]*SGp[L]){
                next # just skip, or see if there is anything useful?
            }

            SGcomplete <- rep(FALSE, nSG[L]) # mark where genotype is known
                                             # for certain
            # set up list to contain genotypes
            alBySG <- rep(list(character(0)), nSG[L])
            # give up on trying to assign?
            stumped <- FALSE
            # Allele assigment loop
           while(sum(!SGcomplete)>0 & !stumped){
                stumped <- TRUE
                # sub genomes that don't yet have definite genotypes:
                SGtoassign <- (1:nSG[L])[!SGcomplete]
            # see if any isoloci have max number of alleles w/o homoplasy.
            for(sg in SGtoassign){
                sgAl <- thisgen[thisA[sg,]==1 &    # Belong to this subgenome,
                                colSums(thisA)==1] # and non-homoplasious.
                if(length(sgAl)>=SGp[L]){
                    alBySG[[sg]] <- sgAl # assign these alleles to this isolocus
                    SGcomplete[sg] <- TRUE # mark as done
                    stumped <- FALSE
                    # Change thisA to reflect that this isolocus can't have
                    # any other alleles.
                    thisA[sg,-match(sgAl, thisgen)] <- 0
                } else {
                # Also check if there is only one allele that could be assigned
                    # to this isolocus.
                    sgAl <- thisgen[thisA[sg,]==1]
                    if(length(sgAl)==1 || # only one allele could be assigned,
                      # or there's no chance of homoplasy (tetrasomic or higher)
                       all(colSums(thisA[,sgAl])==1)){
                         alBySG[[sg]] <- sgAl
                         SGcomplete[sg] <- TRUE
                         stumped <- FALSE
                    }
                }
            }}

            # check if any have >SGp alleles, and whether we are allowing
            # aneuploidy; change ploidy values if necessary.
            nAl <- sapply(alBySG, length)
            if(!all(nAl <= SGp[L])){
                if(allowAneuploidy){
                    newploidies <- nAl
                    while(sum(newploidies) < nSG[L]*SGp[L]){
                newploidies[newploidies==min(newploidies[newploidies!=0])] <-
                  newploidies[newploidies==min(newploidies[newploidies!=0])] + 1
                    }
                    Ploidies(result)[s,Lindex] <- newploidies
                } else {
                    next # leave as missing data if no aneuploidy allowed
                }
            }

            # convert alleles back to numbers if necessary
            if(LisInt){
                alBySG <- lapply(alBySG, as.integer)
            } else {
                if(LisNum){
                    alBySG <- lapply(alBySG, as.numeric)
                }
            }
            # convert zero-length genotypes to missing data
            alBySG[nAl==0] <- Missing(result)
            
            # add genotypes to genambig object
            Genotypes(result, samples=s,loci=Lindex) <- alBySG
        } # end of sample loop
        } # end of "else" clause for if allele assignments need to be done
    } # end of locus loop
    return(result)
} # end of function
                           
# polysat should have its own theme song to the tune of "Smelly Cat" from
# Friends.
