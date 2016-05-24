# Function to find unique alleles at a locus
.unal1loc <- function(object, samples, locus){
  if(!is(object, "genambig"))
    stop("Object must be of class genambig")
  if(missing(samples)) samples <- Samples(object)
  if(missing(locus)) stop("locus argument needed.")
  
  al <- sort(unique(stack(object@Genotypes[samples,locus])$values))
  al <- al[al != Missing(object)]
  return(al)
}

genambig.to.genbinary <- function(object, samples=Samples(object), loci=Loci(object)){
    # set up the object that will ultimately be returned
    objectN <- new("genbinary", samples, loci)

    # fill in the slots that will be identical
    PopNames(objectN) <- PopNames(object)
    if(!all(is.na(PopInfo(object)[samples]))){
      PopInfo(objectN) <- PopInfo(object)[samples]
    }
    Usatnts(objectN) <- Usatnts(object)[loci]
    Description(objectN) <- Description(object)
    if(is(object@Ploidies, "ploidyone")){
      objectN@Ploidies <- object@Ploidies
    }
    if(is(object@Ploidies, "ploidysample")){
      objectN@Ploidies <- new("ploidysample", samples=samples)
      Ploidies(objectN) <- Ploidies(object)[samples]
    }
    if(is(object@Ploidies, "ploidylocus")){
      objectN@Ploidies <- new("ploidylocus", loci=loci)
      Ploidies(objectN) <- Ploidies(object)[loci]
    }
    if(is(object@Ploidies, "ploidymatrix")){
      Ploidies(objectN) <- Ploidies(object)[samples, loci]
    }

    # find all unique alleles for each locus
    locvector<-c()
    allelevector<-c()
    for(L in loci){
      # skip if there are no genotypes
      if(all(isMissing(object, samples, L))) next
      # find the alleles
        thesealleles <- .unal1loc(object, samples, L)
        locvector <- c(locvector, rep(L, length(thesealleles)))
        allelevector <- c(allelevector, thesealleles)
    }

    # Build data frame of locus and allele information
    colinfo <- data.frame(Loci=locvector, Alleles=allelevector,
                          stringsAsFactors=FALSE)
    # Set up matrix to ultimately put in the Genotypes slot
    domdata <- matrix(nrow = length(samples), ncol=dim(colinfo)[1],
                      dimnames=list(samples, paste(locvector, allelevector,
                      sep=".")))
    # Fill the matrix
    for(m in 1:dim(colinfo)[1]){
        for(s in samples){
            if(isMissing(object, s, colinfo[m,1])){
                domdata[s,m] <- Missing(objectN)
            } else {
                if(colinfo[[m,2]] %in% Genotype(object, s, colinfo[m,1])){
                    domdata[s,m] <- Present(objectN)
                } else { domdata[s,m] <- Absent(objectN) }
            }
        }
    }

    # Fill in genotypes slot
    Genotypes(objectN) <- domdata

    # return the converted object
    return(objectN)
}

genbinary.to.genambig <- function(object, samples = Samples(object),
                                  loci = Loci(object)){
    # set up new genambig object
    objectN <- new("genambig", samples, loci)

    # fill in the slots that will be identical
    PopNames(objectN) <- PopNames(object)
    if(!all(is.na(PopInfo(object)[samples]))){
      PopInfo(objectN) <- PopInfo(object)[samples]
    }
    Usatnts(objectN) <- Usatnts(object)[loci]
    Description(objectN) <- Description(object)
    if(is(object@Ploidies, "ploidyone")){
      objectN@Ploidies <- object@Ploidies
    }
    if(is(object@Ploidies, "ploidysample")){
      objectN@Ploidies <- new("ploidysample", samples=samples)
      Ploidies(objectN) <- Ploidies(object)[samples]
    }
    if(is(object@Ploidies, "ploidylocus")){
      objectN@Ploidies <- new("ploidylocus", loci=loci)
      Ploidies(objectN) <- Ploidies(object)[loci]
    }
    if(is(object@Ploidies, "ploidymatrix")){
      Ploidies(objectN) <- Ploidies(object)[samples, loci]
    }

    # Go through matrix one locus at a time to get genotypes
    for(L in loci){
        thesegen <- as.matrix(Genotypes(object, samples, L))
        # Go through each allele
        for(n in 1:dim(thesegen)[2]){
            allele <- strsplit(dimnames(thesegen)[[2]][n], split=".",
                               fixed=TRUE)[[1]][2]
            allele <- as.integer(allele)
            # Look for allele in each sample and write to new object
            for(s in samples){
                if(thesegen[s,n] == Present(object)){
                    if(isMissing(objectN, s, L)){
                        Genotype(objectN, s, L) <- allele
                    } else {
                        Genotype(objectN, s, L) <- c(Genotype(objectN, s, L),
                                                     allele)
                    }
                }
            }
        }
    }

    # return the new genambig object
    return(objectN)
}
