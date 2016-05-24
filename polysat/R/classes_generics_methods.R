## load library just to make things work properly for now
#library(polysat)

## Set classes for containing ploidies
setClass("ploidysuper", contains="VIRTUAL")
setClass("ploidymatrix", representation(pld="matrix"), contains="ploidysuper",
         validity=function(object){
           failures <- character(0)

           if(is.null(dimnames(object@pld)[[1]]))
             failures <- c(failures, "Sample names missing.")
           if(is.null(dimnames(object@pld)[[2]]))
             failrues <- c(failures, "Locus names misisng.")

           if(length(failures)==0){
             return(TRUE)
           } else {
             return(failures)
           }
           })
setClass("ploidysample", representation(pld="integer"), contains="ploidysuper",
         validity=function(object){
           failures <- character(0)

           if(is.null(names(object@pld)))
              failures <- c(failures, "Sample names missing.")

           if(length(failures==0)){
             return(TRUE)
           } else {
             return(failures)
           }
           })
setClass("ploidylocus", representation(pld="integer"), contains="ploidysuper",
         validity=function(object){
           failures <- character(0)

           if(is.null(names(object@pld)))
             failures <- c(failures, "Locus names missing.")

           if(length(failures==0)){
             return(TRUE)
           } else {
             return(failures)
           }
           })
setClass("ploidyone", representation(pld="integer"), contains="ploidysuper",
         validity=function(object){
           failures <- character(0)

           if(!is.null(names(object@pld)))
             failures <- c(failures, "Ploidy vector should not be named.")
           if(length(object@pld)!=1)
             failures <- c(failures, "Ploidy vector should have length 1.")

           if(length(failures==0)){
             return(TRUE)
           } else {
             return(failures)
           }
         })
setMethod("initialize", signature(.Object="ploidymatrix"),
          function(.Object, samples, loci, ...){
            .Object@pld <- matrix(as.integer(NA), nrow=length(samples),
                                 ncol=length(loci), dimnames=list(samples, loci))
            callNextMethod(.Object, ...)})
setMethod("initialize", signature(.Object="ploidysample"),
          function(.Object, samples, loci, ...){
            ploidies <- as.integer(rep(NA, length(samples)))
            names(ploidies) <- samples
            .Object@pld <- ploidies
            callNextMethod(.Object, ...)})
setMethod("initialize", signature(.Object="ploidylocus"),
          function(.Object, samples, loci, ...){
            ploidies <- as.integer(rep(NA, length(loci)))
            names(ploidies) <- loci
            .Object@pld <- ploidies
            callNextMethod(.Object, ...)})
setMethod("initialize", signature(.Object="ploidyone"),
          function(.Object, samples, loci, ...){
            .Object@pld <- as.integer(NA)
            callNextMethod(.Object, ...)})
# generics to get and replace ploidies
setGeneric("pld", function(object, samples, loci) standardGeneric("pld"))
setGeneric("pld<-", function(object, value) standardGeneric("pld<-"))
# methods to get and replace ploidies
setMethod("pld", signature(object = "ploidyone"),
          function(object,samples,loci) return(object@pld))
setMethod("pld", signature(object = "ploidymatrix"),
          function(object,samples,loci){
            if(missing(samples)) samples <- dimnames(object@pld)[[1]]
            if(missing(loci)) loci <- dimnames(object@pld)[[2]]
            return(object@pld[samples,loci, drop=FALSE])
          })
setMethod("pld", signature(object = "ploidysample"),
          function(object,samples,loci){
            if(missing(samples)) samples <- names(object@pld)
            return(object@pld[samples])
          })
setMethod("pld", signature(object = "ploidylocus"),
          function(object,samples,loci){
            if(missing(loci)) loci <- names(object@pld)
            return(object@pld[loci])
          })
setReplaceMethod("pld", "ploidyone", function(object, value){
  if(length(value) > 1) stop("Only one ploidy allowed for this format.")
  object@pld <- unname(as.integer(value))
  object
})
setReplaceMethod("pld", "ploidymatrix", function(object, value){
  dn <- dimnames(object@pld)
  object@pld[dn[[1]],dn[[2]]] <- as.integer(value)
  object
})
setReplaceMethod("pld", "ploidysample", function(object, value){
  object@pld[names(object@pld)] <- as.integer(value)
  object
})
setReplaceMethod("pld", "ploidylocus", function(object, value){
  object@pld[names(object@pld)] <- as.integer(value)
  object
})
# generic and methods to test collapsibility of ploidies
setGeneric("plCollapse",
           function(object, na.rm, returnvalue) standardGeneric("plCollapse"))
setMethod("plCollapse", signature(object="ploidymatrix", na.rm="logical",
                                  returnvalue="logical"),
          function(object, na.rm, returnvalue){
            if(na.rm){
              removeNA <- function(x) return(x[!is.na(x)]) # pulls NA out of vectors
            }
            done <- FALSE # has it been collapsed yet?
            # test to see if there is one ploidy for the whole set
            u <- unique(as.vector(pld(object)))
            if(na.rm) u <- removeNA(u)
            if(length(u)==1){
              done <- TRUE
              if(returnvalue){
                result <- new("ploidyone")
                pld(result) <- u
                return(result)
              } else {
                return(TRUE)
              }
            }
            if(!done){
              # test to see if there is one ploidy per sample
              u <- sapply(data.frame(t(pld(object))), unique, simplify=FALSE)
              if(na.rm) u <- sapply(u, removeNA, simplify=FALSE)
              if(all(sapply(u, length, simplify=TRUE)==1)){
                done <- TRUE
                if(returnvalue){
                  result <- new("ploidysample", samples=dimnames(pld(object))[[1]])
                  pld(result) <- unlist(u)
                  return(result)
                } else {
                  return(TRUE)
                }
              }
            }
            if(!done){
              # test to see if there is one ploidy per locus
              u <- sapply(data.frame(pld(object)), unique, simplify=FALSE)
              if(na.rm) u <- sapply(u, removeNA, simplify=FALSE)
              if(all(sapply(u, length, simplify=TRUE)==1)){
                done <- TRUE
                if(returnvalue){
                  result <- new("ploidylocus", loci=dimnames(pld(object))[[2]])
                  pld(result) <- unlist(u)
                  return(result)
                } else {
                  return(TRUE)
                }
              }
            }
            if(!done) return(FALSE) # cannot be collapsed
          })
setMethod("plCollapse", signature(object="ploidysample", na.rm="logical",
                                  returnvalue="logical"),
          function(object, na.rm, returnvalue){
            u <- unique(pld(object))
            if(na.rm) u <- u[!is.na(u)]
            if(length(u) != 1){
              return(FALSE)
            } else {
              if(returnvalue){
                result <- new("ploidyone")
                pld(result) <- u
                return(result)
              } else {
                return(TRUE)
              }
            }
          })
setMethod("plCollapse", signature(object="ploidylocus", na.rm="logical",
                                  returnvalue="logical"),
          function(object, na.rm, returnvalue){
            u <- unique(pld(object))
            if(na.rm) u <- u[!is.na(u)]
            if(length(u) != 1){
              return(FALSE)
            } else {
              if(returnvalue){
                result <- new("ploidyone")
                pld(result) <- u
                return(result)
              } else {
                return(TRUE)
              }
            }
          })

##                                 gendata Class                             ##
## Contains info about missing data symbol, microsat repeats, sample ploidy, ##
## and populations.                                                          ##
## Is a superclass for other classes containing genotype data.               ##
setClass("gendata", representation(Description="character", Missing="ANY",
                                   Usatnts="integer", Ploidies="ploidysuper",
                                   PopInfo="integer", PopNames="character"),
         # make a virtual class?  (include "VIRTUAL" in contains)
         prototype(Description="Insert dataset description here.",
                   Missing=as.integer(-9)),
         validity=function(object){
             failures <- character(0)
             # check size and names of Ploidies with Usatnts and PopInfo
             if(is(object@Ploidies, "ploidymatrix")){
               np <- dimnames(pld(object@Ploidies))
               if(!identical(np[[1]],names(object@PopInfo)))
                 failures <- c(failures,
                               "Sample names inconsistent between Plodies and PopInfo")
               if(!identical(np[[2]],names(object@Usatnts)))
                 failures <- c(failures,
                               "Locus names inconsistent between Ploidies and Usatnts")
             }
             if(is(object@Ploidies,"ploidysample")){
               np <- names(pld(object@Ploidies))
               if(!identical(np, names(object@PopInfo)))
                 failures <- c(failures,
                        "Names of Ploidies do not match names of PopInfo.")
             }
             if(is(object@Ploidies,"ploidylocus")){
               np <- names(pld(object@Ploidies))
               if(!identical(np, names(object@Usatnts)))
                 failures <- c(failures,
                               "Names of Ploidies do not match names of Usatnts")
             }

             # check to see that Missing only has one element
             if(length(object@Missing) != 1) failures <- c(failures,
                 "Missing needs to have only one element.")

             # check to see that each number in PopInfo has an element in PopNames
             okpopnums <- 1:length(object@PopNames) #numbers that go with a PopName
             usedpopnums <- unique(object@PopInfo) #numbers in PopInfo
             usedpopnums <- usedpopnums[!is.na(usedpopnums)] # ignore NAs
             if(length(usedpopnums) > length(usedpopnums[usedpopnums %in% okpopnums]))
                failures <- c(failures,
                          "PopInfo contains integers that don't index PopNames."
                              )

             # check to see that all sample and locus names are unique
             if(length(unique(names(object@PopInfo))) < length(object@PopInfo)){
                 failures <- c(failures, "Not all sample names are unique.")
             }
             if(length(unique(names(object@Usatnts))) < length(object@Usatnts)){
                 failures <- c(failures, "Not all locus names are unique.")
             }

             # return TRUE or failures
             if(length(failures) == 0){
                 return(TRUE)
             } else {
                 return(failures)
             }
         })

##                                  genambig class                                ##
##                               Subclass of gendata                              ##
## Contains two-dimensional list of vectors to contain genotypes where            ##
## allele copy number is not necessarily known.                                   ##
setClass("genambig", representation(Genotypes="array"), contains="gendata",
          validity = function(object){
    failures <- character(0)

    # do validity testing from gendata
    gendatatest <- getValidity(getClassDef("gendata"))(object)
    if(!identical(TRUE, gendatatest)){
        failures <- c(failures, gendatatest)
    }

    # check to see that the Genotypes list is two-dimensional
    if(length(dim(object@Genotypes)) != 2)
       failures <- c(failures, "Genotypes should be a two-dimensional list.")

    # include a check to see that all list elements are vectors?

    # check to see that sample names match in Genotypes, Ploidies, PopInfo
    ng <- dimnames(object@Genotypes)
    if(is(object@Ploidies,"ploidymatrix")){
      np <- dimnames(pld(object@Ploidies))
      if(!identical(np[[1]],ng[[1]]))
        failures <- c(failures,
                      "Sample names inconsistent between Plodies and Genotypes")
      if(!identical(np[[2]],ng[[2]]))
        failures <- c(failures,
                      "Locus names inconsistent between Ploidies and Genotypes")
    }
    if(is(object@Ploidies,"ploidysample")){
      np <- names(pld(object@Ploidies))
      if(!identical(np, ng[[1]]))
        failures <- c(failures,
               "Names of Ploidies do not match Sample names in Genotypes")
    }
    if(is(object@Ploidies,"ploidylocus")){
      np <- names(pld(object@Ploidies))
      if(!identical(np, ng[[2]]))
        failures <- c(failures,
               "Names of Ploidies do not match Locus names in Genotypes")
    }
    if(!identical(ng[[1]], names(object@PopInfo)))
        failures <- c(failures,
                      "Sample names are not the same in Genotypes and PopInfo."
                      )

    # check to see that locus names match in Genotypes, Usatnts
    if(!identical(ng[[2]], names(object@Usatnts)))
        failures <- c(failures,
                      "Locus names are not the same in Genotypes and Usatnts."
                      )

    # return TRUE or failures
    if(length(failures) == 0){
        return(TRUE)
    } else {
        return(failures)
    }
})

## genbinary class
## subclass of Gendata
## Genotypes stored as a matrix of symbols indicating the presence and
## absence of alleles.
setClass("genbinary", representation(Genotypes = "matrix", Present = "ANY",
                                     Absent = "ANY"), contains="gendata",
#         where = where,
         prototype(Present = as.integer(1), Absent = as.integer(0)),
         validity = function(object){
             failures <- character(0)

             # do validity testing from gendata
             gendatatest <- getValidity(getClassDef("gendata"))(object)
             if(!identical(TRUE, gendatatest)){
                 failures <- c(failures, gendatatest)
             }

             # check that there is only one element for Present and Absent
             if(length(object@Present)!=1) failures <- c(failures,
                      "object@Present must be a vector of length 1.")
             if(length(object@Absent)!=1) failures <- c(failures,
                      "object@Absent must be a vector of length 1.")

             # check that Present, Absent, and Missing are different
             pam <- c(object@Present, object@Absent, object@Missing)
             if(length(unique(pam))!=3) failures <- c(failures,
                      "Different values required for Present, Absent, and Missing.")

             # check that all symbols in Genotypes are in PAM
             matsym <- unique(object@Genotypes)
             if(length(matsym) > 0){
                 if(FALSE %in% (matsym %in% pam)) failures <- c(failures,
                    "All genotype symbols must match Present, Absent, or Missing.")
             }

             # return TRUE or failures
             if(length(failures) == 0){
                 return(TRUE)
             } else {
                 return(failures)
             }
         }
     )

### New generic functions for classes in polysat ###
# generic function to get sample names
setGeneric("Samples", function(object, populations, ploidies)
           standardGeneric("Samples"))
# generic to replace sample names
setGeneric("Samples<-", function(object, value) standardGeneric("Samples<-"))
# generic to get locus names
setGeneric("Loci", function(object, usatnts, ploidies) standardGeneric("Loci"))
# generic to replace locus names
setGeneric("Loci<-", function(object, value) standardGeneric("Loci<-"))
# generics to get and replace population identities
setGeneric("PopInfo", function(object) standardGeneric("PopInfo"))
setGeneric("PopInfo<-", function(object, value) standardGeneric("PopInfo<-"))
# generics to get and replace population names
setGeneric("PopNames", function(object) standardGeneric("PopNames"))
setGeneric("PopNames<-", function(object, value) standardGeneric("PopNames<-"))
# generics to get and replace ploidies
setGeneric("Ploidies", function(object, samples, loci) standardGeneric("Ploidies"))
setGeneric("Ploidies<-", function(object, value) standardGeneric("Ploidies<-"))
# generics to get and replace usatnts
setGeneric("Usatnts", function(object) standardGeneric("Usatnts"))
setGeneric("Usatnts<-", function(object, value) standardGeneric("Usatnts<-"))
# generics to delete samples and loci
setGeneric("deleteSamples", function(object, samples)
           standardGeneric("deleteSamples"))
setGeneric("deleteLoci", function(object, loci) standardGeneric("deleteLoci"))
# generics to get and replace missing data code
setGeneric("Missing", function(object) standardGeneric("Missing"))
setGeneric("Missing<-", function(object, value) standardGeneric("Missing<-"))
# generics to get and replace description
setGeneric("Description", function(object) standardGeneric("Description"))
setGeneric("Description<-", function(object,value)
           standardGeneric("Description<-"))
# generics to get and replace population number by name
setGeneric("PopNum", function(object, popname) standardGeneric("PopNum"))
setGeneric("PopNum<-", function(object, popname, value)
           standardGeneric("PopNum<-"))
# generics to get and edit genotype data
setGeneric("Genotype", function(object, sample, locus) standardGeneric("Genotype"))
setGeneric("Genotype<-", function(object, sample, locus, value)
           standardGeneric("Genotype<-"))
setGeneric("Genotypes", function(object, samples = Samples(object), loci = Loci(object))
           standardGeneric("Genotypes"))
setGeneric("Genotypes<-", function(object, samples = Samples(object),
                                   loci = Loci(object), value)
           standardGeneric("Genotypes<-"))
# generic to determine if a genotype is missing
setGeneric("isMissing", function(object, samples = Samples(object),
                                 loci = Loci(object))
           standardGeneric("isMissing"))
# generic to print genotypes to console
setGeneric("viewGenotypes",
           function(object, samples = Samples(object), loci = Loci(object))
           standardGeneric("viewGenotypes"))
# generic to edit genotypes
setGeneric("editGenotypes",
           function(object, maxalleles = max(Ploidies(object)),
                    samples = Samples(object),
                    loci = Loci(object)) standardGeneric("editGenotypes"))
# generic to estimate ploidy
setGeneric("estimatePloidy",
           function(object, extrainfo, samples = Samples(object), loci = Loci(object))
           standardGeneric("estimatePloidy"))
# generics to get and replace symbols for allele present/absent
setGeneric("Present", function(object) standardGeneric("Present"))
setGeneric("Present<-", function(object, value) standardGeneric("Present<-"))
setGeneric("Absent", function(object) standardGeneric("Absent"))
setGeneric("Absent<-", function(object, value) standardGeneric("Absent<-"))

#### gendata methods
# initialization for a gendata object
setMethod("initialize",
    signature(.Object = "gendata"),
    function (.Object, samples, loci, Missing, ...)
    {
        if(missing(samples)) samples <- c("ind1","ind2")
        if(missing(loci)) loci <- c("loc1","loc2")

        # make a vector to contain repeat lengths
        usatnts <- as.integer(rep(NA, length(loci)))
        names(usatnts) <- loci
        .Object@Usatnts <- usatnts

        # make a matrix to contain ploidy
        .Object@Ploidies <- new("ploidymatrix",samples=samples,loci=loci)

        # make a vector to contain population identity
        popinfo <- as.integer(rep(NA, length(samples)))
        names(popinfo) <- samples
        .Object@PopInfo <- popinfo

        # put in the missing data symbol
        if(missing(Missing)) Missing <- as.integer(-9)
        .Object@Missing <- Missing

        # call default method
        callNextMethod(.Object, ...)
    }
)

## methods to get sample names
# All samples if only gendata object is given
setMethod("Samples", signature(object = "gendata", populations = "missing",
                               ploidies = "missing"),
          function(object){return(names(object@PopInfo))})
# Just for a subset of population names
setMethod("Samples", signature(object = "gendata", populations = "character",
                               ploidies = "missing"),
          function(object, populations){
              pops <- match(populations, object@PopNames)
              return(names(object@PopInfo)[object@PopInfo %in% pops])
              })
# Population names and ploidies
setMethod("Samples", signature(object = "gendata", populations = "character",
                               ploidies = "numeric"),
          function(object, populations, ploidies){
            if(!is(object@Ploidies,"ploidysample"))
        stop("Ploidies argument only valid if ploidies in dataset are indexed by sample.")
              ploidies <- as.integer(ploidies)
              pops <- match(populations, object@PopNames)
              popsam <- names(object@PopInfo)[object@PopInfo %in% pops]
              ploisam <- names(pld(object@Ploidies))[pld(object@Ploidies) %in% ploidies]
              return(popsam[popsam %in% ploisam])
              })
# Just a subset of population numbers
setMethod("Samples", signature(object = "gendata", populations = "numeric",
                               ploidies = "missing"),
          function(object, populations){
              populations <- as.integer(populations)
              return(names(object@PopInfo)[object@PopInfo %in% populations])
              })
# Population numbers and ploidies
setMethod("Samples", signature(object = "gendata", populations = "numeric",
                               ploidies = "numeric"),
          function(object, populations, ploidies){
            if(!is(object@Ploidies,"ploidysample"))
        stop("Ploidies argument only valid if ploidies in dataset are indexed by sample.")

              ploidies <- as.integer(ploidies)
              populations <- as.integer(populations)
              popsam <- names(object@PopInfo)[object@PopInfo %in% populations]
              ploisam <- names(pld(object@Ploidies))[pld(object@Ploidies) %in% ploidies]
              return(popsam[popsam %in% ploisam])
              })
# Just ploidies
setMethod("Samples", signature(object = "gendata", populations = "missing",
                               ploidies = "numeric"),
          function(object, ploidies){
            if(!is(object@Ploidies,"ploidysample"))
        stop("Ploidies argument only valid if ploidies in dataset are indexed by sample.")

            ploidies <- as.integer(ploidies)
            return(names(pld(object@Ploidies))[pld(object@Ploidies) %in% ploidies])
              })

##                     Replacement method for sample names                   ##
setReplaceMethod("Samples", "gendata", function(object, value){
    # check that all sample names are unique
    if(length(unique(value)) < length(value)){
        stop("Not all sample names are unique")
    }
  if(is(object@Ploidies,"ploidysample"))
    names(object@Ploidies@pld) <- value
  if(is(object@Ploidies,"ploidymatrix"))
    dimnames(object@Ploidies@pld)[[1]] <- value
  names(object@PopInfo) <- value
  object
})

## Methods to get locus names.
# Return all locus names
setMethod("Loci", signature(object = "gendata", usatnts = "missing", ploidies="missing"),
          function(object){
              return(names(object@Usatnts))
              })
# Return locus names for a certain subset of repeat types
setMethod("Loci", signature(object = "gendata", usatnts = "numeric", ploidies="missing"),
          function(object, usatnts){
              usatnts <- as.integer(usatnts)
              return(names(object@Usatnts)[object@Usatnts %in% usatnts])
          })
# Return locus names for a certain subset of ploidies
setMethod("Loci", signature(object = "gendata", usatnts = "missing", ploidies="numeric"),
          function(object, ploidies){
            if(!is(object@Ploidies,"ploidylocus"))
          stop("Ploidies argument only valid if ploidies in dataset are indexed by locus")

            ploidies <- as.integer(ploidies)
            loci <- names(pld(object@Ploidies))[pld(object@Ploidies) %in% ploidies]

            return(loci)
          })
# Return locus names for a subset of ploidies and repeat types
setMethod("Loci", signature(object = "gendata", usatnts = "numeric", ploidies="numeric"),
          function(object, usatnts, ploidies){
            if(!is(object@Ploidies,"ploidylocus"))
          stop("Ploidies argument only valid if ploidies in dataset are indexed by locus")

            ploidies <- as.integer(ploidies)
            ploci <- names(pld(object@Ploidies))[pld(object@Ploidies) %in% ploidies]

            usatnts <- as.integer(usatnts)
            uloci <- names(object@Usatnts)[object@Usatnts %in% usatnts]

            return(uloci[uloci %in% ploci])
          })

## Replacement method for locus names
setReplaceMethod("Loci", "gendata", function(object, value){
    # check that all locus names are unique
    if(length(unique(value)) < length(value)){
        stop("Not all locus names are unique")
    }
  # Replace Ploidies names if necessary
  if(is(object@Ploidies,"ploidylocus"))
    names(object@Ploidies@pld) <- value
  if(is(object@Ploidies,"ploidymatrix"))
    dimnames(object@Ploidies@pld)[[2]] <- value

  # Replace Usatnts names
  names(object@Usatnts) <- value
  object
})

## Methods to get and replace population identities
setMethod("PopInfo", signature(object = "gendata"),
          function(object) return(object@PopInfo))
setReplaceMethod("PopInfo", "gendata", function(object, value){
    # add values to PopInfo slot
    object@PopInfo[names(object@PopInfo)] <- as.integer(value)
    # make sure there are enough PopNames
    if(!all(is.na(object@PopInfo)) &&
        length(object@PopNames) < max(object@PopInfo, na.rm=TRUE)){
        missingnames <- (length(object@PopNames) + 1):max(object@PopInfo,
                                                          na.rm=TRUE)
        object@PopNames[missingnames] <- paste("Pop", missingnames,sep="")
    }
    # return object
    object
})

## Methods to get and replace population names
setMethod("PopNames", signature(object = "gendata"),
          function(object) return(object@PopNames))
setReplaceMethod("PopNames", "gendata", function(object, value){
    object@PopNames <- value
    object
})

## Methods to get and replace ploidies
setMethod("Ploidies", signature(object = "gendata", samples="ANY", loci="ANY"),
          function(object, samples, loci){
            return(pld(object@Ploidies, samples=samples, loci=loci))
            })
setReplaceMethod("Ploidies", "gendata", function(object, value){
    if(!all(is.na(as.integer(value)) | as.integer(value) >= 0))
        stop("Ploidy less than zero not allowed.")
    if(!all(is.na(as.integer(value)) | as.integer(value) <= 12))
        warning("Large ploidies (>12) detected.")
    pld(object@Ploidies) <- value
    object
})
# Function to change the ploidy format between the four types
reformatPloidies <- function(object, output="collapse", na.rm=FALSE,
                             erase=FALSE){

  if(!output %in% c("matrix", "sample", "locus", "one", "collapse"))
    stop("Unrecognized output argument.  See reformatPloidies documentation.")

  if(erase){ # just make a new ploidies object filled with NA
    if(output=="collapse") output <- "one"
    object@Ploidies <- new(paste("ploidy",output,sep=""),
                           samples=Samples(object),
                           loci=Loci(object))
  } else {

    ## If you start with ploidymatrix - collapse or do nothing
  if(is(object@Ploidies,"ploidymatrix") && output %in% c("sample","locus","one",
                                                         "collapse")){
    pl <- plCollapse(object@Ploidies, na.rm=na.rm, returnvalue=TRUE)
    if(is(pl, "ploidysuper")){ # if the ploidy was collapsible
      if(output=="collapse" || is(pl, paste("ploidy",output,sep=""))){
        # if collapsed ploidy is the right format, we're done
        object@Ploidies <- pl
      } else { # if collapsed ploidy is not right format
        if(is(pl, "ploidyone")){ # expand from ploidyone
          pl2 <- new(paste("ploidy",output,sep=""),
                     samples=Samples(object),
                     loci=Loci(object))
          pld(pl2) <- pld(pl)
          object@Ploidies <- pl2
        } else {
          stop("Ploidies not collapsible to desired format.")
        }
      }
    } else { # if ploidy was not collapsible at all
      if(output != "collapse") stop("Ploidies not collapsible to desired format.")
    }
  }

  # if you start with ploidysample and want to collapse
  if(is(object@Ploidies,"ploidysample") && output %in% c("one","collapse","locus")){
    pl <- plCollapse(object@Ploidies, na.rm=na.rm, returnvalue=TRUE)
    if(is(pl, "ploidyone")){ # if it is collapsible
      if(output %in% c("collapse","one")){
        object@Ploidies <- pl
      }
      if(output=="locus"){
        pl2 <- new("ploidylocus", loci=Loci(object))
        pld(pl2) <- pld(pl)
        object@Ploidies <- pl2
      }
    } else { # if it is not collapsible
      if(output != "collapse") stop("Ploidies not collapsible to desired format.")
    }
  }

  # if you start with ploidysample and want to expand
  if(is(object@Ploidies,"ploidysample") && output=="matrix"){
    pl <- new("ploidymatrix", samples=Samples(object), loci=Loci(object))
    pld(pl) <- pld(object@Ploidies)
    object@Ploidies <- pl
  }

  # if you start with ploidylocus and want to collapse
  if(is(object@Ploidies,"ploidylocus") && output %in% c("sample","one","collapse")){
    pl <- plCollapse(object@Ploidies, na.rm=na.rm, returnvalue=TRUE)
    if(is(pl, "ploidyone")){ # if it is collapsible
      if(output %in% c("collapse","one")){
        object@Ploidies <- pl
      }
      if(output=="sample"){
        pl2 <- new("ploidysample", samples=Samples(object))
        pld(pl2) <- pld(pl)
        object@Ploidies <- pl2
      }
    } else { # if it is not collapsible
      if(output != "collapse") stop("Ploidies not collapsible to desired format.")
    }
  }

  # if you start with ploidylocus and want to expand
  if(is(object@Ploidies,"ploidylocus") && output == "matrix"){
    pl <- new("ploidymatrix", samples=Samples(object), loci=Loci(object))
    pld(pl) <- rep(pld(object@Ploidies), each=length(Samples(object)))
    object@Ploidies <- pl
  }

  # if you start with ploidyone and want to expand
  if(is(object@Ploidies,"ploidyone") && output %in% c("matrix","sample","locus")){
    pl <- new(paste("ploidy",output,sep=""),
                     samples=Samples(object), loci=Loci(object))
    pld(pl) <- pld(object@Ploidies)
    object@Ploidies <- pl
  }
  }

  return(object)
}


## Methods to get and replace repeat lengths
setMethod("Usatnts", signature(object = "gendata"),
          function(object) return(object@Usatnts))
setReplaceMethod("Usatnts", "gendata", function(object, value){
    object@Usatnts[names(object@Usatnts)] <- as.integer(value)
    object
})

## summary method for gendata
setMethod("summary", "gendata", function(object){
    cat(paste(length(Samples(object)), "samples,", length(Loci(object)), "loci."),
        paste(length(unique(PopInfo(object))), "populations."),
        sep="\n")
    cat("Ploidies:", unique(as.vector(Ploidies(object))))
    cat("\n", sep="")
    cat("Length(s) of microsatellite repeats:", unique(Usatnts(object)))
    cat("\n", sep="")
})

## Methods to delete samples and loci
setMethod("deleteSamples", "gendata", function(object, samples){
    samtouse <- Samples(object)[!Samples(object) %in% samples]
    object@PopInfo <- object@PopInfo[samtouse]
    if(is(object@Ploidies,"ploidysample")){
      object@Ploidies@pld <- object@Ploidies@pld[samtouse]
    }
    if(is(object@Ploidies,"ploidymatrix")){
      object@Ploidies@pld <- object@Ploidies@pld[samtouse,, drop=FALSE]
    }
    return(object)
})

setMethod("deleteLoci", "gendata", function(object, loci){
    loctouse <- Loci(object)[!Loci(object) %in% loci]
    object@Usatnts <- object@Usatnts[loctouse]
    if(is(object@Ploidies,"ploidylocus")){
      object@Ploidies@pld <- object@Ploidies@pld[loctouse]
    }
    if(is(object@Ploidies,"ploidymatrix")){
      object@Ploidies@pld <- object@Ploidies@pld[,loctouse, drop=FALSE]
    }
    return(object)
})

## Subscripting method
setMethod("[", signature(x="gendata", i="ANY", j="ANY"), function(x, i, j){
    # i is samples, j is loci
    x@PopInfo <- x@PopInfo[i]
    if(is(x@Ploidies, "ploidymatrix")){
      x@Ploidies@pld <- x@Ploidies@pld[i,j, drop=FALSE]
    }
    if(is(x@Ploidies, "ploidysample")){
      x@Ploidies@pld <- x@Ploidies@pld[i]
    }
    if(is(x@Ploidies, "ploidylocus")){
      x@Ploidies@pld <- x@Ploidies@pld[j]
    }
    x@Usatnts <- x@Usatnts[j]
    return(x)
})

## Methods to get and replace missing data code
setMethod("Missing", "gendata", function(object) return(object@Missing))
setReplaceMethod("Missing", "gendata", function(object, value){
    if(length(value) != 1) stop("Assigned value should only have one element.")
    object@Missing <- value
    object
})

## Methods to get and replace description
setMethod("Description", "gendata", function(object) return(object@Description))
setReplaceMethod("Description", "gendata", function(object, value){
    object@Description <- value
    object
})

## Methods to get and replace population numbers by name
setMethod("PopNum", signature(object="gendata", popname="character"),
          function(object, popname){
              return(match(popname, PopNames(object)))
          })
setReplaceMethod("PopNum", signature(object="gendata", popname="character"),
                 function(object, popname, value){
                     value <- as.integer(value)
                     if(length(Samples(object, populations=value))>0 &&
                        length(Samples(object, populations=popname))>0 &&
                        !identical(Samples(object, populations=value),
                                   Samples(object, populations=popname)))
                         cat(paste("Samples already present in population ",
                                   value,".  These have been merged into",
                                   popname," .", sep=""), sep="\n")
                     PopInfo(object)[Samples(object, populations=popname)] <-
                         value
                     PopNames(object)[PopNum(object, popname)] <- NA
                     PopNames(object)[value] <- popname
                     object
                     })

## Method for merging objects
setMethod("merge", signature(x="gendata", y="gendata"),
          function(x, y, objectm, samples, loci, overwrite){
            # error if ploidies are not in the same format
            if(!identical(class(x@Ploidies),class(y@Ploidies)))
              stop("x and y must have Ploidies in the same format.")

              # set up new gendata object if this wasn't called from the method
              # of a subclass.
              if(missing(objectm)){
                  if(missing(samples)) samples <- unique(c(Samples(x), Samples(y)))
                  if(missing(loci)) loci <- unique(c(Loci(x), Loci(y)))
                  objectm <- new("gendata", samples=samples, loci=loci)
              }

              # determine what to do with conflicting data: overwrite or give error
              if(missing(overwrite)) overwrite <- "empty"
              owerror <- ifelse(overwrite != "x" && overwrite != "y", TRUE, FALSE)
              if(overwrite=="x"){
                  objectB <- x
                  objectT <- y
              } else {
                  objectB <- y
                  objectT <- x
              }

            # function definition for merging vectors
            # B and T are vectors, m is final index (Samples(objectm) or Loci(objectm)), x is output
            vectmerge <- function(m, B, T, type){
              x <- rep(NA, length(m)) # set up an empty vector
              names(x) <- m

              b <- m[m %in% names(B)] # subset indexing vectors
              t <- m[m %in% names(T)]

              x[b] <- B[b]
              x[t] <- T[t]

              if(owerror && !identical(x[b], B[b]))
                stop(paste("Conflicting", type,"data; set overwrite to \"x\" or \"y\"."))

              return(x)
            }

              # merge Usatnts
            Usatnts(objectm) <- vectmerge(Loci(objectm), Usatnts(objectB), Usatnts(objectT), "Usatnts")

              # merge Ploidies
            if(is(objectT@Ploidies,"ploidysample")){
              objectm@Ploidies <- new("ploidysample", samples=Samples(objectm))
              Ploidies(objectm) <- vectmerge(Samples(objectm), Ploidies(objectB), Ploidies(objectT),
                                             "Ploidies")
            }
            if(is(objectT@Ploidies,"ploidylocus")){
              objectm@Ploidies <- new("ploidylocus", loci=Loci(objectm))
              Ploidies(objectm) <- vectmerge(Loci(objectm), Ploidies(objectB), Ploidies(objectT),
                                             "Ploidies")
            }
            if(is(objectT@Ploidies,"ploidyone")){
              if(owerror && Ploidies(objectB) != Ploidies(objectT))
                stop("Conflicting Ploidies data. Set overwrite to \"x\" or \"y\".")
              objectm@Ploidies <- new("ploidyone")
              Ploidies(objectm) <- Ploidies(objectT)
            }
            if(is(objectT@Ploidies,"ploidymatrix")){
              a <- Samples(objectm)[Samples(objectm) %in% Samples(objectB)]
              b <- Loci(objectm)[Loci(objectm) %in% Loci(objectB)]
              c <- Samples(objectm)[Samples(objectm) %in% Samples(objectT)]
              d <- Loci(objectm)[Loci(objectm) %in% Loci(objectT)]

              Ploidies(objectm)[a,b] <- Ploidies(objectB)[a,b]
              Ploidies(objectm)[c,d] <- Ploidies(objectT)[c,d]

              if(owerror && !identical(Ploidies(objectm)[a,b], Ploidies(objectB)[a,b]))
                stop("Conflicting Ploidies data. Set overwrite to \"x\" or \"y\".")
            }

              # take description from "top" object or both
              if(owerror && Description(objectB)[1] != Description(objectT)[1]){
                  Description(objectm) <- c(Description(objectT), Description(objectB))
              } else {
                  Description(objectm) <- Description(objectT)
              }

              ## merge PopInfo and PopNames
              # Fill vector of object names
              PopNames(objectm) <- unique(c(PopNames(objectT),PopNames(objectB)))
              # Scenerio for identical sample sets
              if(identical(Samples(objectB), Samples(objectT))){
                  if(owerror && !identical(PopInfo(objectB),PopInfo(objectT))){
                stop("Conflicting data in PopInfo. Set overwrite to \"x\" or \"y\".")
                  } else {
                      PopInfo(objectm) <- PopInfo(objectT)[Samples(objectm)]
                  }
              } else {
              # Scenario for non-overlapping or partially overlapping sample sets
                  for(p in PopNames(objectm)){
                      # fill in values from "bottom" object from this population
                      if(p %in% PopNames(objectB)){
                          samtofill <- Samples(objectB,populations=p)
                          samtofill <- samtofill[samtofill %in% Samples(objectm)]
                          # give error if something would be overwritten
                          if(owerror && (FALSE %in% is.na(PopInfo(objectm)[samtofill])))
                      stop("Conflicting data in PopInfo.  Set overwrite to \"x\" or \"y\".")
                          # only fill in values that are currently NA
                          PopInfo(objectm)[samtofill][is.na(PopInfo(objectm)[samtofill])] <-
                              PopNum(objectm, p)
                      }
                      if(p %in% PopNames(objectT)){
                          samtofill <- Samples(objectT, populations=p)
                          samtofill <- samtofill[samtofill %in% Samples(objectm)]
                          # give error if something would be overwritten
                          alreadyfilled <- samtofill[!is.na(PopInfo(objectm)[samtofill])]
                          if(owerror && !identical(PopInfo(objectm)[alreadyfilled],
                                                   PopInfo(objectT)[alreadyfilled])) stop(
                            "Conflicting data in PopInfo.  Set overwrite to \"x\" or \"y\".")
                          # fill in the population number
                          PopInfo(objectm)[samtofill] <-
                              PopNum(objectm, p)
                      }
                  }
              }
              #
              # easiest scenario- PopNames and PopInfo numbers are the same, do like above
              # more complicated- different sets of populations, but numbers are consistent
              # most complicated- numbers have different meanings in the two objects

              # return the merged object.
              return(objectm)
          })

#### genambig methods
# initialization for a genambig object
setMethod("initialize", signature(.Object = "genambig"),
          function(.Object, samples, loci, Missing, ...){
              # check to see that arguments are all there
              if(missing(samples)) samples <- c("ind1","ind2")
              if(missing(loci)) loci <- c("loc1","loc2")
              if(missing(Missing)) Missing <- as.integer(-9)

              # add empty genotype list
    .Object@Genotypes <- array(list(Missing), dim=c(length(samples),
                                                length(loci)),
                                 dimnames=list(samples, loci))

              callNextMethod(.Object, samples = samples, loci = loci, Missing=Missing,...)
              })

# replacement method for sample names
setReplaceMethod("Samples", "genambig", function(object, value){
    dimnames(object@Genotypes)[[1]] <- value
    callNextMethod(object, value)
})

# replacement method for locus names
setReplaceMethod("Loci", "genambig", function(object, value){
    dimnames(object@Genotypes)[[2]] <- value
    callNextMethod(object, value)
})

## Methods to get and edit genotype data
setMethod("Genotype", signature(object = "genambig", sample = "ANY",
                                locus = "ANY"),
          function(object, sample, locus){
              if(length(sample) != 1 | length(locus) != 1){
                  stop("sample and locus arguments must have only one element.")
              }
              return(object@Genotypes[[sample, locus]])
          })
setReplaceMethod("Genotype", "genambig", function(object, sample, locus, value){
    if(length(sample) != 1 | length(locus) != 1){
        stop("sample and locus arguments must have only one element.")
    }
    if(!is.vector(value)) stop("Assigned value must be a vector.")
    object@Genotypes[[sample, locus]] <- value
    object
})
setMethod("Genotypes", signature(object = "genambig", samples = "ANY",
                                 loci = "ANY"),
          function(object, samples, loci){
              return(object@Genotypes[samples, loci, drop=FALSE])
          })
setReplaceMethod("Genotypes", "genambig",
                 function(object, samples,
                          loci, value){
                     if(length(dim(value))> 2)
                         stop("Assigned value has too many dimensions.")
                     if(FALSE %in% mapply(is.vector, value))
                         stop("All array elements must be vectors.")
    object@Genotypes[samples, loci] <- value
    object
})

## Method to determine if a genotype is missing
setMethod("isMissing", "genambig", function(object, samples, loci){
    if(length(samples) == 1 && length(loci) == 1){
        return(Genotype(object, samples, loci)[1] == Missing(object))
    } else {
        result <- array(FALSE, dim=c(length(samples), length(loci)),
                        dimnames=list(samples, loci))
        for(s in samples){
            for(L in loci){
                result[s,L] <- isMissing(object, s, L)
            }
        }
        return(result)
    }
})

# summary method for genambig
setMethod("summary", "genambig", function(object){
    nummissing <- 0
    for(L in Loci(object)){
        for(s in Samples(object)){
            if(isMissing(object, s, L)) nummissing <- nummissing + 1
        }
    }
    cat("Dataset with allele copy number ambiguity.",
        Description(object),
        paste("Number of missing genotypes:", nummissing), sep="\n")
    callNextMethod(object)
})

# Method to print the genotypes to the console
setMethod("viewGenotypes", signature(object = "genambig", samples = "ANY",
                                     loci = "ANY"),
          function(object, samples, loci){
                  # print a header
                  cat("Sample\tLocus\tAlleles", sep="\n")

                  # loop through genotypes and print them
                  for(L in loci){
                      for(s in samples){
                          cat(s, L, Genotype(object, s, L), sep="\t")
                          cat("\n", sep="")
                      }
                  }
          })

# method for deleting samples
setMethod("deleteSamples", "genambig", function(object, samples){
    samtouse <- Samples(object)[!Samples(object) %in% samples]
    object@Genotypes <- as.array(object@Genotypes[samtouse, Loci(object)])
    callNextMethod(object, samples)
})

# method for deleting loci
setMethod("deleteLoci", "genambig", function(object, loci){
    loctouse <- Loci(object)[!Loci(object) %in% loci]
    object@Genotypes <- as.array(object@Genotypes[Samples(object), loctouse])
    callNextMethod(object, loci)
})

# subscripting method
setMethod("[", signature(x = "genambig", i="ANY", j="ANY"), function(x, i, j){
    x@Genotypes <- x@Genotypes[i,j, drop=FALSE]
    x@PopInfo <- x@PopInfo[i]
    if(is(x@Ploidies, "ploidymatrix")){
      x@Ploidies@pld <- x@Ploidies@pld[i,j, drop=FALSE]
    }
    if(is(x@Ploidies, "ploidysample")){
      x@Ploidies@pld <- x@Ploidies@pld[i]
    }
    if(is(x@Ploidies, "ploidylocus")){
      x@Ploidies@pld <- x@Ploidies@pld[j]
    }
    x@Usatnts <- x@Usatnts[j]
    return(x)
})

# method to change code for missing genotypes
setReplaceMethod("Missing", "genambig", function(object, value){
    if(length(value) != 1) stop("Assigned value should have only one element.")
    for(L in Loci(object)){
        for(s in Samples(object)){
            if(isMissing(object, s, L)) Genotype(object, s, L) <- value
        }
    }
    callNextMethod(object, value)
})

# Method to edit genotypes in a spreadsheet-like format
setMethod("editGenotypes", "genambig",
          function(object, maxalleles, samples, loci){
              if(is.na(maxalleles))
                  stop("Enter ploidies first or give argument for maxalleles.")
              # take the selected genotypes and convert to a data frame
              dummyarray <- matrix(NA, nrow=length(samples)*length(loci), ncol=maxalleles,
                                   dimnames=list(NULL,
                                   paste("Allele", 1:maxalleles, sep="")))
              dummyrow <- 1
              for(L in loci){
                  for(s in samples){
                      numalleles <- length(Genotype(object, s, L))
                      if(numalleles > maxalleles)
                          stop(paste("Set maxalleles to at least", numalleles))
                      dummyarray[dummyrow, 1:numalleles] <- Genotype(object, s, L)
                      dummyrow <- dummyrow + 1
                  }
              }
              samvect <- rep(samples, times=length(loci))
              locvect <- rep(loci, each=length(samples))
              genframe <- data.frame(Samples=samvect, Loci=locvect, dummyarray,
                                     stringsAsFactors=FALSE)
              # open data frame for editing by user
              cat("Edit the alleles, then close the data editor window.", sep="\n")
              genframe <- edit(genframe)
              # overwrite genotypes from those in data frame
              for(r in 1:dim(genframe)[1]){
                  thesealleles <- genframe[r, 3:dim(genframe)[2]]
                  thesealleles <- unique(thesealleles[!is.na(thesealleles)])
                  Genotype(object, genframe[r, "Samples"],
                           genframe[r, "Loci"]) <- thesealleles
              }

              # return object
              return(object)
          })

# Generic function and method to estimate ploidy
setMethod("estimatePloidy", "genambig",
          function(object, extrainfo, samples, loci){
            # get Ploidies into ploidysample format if necessary
              if(!is(object@Ploidies, "ploidysample")){
                object <- reformatPloidies(object, output="sample", erase=TRUE)
              }
              # set up array to contain the maximum and average number of alleles
              ploidyinfo <- array(dim=c(length(samples),2),
                                  dimnames=list(samples,
                                  c("max.alleles","mean.alleles")))
              # fill the array
              for(s in samples){
                  genotypes <- Genotypes(object, s, loci)[!isMissing(object, s, loci)]
                  if(length(genotypes) ==0){
                      ploidyinfo[s,] <- c(NA, NA)
                  } else {
                      numalleles <- mapply(length, genotypes)
                      ploidyinfo[s, "max.alleles"] <- max(numalleles)
                      ploidyinfo[s, "mean.alleles"] <- mean(numalleles)
                  }
              }

              ## Build data frame
              ploidyinfo <- as.data.frame(ploidyinfo)
              # Add in extrainfo, allowing for named or unnamed vectors, arrays, and
              # data frames.
              if(!missing(extrainfo)){
                  if(is.vector(extrainfo)){
                      if(!is.null(names(extrainfo))) extrainfo <- extrainfo[samples]
                      ploidyinfo <- data.frame(extrainfo = extrainfo, ploidyinfo)
                  } else {
                      if(is.array(extrainfo) | is.data.frame(extrainfo)){
                          extrainfo <- as.data.frame(extrainfo)
                          if(!identical(row.names(extrainfo),
                                        as.character(1:dim(extrainfo)[1])))
                              extrainfo <- extrainfo[samples,]
                          ploidyinfo <- data.frame(extrainfo, ploidyinfo)
                      }
                  }
              }
              # Add in PopInfo from the object
              if(length(unique(PopInfo(object))) >1)
                  ploidyinfo <- data.frame(population = PopInfo(object)[samples],
                                           ploidyinfo)
              # Add in Ploidies from the object
              if(length(unique(Ploidies(object))) >1)
                  ploidyinfo <- data.frame(ploidyinfo,
                                           current.ploidy = Ploidies(object)[samples])
              # Make a column for new ploidies
              ploidyinfo <- data.frame(ploidyinfo, new.ploidy = ploidyinfo$max.alleles)

              # let user edit data frame
              cat("Edit the new.ploidy values, then close the data editor window.",
                  sep="\n")
              ploidyinfo <- edit(ploidyinfo)

              # write new.ploidy to object@Ploidies
              Ploidies(object)[samples] <- ploidyinfo$new.ploidy

              # return the object with edited ploidies
              return(object)
          })

setMethod("merge", signature(x="genambig", y="genambig"),
          function(x, y, objectm, samples, loci, overwrite){
              # set up new genambig object if this wasn't called from the method
              # of a subclass.
              if(missing(objectm)){
                  if(missing(samples)) samples <- unique(c(Samples(x), Samples(y)))
                  if(missing(loci)) loci <- unique(c(Loci(x), Loci(y)))
                  objectm <- new("genambig", samples=samples, loci=loci)
              }

              # determine what to do with conflicting data: overwrite or give error
              if(missing(overwrite)) overwrite <- "empty"
              owerror <- ifelse(overwrite != "x" && overwrite != "y", TRUE, FALSE)
              if(overwrite=="x"){
                  objectB <- x
                  objectT <- y
              } else {
                  objectB <- y
                  objectT <- x
              }


              # make sure missing data values are the same, change if appropriate
              if(Missing(x) != Missing(y)){
                  if(owerror) stop(
                   "Different missing data symbols used; set overwrite to \"x\" or \"y\".")
                  Missing(objectB) <- Missing(objectT)
              }
              Missing(objectm) <- Missing(objectT)

              # merge genotypes into objectm
              # put in objectB first
              samtofillB <- Samples(objectm)[Samples(objectm) %in% Samples(objectB)]
              loctofillB <- Loci(objectm)[Loci(objectm) %in% Loci(objectB)]
              Genotypes(objectm, samtofillB, loctofillB) <-
                  Genotypes(objectB, samtofillB, loctofillB)
              # then put in objectT
              samtofillT <- Samples(objectm)[Samples(objectm) %in% Samples(objectT)]
              loctofillT <- Loci(objectm)[Loci(objectm) %in% Loci(objectT)]
              Genotypes(objectm, samtofillT, loctofillT) <-
                  Genotypes(objectT, samtofillT, loctofillT)
              # check to see if objectB overwritten if there wasn't supposed to be
              # conflicting data
              if(owerror && !identical(Genotypes(objectm, samtofillB, loctofillB),
                                       Genotypes(objectB, samtofillB, loctofillB)))
                 stop("Conflicting genotype data.  Set overwrite to \"x\" or \"y\".")

              # call gendata method to merge popinfo etc
              callNextMethod(x, y, objectm, overwrite=overwrite)
          })

#### genbinary methods
setMethod("initialize", signature(.Object = "genbinary"),
          function(.Object, samples, loci, ...){
              # fill in empty arguments if necessary
              if(missing(samples)) samples <- c("ind1", "ind2")
              if(missing(loci)) loci <- c("loc1", "loc2")
#              if(missing(Missing)) Missing <- as.integer(-9)
#              if(missing(Present)) Present <- as.integer(1)
#              if(missing(Absent)) Absent <- as.integer(0)

              # add empty genotype matrix
              .Object@Genotypes <- matrix(nrow=length(samples), ncol=0,
                                          dimnames=list(samples, NULL))
              # go to the gendata method
              callNextMethod(.Object, samples=samples, loci=loci, ...)
          })

# access and replacement methods for Present and Absent values
setMethod("Present", "genbinary", function(object) return(object@Present))
setReplaceMethod("Present", signature(object="genbinary"), function(object, value){
    # potential errors
    if(length(value) != 1) stop("Value must be vector of length 1.")
    if(value == object@Absent) stop("Value must be different from Absent.")
    if(value == object@Missing) stop("Value must be different from Missing.")

    # replace the value in the genotype matrix
    object@Genotypes[object@Genotypes == object@Present] <- value

    # replace the value in the Present slot
    object@Present <- value

    return(object)
})
setMethod("Absent", "genbinary", function(object) return(object@Absent))
setReplaceMethod("Absent", signature(object="genbinary"), function(object, value){
    # potential errors
    if(length(value) != 1) stop("Value must be vector of length 1")
    if(value == object@Present) stop("Value must be different from Present.")
    if(value == object@Missing) stop("Value must be different from Missing.")

    # replace the value in the genotype matrix
    object@Genotypes[object@Genotypes == object@Absent] <- value

    # replace the value in the Absent slot
    object@Absent <- value

    return(object)
})

# replacement method for Missing
setReplaceMethod("Missing", signature(object="genbinary"), function(object, value){
    # potential errors
    if(length(value) != 1) stop("Value must be vector of length 1.")
    if(value == object@Present) stop("Value must be different from Present")
    if(value == object@Absent) stop("Value must be different from Absent")

    # replace the value in the genotype matrix
    object@Genotypes[object@Genotypes == object@Missing] <- value

    # go to the gendata method to fill in Missing slot
    callNextMethod(object, value)
})

# methods for retrieving and replacing genotypes
setMethod("Genotypes", "genbinary", function(object, samples, loci){
    # set up vector to index columns for these loci
    loccolumns <- integer(0)
    for(L in loci){
        loccolumns <- c(loccolumns, grep(paste(L,".",sep=""),
                                         dimnames(object@Genotypes)[[2]],
                                         fixed = TRUE))
    }

    # return subset of the matrix
    return(object@Genotypes[samples, loccolumns, drop=FALSE])
})
setMethod("Genotype", "genbinary", function(object, sample, locus){
    return(Genotypes(object, sample, locus))
})

setReplaceMethod("Genotypes", "genbinary",
                 function(object, samples, loci, value){
                     # make data frame of loci and alleles by column
                     tempcol <- strsplit(dimnames(value)[[2]], split=".", fixed=TRUE)
                     colinfo <- data.frame(Locus=rep("blank", length(tempcol)),
                                           Allele=rep("blank", length(tempcol)),
                                           stringsAsFactors=FALSE)
                     for(i in 1:length(tempcol)){
                         if(length(tempcol[[i]]) != 2)
                             stop("Column names should be of form \"locus.allele\".")
                         colinfo[i,]<-tempcol[[i]]
                     }
                     # check that loci are in the object
                     if(FALSE %in% (unique(colinfo$Locus) %in% loci))
                         stop("Locus names in column names should match loci in object.")
                     # check that length of first dimension matches # samples
                     if(dim(value)[1] != length(samples))
                         stop("Number of rows should match number of samples to replace.")
                     # check that only valid symbols are used
                     if(FALSE %in% (unique(value) %in% c(Absent(object), Present(object),
                                                         Missing(object))))
                         stop("Symbols in matrix must match Absent(object), Present(object), or Missing(object).")

                     # get samples that aren't being filled in (for indexing blank
                     # spaces below)
                     samblank <- Samples(object)[!Samples(object) %in% samples]

                     # fill in the values
                     for(i in 1:dim(value)[2]){
                         # get the name of this column
                         colname <- dimnames(value)[[2]][i]
                     # for allele columns that already exist, replace
                         if(colname %in% dimnames(object@Genotypes)[[2]]){
                             object@Genotypes[samples, colname] <-
                                 value[,i]
                         } else {
                     # for allele columns that don't exist, add them
                             object@Genotypes <- cbind(object@Genotypes,
                                  matrix(nrow=length(Samples(object)), ncol=1,
                                         dimnames=list(NULL, colname)))
                             object@Genotypes[samples, colname] <-
                                 value[,i]
                     # fill any new blank spaces with absent or missing
                             for(s in samblank){
                                 g <- Genotype(object, s, colinfo[i,"Locus"])
                                 if(all(is.na(g)) || all(g[!is.na(g)]==Missing(object))){
                                     object@Genotypes[s,colname] <- Missing(object)
                                 } else {
                                     object@Genotypes[s,colname] <- Absent(object)
                                 }
                             }
                         }
                     }

                     return(object)
                 })

# replacement methods for samples and loci
setReplaceMethod("Samples", "genbinary", function(object, value){
    # change the names in the Genotypes slot
    dimnames(object@Genotypes)[[1]] <- value

    # go to the gendata method to change sample names in other slots.
    callNextMethod(object, value)
})
setReplaceMethod("Loci", "genbinary", function(object, value){
    # get the original locus names and systematically replace them
    oldloci <- Loci(object)
    for(i in 1:length(oldloci)){
        dimnames(object@Genotypes)[[2]] <- gsub(paste(oldloci[i],".",sep=""),
                                                paste(value[i],".",sep=""),
                                                dimnames(object@Genotypes)[[2]],
                                                fixed=TRUE)
    }
    # go to gendata method to change locus names in Usatnts
    callNextMethod(object, value)
})

# method to determine if a genotype is missing
setMethod("isMissing", "genbinary", function(object, samples, loci){
    if(length(samples) == 1 && length(loci) == 1){
        return(TRUE %in% (Genotypes(object, samples, loci) == Missing(object)))
    } else {
        result <- array(FALSE, dim=c(length(samples), length(loci)),
                        dimnames=list(samples, loci))
        for(s in samples){
            for(L in loci){
                result[s,L] <- isMissing(object, s, L)
            }
        }
        return(result)
    }
})

# summary method for genbinary
setMethod("summary", "genbinary", function(object){
    # tally missing genotypes
    nummissing <- 0
    for(L in Loci(object)){
        for(s in Samples(object)){
            if(isMissing(object, s, L)) nummissing <- nummissing + 1
        }
    }
    # print stuff about dataset
    cat("Binary presence/absence dataset.",
        Description(object),
        paste("Number of missing genotypes:", nummissing), sep="\n")
    # go to the summary method for gendata
    callNextMethod(object)
})

# method for viewing genotypes
setMethod("viewGenotypes", "genbinary", function(object, samples, loci){
    for(L in loci){
        print(Genotypes(object, samples, L))
        cat("\n", sep="")
    }
})

# method for editing genotypes
setMethod("editGenotypes", "genbinary", function(object, maxalleles, samples, loci){
    Genotypes(object, samples) <- edit(Genotypes(object, samples, loci))
    return(object)
})

# methods for deleting samples and loci
setMethod("deleteSamples", "genbinary", function(object, samples){
    samtouse <- Samples(object)[!Samples(object) %in% samples]
    object@Genotypes <- object@Genotypes[samtouse,]
    callNextMethod(object, samples)
})
setMethod("deleteLoci", "genbinary", function(object, loci){
    loccolumns <- integer(0)
    for(L in loci){
        loccolumns <- c(loccolumns, grep(paste(L,".",sep=""),
                                         dimnames(object@Genotypes)[[2]],
                                         fixed = TRUE))
    }
    object@Genotypes <- object@Genotypes[,-loccolumns]
    callNextMethod(object, loci)
})

# subscripting methods
setMethod("[", signature(x = "genbinary", i = "ANY", j= "ANY"), function(x, i, j){
    x@Genotypes <- Genotypes(x, i, j)
    x@PopInfo <- x@PopInfo[i]
    if(is(x@Ploidies, "ploidymatrix")){
      x@Ploidies@pld <- x@Ploidies@pld[i,j, drop=FALSE]
    }
    if(is(x@Ploidies, "ploidysample")){
      x@Ploidies@pld <- x@Ploidies@pld[i]
    }
    if(is(x@Ploidies, "ploidylocus")){
      x@Ploidies@pld <- x@Ploidies@pld[j]
    }
    x@Usatnts <- x@Usatnts[j]
    return(x)
})

# method for estimating ploidy
setMethod("estimatePloidy", "genbinary",
          function(object, extrainfo, samples, loci){
            # get Ploidies into ploidysample format if necessary
              if(!is(object@Ploidies, "ploidysample")){
                object <- reformatPloidies(object, output="sample", erase=TRUE)
              }
              # set up array to contain the maximum and average number of alleles
              ploidyinfo <- array(dim=c(length(samples),2),
                                  dimnames=list(samples,
                                  c("max.alleles","mean.alleles")))
              # Convert Present and Absent to 1's and 0's so
              # simple addition can be done
              objectA <- object
              Missing(objectA) <- as.integer(-9)
              Absent(objectA) <- as.integer(0)
              Present(objectA) <- as.integer(1)

              # fill the array
              for(s in samples){
                  loctouse <- loci[!isMissing(objectA, s, loci)]
                  if(length(loctouse) ==0){
                      ploidyinfo[s,] <- c(NA, NA)
                  } else {
                      numalleles <- rep(0, times=length(loctouse))
                      names(numalleles) <- loctouse
                      for(L in loctouse){
                          numalleles[L] <- sum(Genotypes(objectA,s,L))
                      }
                      ploidyinfo[s, "max.alleles"] <- max(numalleles)
                      ploidyinfo[s, "mean.alleles"] <- mean(numalleles)
                  }
              }

              ## Build data frame
              ploidyinfo <- as.data.frame(ploidyinfo)
              # Add in extrainfo, allowing for named or unnamed vectors, arrays, and
              # data frames.
              if(!missing(extrainfo)){
                  if(is.vector(extrainfo)){
                      if(!is.null(names(extrainfo))) extrainfo <- extrainfo[samples]
                      ploidyinfo <- data.frame(extrainfo = extrainfo, ploidyinfo)
                  } else {
                      if(is.array(extrainfo) | is.data.frame(extrainfo)){
                          extrainfo <- as.data.frame(extrainfo)
                          if(!identical(row.names(extrainfo),
                                        as.character(1:dim(extrainfo)[1])))
                              extrainfo <- extrainfo[samples,]
                          ploidyinfo <- data.frame(extrainfo, ploidyinfo)
                      }
                  }
              }
              # Add in PopInfo from the object
              if(length(unique(PopInfo(object))) >1)
                  ploidyinfo <- data.frame(population = PopInfo(object)[samples],
                                           ploidyinfo)
              # Add in Ploidies from the object
              if(length(unique(Ploidies(object))) >1)
                  ploidyinfo <- data.frame(ploidyinfo,
                                           current.ploidy = Ploidies(object)[samples])
              # Make a column for new ploidies
              ploidyinfo <- data.frame(ploidyinfo, new.ploidy = ploidyinfo$max.alleles)

              # let user edit data frame
              cat("Edit the new.ploidy values, then close the data editor window.",
                  sep="\n")
              ploidyinfo <- edit(ploidyinfo)

              # write new.ploidy to object@Ploidies
              Ploidies(object)[samples] <- ploidyinfo$new.ploidy

              # return the object with edited ploidies
              return(object)
          })

# method for merging objects
setMethod("merge", signature(x = "genbinary", y="genbinary"),
          function(x, y, objectm, samples, loci, overwrite){
    # set up new genbinary object if this wasn't called from the method of
    # a subclass
    if(missing(objectm)){
        if(missing(samples)) samples <- unique(c(Samples(x), Samples(y)))
        if(missing(loci)) loci <- unique(c(Loci(x), Loci(y)))
        objectm <- new("genbinary", samples, loci)
    }

    # determine what to do with conflicting data
    if(missing(overwrite)) overwrite <- "empty"
    owerror <- (!overwrite %in% c("x", "y"))
    if(overwrite=="x"){
        objectB <- x
        objectT <- y
    } else {
        objectB <- y
        objectT <- x
    }

    # make sure Missing, Present, and Absent are the same
    if(Missing(x) != Missing(y)){
        if(owerror) stop("Different missing data symbols used; set overwrite to \"x\" or \"y\".")
        Missing(objectB) <- Missing(objectT)
    }
    Missing(objectm) <- Missing(objectT)
    if(Present(x) != Present(y)){
        if(owerror) stop("Different values for Present; set overwrite to \"x\" or \"y\".")
        Present(objectB) <- Present(objectT)
    }
    Present(objectm) <- Present(objectT)
    if(Absent(x) != Absent(y)){
        if(owerror) stop("Different values for Absent; set overwrite to \"x\" or \"y\".")
        Absent(objectB) <- Absent(objectT)
    }
    Absent(objectm) <- Absent(objectT)

    # merge genotypes into objectm
    # put in objectB first
    samtofillB <- Samples(objectm)[Samples(objectm) %in% Samples(objectB)]
    loctofillB <- Loci(objectm)[Loci(objectm) %in% Loci(objectB)]
    Genotypes(objectm, samtofillB, loctofillB) <-
    Genotypes(objectB, samtofillB, loctofillB)
    # then put in objectT
    samtofillT <- Samples(objectm)[Samples(objectm) %in% Samples(objectT)]
    loctofillT <- Loci(objectm)[Loci(objectm) %in% Loci(objectT)]
    Genotypes(objectm, samtofillT, loctofillT) <-
    Genotypes(objectT, samtofillT, loctofillT)
    # some missing data symbols may have been introduced; get rid of them.
    # (this can occur if B has alleles that T doesn't)
    for(s in samtofillT){
        for(L in loctofillT){
            if(isMissing(objectm, s, L) && !isMissing(objectT, s, L)){
                Ballele <- names(Genotypes(objectm,s,L))[Genotypes(objectm,s,L)
                                                         == Missing(objectm)]
                objectm@Genotypes[s,Ballele] <- Absent(objectm)
            }
        }
    }
    # check to see if objectB overwritten if there wasn't supposed to be
    # conflicting data
    if(owerror && !identical(Genotypes(objectm, samtofillB, loctofillB),
                             Genotypes(objectB, samtofillB, loctofillB)))
       stop("Conflicting genotype data.  Set overwrite to \"x\" or \"y\".")

    # call gendata method to merge popinfo etc
    callNextMethod(x, y, objectm, overwrite=overwrite)
})

#}
