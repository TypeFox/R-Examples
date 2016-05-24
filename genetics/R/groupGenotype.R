### groupGenotype.R
###------------------------------------------------------------------------
### What: Group genotype values code
### Time-stamp: <2007-07-21 12:04:16 ggorjan>
###------------------------------------------------------------------------

groupGenotype <- function(x, map, haplotype=FALSE, factor=TRUE,
                          levels=NULL, verbose=FALSE)
{
  if(!is.genotype(x))
    stop("'x' must be of a genotype or haplotype class")
  if(any(names(map) == ""))
    stop("all list components in 'map' need to have a name")

  alleles <- allele.names(x)

  ## Put else at the end and change it to */*
  elseF <- FALSE
  elsePos <- sapply(map, function(x) length(x) == 1 && x == ".else")
  if(any(elsePos)) {
    elseF <- TRUE
    map <- c(map[!elsePos], map[elsePos])
    map[elsePos] <- "*/*"
  }

  ## Extend the map
  for(i in seq(along=map)) {
    map[[i]] <- unlist(.matchGenotype(alleles=alleles,
                                                 pattern=map[[i]],
                                                 haplotype=haplotype),
                       use.names=FALSE)
  }

  ## Remove duplicates sequentially over all map
  nM <- length(map)
  if(nM > 1) {
    for(i in 2:nM) {
      test <- map[[i]] %in% unlist(map[1:(i - 1)], use.names=FALSE)
      map[[i]] <- map[[i]][!test]
    }
  }

  ## Show matches
  if(verbose) print(map)

  ## Group
  x <- as.factor(x)
  if(!is.null(levels)) {
    if(length(map) != length(levels))
      warning("length of 'map' and 'levels' does not match")
    map <- map[levels]
  }
  mapLevels(x) <- as.levelsMap(map)

  ## Factor?
  if(!factor) x <- as.character(x)

  ## Return
  x
}

.matchGenotype <- function(x, alleles=allele.names(x), pattern, haplotype=FALSE)
{
  ## Internal function
  ##
  ## Finds genotype matches according to given patterns out of possible
  ## genotype values that might appear in genotype data
  ##
  ## Arguments:
  ## x - genotype or haplotype
  ## alleles, character, allele names
  ## pattern - character, pattern in form of "A/A", "A/B", "A/*", "*/A" or "*/*"
  ## haplotype - logical, should order of alleles in the pattern matter
  ##
  ## Value:
  ## A list of length equal to length of pattern values. Each list
  ## component is named with pattern and has genotype values that match a
  ## pattern
  ##
  ## Details:
  ## Internally, \code{\link{genotype}} can store values as "A/B" or "B/A",
  ## so output for pattern="A/*" holds both "A/B" and "B/A", when
  ## haplotype=FALSE and there are two alleles (A and B).
  ##
  ## Example:
  ## pattern <- c("A/*", "A/B", "*/*", "B/A")
  ## .matchGenotype(alleles=c("B", "A"), pattern=pattern)
  ## $`A/*`
  ## [1] "A/B" "B/A" "A/A"
  ## $`A/B`
  ## [1] "A/B" "B/A"
  ## $`*/*`
  ## [1] "B/B" "B/A" "A/B" "A/A"
  ## $`B/A`
  ## [1] "B/A" "A/B"
  ##
  ## .matchGenotype(alleles=c("B", "A"), pattern=pattern, haplotype=TRUE)
  ## $`A/*`
  ## [1] "A/B" "A/A"
  ## $`A/B`
  ## [1] "A/B"
  ## $`*/*`
  ## [1] "B/B" "B/A" "A/B" "A/A"
  ## $`B/A`
  ## [1] "B/A"

  if(!missing(x)) {
    if(!is.genotype(x))
      stop("'x' must be of a genotype or haplotype class")
  } else {
    if(missing(alleles))
      stop("at least one of 'x' or 'alleles' must be given")
  }
  if(missing(pattern))
    stop("'pattern' must be given")

  nP <- length(pattern)
  ret <- vector(mode="list", length=nP)
  names(ret) <- pattern

  ## Change * with allele names setup
  parts <- .genotype2Allele(x=pattern)
  parts <- cbind(parts, 1:nP)
  testStar <- parts == "*"
  nA <- length(alleles)

  ## Expand A/* to A/{alleles} etc.
  for(i in 1:nrow(parts)) {
    ## Allele beside *
    whichStar <- which(!testStar[i, 1:2])
    nWS <- length(whichStar)
    if(nWS < 2) {
      if(nWS == 1) { # A/* or */A
        a <- parts[i, whichStar]
        ## Create possible genotypes
        if(whichStar == 1) {
          parts <- rbind(parts, cbind(a, alleles, i))
        } else {
          parts <- rbind(parts, cbind(alleles, a, i))
        }
      } else {       # */*
        tmp <- expectedHaplotypes(alleles=alleles)
        tmp <- .genotype2Allele(x=tmp)
        parts <- rbind(parts, cbind(tmp, i))
      }
    }
  }

  ## Remove *
  testStar <- rowSums(parts == "*") > 0
  parts <- parts[!testStar, , drop=FALSE]

  ## Order by pattern and create genotypes
  parts <- parts[order(parts[, 3, drop=FALSE]), , drop=FALSE]
  parts <- cbind(paste(parts[, 1], parts[, 2], sep="/"), parts[, 3])

  ## Fill the return
  patternId <- unique(parts[, 2])
  for(i in 1:nP) {
    ret[[i]] <- parts[parts[, 2] == patternId[i], 1]
  }

  ## For genotype treat A/* the same as */A and A/B as B/A
  if(!haplotype) ret <- lapply(ret, .genotype2Haplotype)

  ## Return
  ret
}

###------------------------------------------------------------------------
### groupGenotype.R ends here
