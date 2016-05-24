
## Functions: pedgene 
## wrapper for computing retrospective likelihood stat on pedigrees for rare
## variants over multiple genes
## Authors: Jason Sinnwell and Dan Schaid

pedgene <- function(ped, geno, map=NULL, male.dose=2, checkpeds=TRUE, verbose.return=FALSE,
                    weights=NULL, weights.beta=c(1,25), weights.mb=FALSE,
                    method="kounen", acc.davies=1e-5) {

##Arguments:
##  
##  Required:
##  ped: data.frame with columns needed to create the pedigree
##  geno: data.frame with ped, person ids in the first two columns, and numeric columns
##       with minor allele count (0/1/2) for markers (columns) and subjects (rows) in
##       the ped object
##  
##  Optional:
##  map: data.frame with columns chrom, position, and gene.
##        gene can be gene symbol or geneid (entrez or ensemble);
##        it is not matched with other annotation, just used for marker groups
##        If not passed, assume all variants from the same gene 
##  male.dose: When doing X-chrom, define how male genotypes should be
##        analyzed. male.dose can be between 0 and 2, but usually 1 or 2
##  checkpeds, perform basic pedigree structure checks. Time-consuming if peds are
##             are already validated
##  method: either Kounen or Davies method to calculate kernel test p-values
##  verbose.return: similar to glm in R, return geno and response data that is used in calculations
##  weights: allow user-specified weights. If Null, use beta distribution weighting
##  weights.mb: if weights=NULL and weights.MB=TRUE, do Madsen-Browning weights
##  weights.beta: beta distribution coefficients used for weighting.
##     By default, beta weights are used
##  acc.davies: numerical accuracy for davies method to choose eigven values
##       and to determine p-value
  
## Steps
## 1) verify ped columns, map and geno match dimensions
## 2) Create kinship matrices for autosomes and X for all subjects
## 3) run pedgene.stats on each gene given in map

  verbose=FALSE
  
  ## save the call
  call <- match.call() 
  ## save options before setting stringsAsFactors for just this function
  saveOpt <- options()
  options(stringsAsFactors=FALSE)

  ## check method, must be davies or kounen
  method=casefold(method[1])
  method=c("kounen","davies")[pmatch(method, c("kounen", "davies"))]
  if(is.na(method)) {
    warning("method not kounen or davies, setting to kounen\n")
    method="kounen"
  }  
  
  ## require kinship function to be recent
  kin2v <- sessionInfo()$otherPkgs$kinship2$Version
  if(is.null(kin2v)) {
    kin2v <- sessionInfo()$loadedOnly$kinship2$Version
  }
  if(as.numeric(substring(kin2v, 1, nchar(kin2v)-2)) < 1.5) {
    stop("kinship2 needs to be version 1.5.3 or later\n")
  }
  
  ## if no map, create one, assuming all one gene
  if(is.null(map)) {
    map <- data.frame(chrom=rep("unknown", ncol(geno)-2),
                      gene=rep("unknown", ncol(geno)-2))
  }

  ## unify names of ped and map to lowercase
  names(map) <- casefold(names(map))
  names(ped) <- casefold(names(ped))

  ## verify map data.frame ###
  ## i. check column names
  if(any(!(c("chrom", "gene") %in% names(map)))) {
    stop("map requires columns for chrom and gene")
  }
  ## ii. recode chrom 23 to X
  map$chrom[map$chrom==23] <- "X"

  ## verify geno matrix, and that it matches map
  if(any(!(c("ped", "person") %in% names(geno)))) {
    stop("geno requires columns 'ped' and 'person' ids")
  }
  
  ## get indices of ped/person of geno to match ped, then strip off those columns
  keepped <- match(paste(geno$ped, geno$person, sep="-"),
                   paste(ped$ped, ped$person, sep="-"))
  tblkeep <- table(keepped)

  ## check for multiple subject entries and not matching a subject
  if(any(tblkeep > 1)) {
    warning(paste("subject with multiple entries, only the first is used: ", names(tblkeep)[which(tblkeep>1)], ".\n", sep=""))
    geno <- geno[!duplicated(paste(geno$ped, geno$person, sep="-")),,drop=FALSE]
    keepped <- match(paste(geno$ped, geno$person, sep="-"),
                     paste(ped$ped, ped$person, sep="-"))
  }

  if(any(is.na(keepped))) { 
    warning("removing subject in genotype matrix who is not in pedigree \n")
    geno <- geno[!is.na(keepped),,drop=FALSE]
    keepped <- match(paste(geno$ped, geno$person, sep="-"),
                     paste(ped$ped, ped$person, sep="-"))
  }
  
  ## after matching, get rid of id cols  
  geno <- geno[,!(names(geno) %in% c("ped", "person")),drop=FALSE]

  ## rm subjects in geno who are not in ped  
  
  if(nrow(map) != (ncol(geno))) {
    stop(paste("map rows (", nrow(map), ") and geno columns (", ncol(geno),
               ") do not match \n",sep=""))
  }  

  ## Check that geno for males on X should only have 0 and 1 dosages
  xidx <- which(map$chrom=="X" | map$chrom=="x")
  if(length(xidx)) {
    xdosemale <- geno[ped$sex[keepped]==1,xidx, drop=TRUE]
    if(sum(xdosemale>1, na.rm=TRUE)) {
      stop("All male dose on X chromosome should be <= 1")
    }
  }
  
  ## verify ped data.frame has expected column names
  if(any(!(c("ped", "person", "father", "mother", "sex", "trait")
           %in% names(ped)))) {
    stop("Error: ped requires columns: ped, person, father, mother, sex, trait")
  }
 
  #############################################################################
  ## this is where to do trait.adjusted if we want it on all people that have trait 
  ## this caused different results when flipping major/minor alleles, so moved later
  #############################################################################
  
  ## check weights parameters
  ## verify user-passed weights, match ncol(geno)
  if(!is.null(weights)) {
    ## by default, do Beta weights, implemented in ped.gene.stats
    ## otherwise, these are user-specified, check length
    if(length(weights) != ncol(geno)) {
       stop(paste("Error: should have weights(", length(weights),
                  ") for every variant position(", ncol(geno), ")", sep=""))
    }
  } else {  ## no user-given weights
    if(weights.mb==FALSE) {
      ## verify weights.beta
      if(length(weights.beta) != 2 | any(weights.beta < 0)) {
        warning("weights.beta should be two positive numbers, setting to (1,25)\n")
        weights.beta=c(1,25)
      }
    }  ## m-b weights, nothing to check except that weights.mb is true/false
  }
  
  ## perform simple pedigree checks
  if(checkpeds) {
    uped <- unique(ped$ped)
    nped <- length(uped)
    
    for(i in 1:nped) {      
      iped <- uped[i]      
      temp.ped <- ped[ped$ped == iped,, drop=FALSE]      
      if(nrow(temp.ped) > 1) {      
        ## simple checks on pedigree
        pedigreeChecks(temp.ped, male.code=1, female.code=2)
      }
    }
  }
  ## additional checks <could> be done on peds when creating pedlist object,
  ## which could be used to create kinmat.
  # pedall <- with(ped, kinship2::pedigree(id=person, dadid=father, momid=mother,
  #                          sex=sex, famid=ped, missid=missid))
  #  kinmat <- kinship2::kinship(pedall, chrtype="auto")

  ## We rather created it directly from ped  
  ## create kinship matrix, also for X if any genes on X chrom
  ## subset to only those with geno rows
  kinmat <- Matrix(with(ped, kinship(id=paste(ped,person,sep="-"),
                      dadid=ifelse(father>0,paste(ped,father,sep="-") , as.character(father)),
                      momid=ifelse(mother>0, paste(ped,mother,sep="-"), as.character(mother)),
                      sex=sex, chrtype="autosome")))
  kinmat <- kinmat[keepped, keepped]
 
  if(any(map$chrom=="X")) {
    kinmatX <- Matrix(with(ped, kinship(id=paste(ped,person,sep="-"),
                      dadid=ifelse(father>0,paste(ped,father,sep="-") , as.character(father)),
                      momid=ifelse(mother>0, paste(ped,mother,sep="-"), as.character(mother)),
                      sex=sex, chrtype="X")))
    kinmatX <- kinmatX[keepped, keepped]
  } else {
    kinmatX <- NULL
  }
  ped <- ped[keepped,]
 
  
  ## subset pedgeno kinmat, kinmatX to only subject who have genotype data
  missidx <- is.na(ped$trait) | apply(is.na(geno), 1, all) 
  if("trait.adjusted" %in% names(ped)) missidx <- missidx | is.na(ped$trait.adjusted)
  if(sum(missidx)>0) {
    ped <- ped[!missidx,]
    kinmat <- kinmat[!missidx, !missidx]
    kinmatX <- kinmatX[!missidx, !missidx]
    geno <- geno[!missidx,,drop=FALSE]
  }

## this is where trait.adjusted should be calculated if its on everyone with genotype data
  ## add trait.adjusted if not already there
  if(!("trait.adjusted" %in% names(ped))) {
    ped$trait.adjusted <- mean(ped$trait, na.rm=TRUE)      
  }
  
  gvec <- chromvec <- nvariant <- noninform <- kstat <- kpval <- bstat <- bpval <- NULL
  
  for(g in unique(map$gene)) {
    if(verbose) {
      cat("test on gene ", g, "\n")
    }
    gidx <- which(map$gene==g)
    ## drop=FALSE for 1-marker gene
    genosub <- geno[,gidx,drop=FALSE]

    resid <- ped$trait - ped$trait.adjusted
    sex <- ped$sex
    chrom <- map$chrom[gidx[1]]
    
    c.factor <- quadfactor(
             if(chrom=="X") kinmatX else kinmat,
             chrom, resid, sex, male.dose)

    pgstat <- pedgene.stats(genosub, as.vector(c.factor), map$chrom[gidx[1]], male.dose, sex, resid,
                    weights=weights[gidx], weights.beta=weights.beta, weights.mb=weights.mb,
                    method=method, acc.davies=acc.davies)
    if(pgstat$nvariant==0) {
      cat("gene: '", g, "' has no markers after removing markers with all same genotype\n")     
    }
    gvec <- c(gvec, g)
    chromvec <- c(chromvec, chrom)
    nvariant <- c(nvariant,pgstat$nvariant)
    noninform <- c(noninform, pgstat$noninform)
    kstat <- c(kstat, pgstat$stat.kernel)
    kpval <- c(kpval, pgstat$pval.kernel)
    bstat <- c(bstat, pgstat$stat.burden)
    bpval <- c(bpval, pgstat$pval.burden)   
  }
  
  pgdf <- data.frame(gene=gvec, chrom=chromvec, n.variant=nvariant,
                     n.noninform=noninform,
                     stat.kernel=kstat, pval.kernel=kpval,
                     stat.burden=bstat, pval.burden=bpval)
  
  # re-set options
  options(saveOpt)
  if(verbose.return) {
    save <- list(geno=geno, ped=ped, map=map)
  } else {
    save=NULL
  }
  pglist <- list(pgdf=pgdf, call=call, save=save)
  class(pglist) <- "pedgene"
  return(pglist)
}

## print and summary methods for pedgene S3 class

print.pedgene <- function(x, ...) {
## suggest digits=4  
  print.data.frame(x$pgdf, ...)

  invisible()
}
summary.pedgene <- function(object, ...) {
## suggest digits=4 or 5  
  cat("\nSummary for pedgene object: \n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n\n")

  ## invoke print method
  print(object,  ...)
  
  invisible()
}
