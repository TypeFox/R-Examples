## helper functions, not exported
check.args <- function(X,a1,a2) {
  if (!is(X, "snp.matrix") & !is(X, "SnpMatrix"))
    stop("X argument must be a snp.matrix or SnpMatrix object")
  if(any(a1=="" | is.na(a1) | a2=="" | is.na(a2) | a1==a2))
    stop("alleles must be non-missing characters, and each pair must contain different bases")
  if(length(a1)!=ncol(X) | length(a1)!=length(a2))
    stop("require ncol(X) == length(a1) == length(a2)")
}
get.eol <- function() {
  eol <- "\n"
  if(.Platform$OS.type == "windows")
    eol <- "\r\n"
  return(eol)
}
#'Fast and flexible writing of snpStats objects to flat files
#'
#'Different genetics phasing and analysis programs (beagle, mach,
#' impute, snptest, phase/fastPhase, snphap, etc) have different requirements
#' for input files.  These functions aim to make creating these files
#' from a SnpMatrix object straightfoward.
#'
#'It's written in C, so should be reasonably fast even for large datasets.
#'
#'\code{write.simple} is the most flexible function.  It should be able to
#'write most rectangular based formats.
#'
#'Additional functions are available tailored to software that require a bit
#'more than a rectangular format: \code{\link{write.beagle}},
#'\code{\link{write.impute}}, \code{\link{write.mach}}, \code{\link{write.phase}}.
#'
#'@aliases write.simple
#'@export
#'@param X SnpMatrix object
#'@param a1 vector of first allele at each SNP
#'@param a2 vector of second allele at each SNP
#'@param bp vector of base pair positions for each SNP
#'@param fsep,gsep Field and genotype separators.
#'@param nullallele Character to use for missing alleles
#'@param file Output file name.
#'@param write.header Write a header line
#'@param transpose Output SNPs as rows, samples as columns if \code{TRUE}.  The
#'default is samples as rows, SNPs as columns, as represented internally by
#'snpStats/SnpMatrix.
#'@param write.sampleid Output sample ids
#'@param num.coding Use alleles 1 and 2 instead of supplying allele vectors.
#'@return No return value, but has the side effect of writing specified output
#'files.
#'@note This has been tested with \code{SnpMatrix} objects from the package
#'\code{snpStats} but should also work with \code{snp.matrix} objects from the
#'package \code{snpMatrix}.
#'@section Warning: Any uncertain genotypes (stored by snpStats as raw codes 4
#'to 253) are output as missing.
#'
#'The functions use "\\n" as an end of line character, unless
#'\code{.Platform$OS.type == "windows"}, when eol is "\\r\\n".  I only have
#'access to linux machines for testing.
#'
#'I have tested these functions with my own data, but it is always possible
#'that your data may contain quirks mine don't, or that input formats could
#'change for any program mentioned here.  Please do have a quick check on a
#'small subset of data (eg, as in the example below), that the output for your
#'exact combination of options looks sensible and matches the specified input
#'format.
#'@author Chris Wallace
#'@references David Clayton (2012). snpStats: SnpMatrix and XSnpMatrix classes
#'and methods. R package version 1.6.0.  http://www-gene.cimr.cam.ac.uk/clayton
#'
#'phase/fastPhase: \url{http://stephenslab.uchicago.edu/software.html}
#'
#'beagle: \url{http://faculty.washington.edu/browning/beagle/beagle.html}
#'
#'IMPUTE: \url{http://mathgen.stats.ox.ac.uk/impute/impute_v2.html}
#'
#'MACH: \url{http://www.sph.umich.edu/csg/abecasis/MACH}
#'
#'snphap:
#'\url{https://www-gene.cimr.cam.ac.uk/staff/clayton/software/snphap.txt}
#'@keywords manip
#'@examples
#'
#'data(testdata,package="snpStats")
#'A.small <- Autosomes[1:6,1:10]
#'f <- tempfile()
#'## write in suitable format for snphap
#'nsnps <- ncol(A.small)
#'write.simple(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), gsep=" ",
#'              nullallele='0', file=f,
#'                 write.sampleid=FALSE)
#'unlink(f)
#'
write.simple <- function(X,a1,a2,
                         file,fsep="\t",gsep="",nullallele='N',write.header=TRUE,
                         transpose=FALSE,
                         write.sampleid=TRUE,bp=NULL,num.coding=FALSE) {
 ## standard function, can be tailored to most output
## a1 -> uncounted allele
## a2 -> counted allele
## write.simple(snp.matrix, a1, a2, output.file)
  check.args(X,a1,a2)
  if(!is.null(bp) & length(bp)!=length(a1))
    stop("require ncol(X) == length(bp)")
  if(nchar(fsep)>1 || nchar(gsep)>1)
    stop("fsep and gsep must be \"\" or single character strings\n")
  unlink(file)
  bp.do <- if(!is.null(bp)) 1 else 0 # 1 if include bp column after snp id
  if(transpose) { write.header <- FALSE }
  if(write.header) {
    if(write.sampleid)
      cat("sampleid",fsep,sep="",file=file)
    cat(colnames(X),file=file,sep=fsep,append=TRUE)
    cat("\n",file=file,append=TRUE)
  }
  res <- .C("write_simple", X@.Data,
            as.character(a1), as.character(a2),
            as.integer(bp),as.integer(bp.do),
            as.character(fsep),as.character(gsep),
            as.integer(num.coding),as.character(nullallele),as.integer(transpose),
            as.character(file), as.integer(nrow(X)),
            as.integer(ncol(X)), rownames(X),
            colnames(X), get.eol(),PACKAGE="snpStatsWriter")
  return(c(nrow(X), ncol(X)))
}

##' Simple wrapper to write.simple to write files in SNPHAP format
##'
##' If not allele codes are given, a1 and a2 will be set to 1 and 2 for all SNPs
##' @title Write SNPHAP files
##' @inheritParams write.simple
##' @return No return value, but has the side effect of writing specified output
##' file.
##' @author Chris Wallace
#'@export
#'@examples
#'
#'data(testdata,package="snpStats")
#'A.small <- Autosomes[1:6,1:10]
#'f <- tempfile()
#'## write in suitable format for snphap
#'write.snphap(A.small, file=f)
#'unlink(f)
#'
write.snphap <- function(X, a1=NULL, a2=NULL, file) {  
  if(is.null(a1) || is.null(a2)) {
    nsnps <- ncol(X)
    a1 <- rep("1",nsnps)
    a2 <- rep("2",nsnps)
  }
  valid.num <- as.character(1:2)
  valid.nuc <- c("A","C","G","T")
  coding <- "numeric"
  a1 <- as.character(a1)
  a2 <- as.character(a2)
  alleles <- c(unique(a1),unique(a2))
  if(!all(alleles %in% valid.num)) {
    coding <- "nucleotide"
    if(!all(alleles %in% valid.nuc)) {
      warning("detected nucleotide coding, but invalid alleles found.  Will recode to A/T.  To avoid this, please use write.simple(), but note that snphap recognises only 0/1/2 or A/C/G/T/0 coding.")
      which.bad <- which(!(a1 %in% valid.nuc | a2 %in% valid.nuc))
      anames <- colnames(X)
      for(i in which.bad) {
        cat("recoding SNP",i,":",anames[i],"from",a1[i],"/",a2[i]," -> A / T.\n")
        a1[i] <- "A"
        a2[i] <- "T"
      }
    }
  }  
  write.simple(X, a1=a1, a2=a2, gsep=" ",
               nullallele='0', file=file,
               write.sampleid=FALSE)
}

################################################################################

## mach requires 3 output files
##' Write a snpStats object in mach format
##'
##' see \code{\link{write.simple}} for general information
##' 
#' @export
##' @inheritParams write.simple
##' @param pedfile Output pedigree file name.
##' @param mfile Output marker file name.  
##' @param pedigree Optional pedigree/member/father/mother/sex indentifier vectors, same order as rows in snpStats object.  If missing, pedigree is set to rownames(X) and the others default to unrelated males
##' @param member See pedigree
##' @param father See pedigree
##' @param mother See pedigree
##' @param sex See pedigree
##' @param snp.names optional SNP names to include in the marker map file.  Defaults to colnames(X).
##' 
##' @return No return value, but has the side effect of writing specified output
##' files.
##' @author Chris Wallace
#'@keywords manip
#'@examples
#'
#'data(testdata,package="snpStats")
#'A.small <- Autosomes[1:6,1:10]
#'pf <- tempfile()
#'mf <- tempfile()
#'
#'## write in suitable format for MACH
#'nsnps <- ncol(A.small)
#'write.mach(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), pedfile=pf, mfile=mf)
#'unlink(pf)
#'unlink(mf)
#'
write.mach <- function(X,a1,a2,pedfile,mfile,
                       pedigree=rownames(X),
                       member=rep(1,nrow(X)),
                       father=rep(0,nrow(X)),
                       mother=rep(0,nrow(X)),
                       sex=rep("M",nrow(X)),
                       snp.names=colnames(X)) {
  check.args(X,a1,a2)
  if(length(snp.names) != ncol(X))
    stop("snp.names must have length == ncol(X)")
  res <- .C("write_mach", X@.Data, as.character(a1), as.character(a2),
            as.character(pedfile),
            as.character(pedigree), as.character(member), as.character(father),
            as.character(mother), as.character(sex),
            as.integer(nrow(X)),
            as.integer(ncol(X)), rownames(X),
            colnames(X), get.eol(),PACKAGE="snpStatsWriter")
  cat(paste("M",snp.names),file=mfile,sep="\n")
  return(c(nrow(X), ncol(X)))
}

################################################################################

## impute

##' Write a snpStats object in IMPUTE format
##'
##' see \code{\link{write.simple}} for general information
##' 
#' @export
##' @inheritParams write.simple
##' @param pedfile Output file name. 
##'@param snp.id vector of snp ids
##' 
##' @return No return value, but has the side effect of writing specified output
##' files.
##' @author Chris Wallace
#'@keywords manip
#'@examples
#'
#'data(testdata,package="snpStats")
#'A.small <- Autosomes[1:6,1:10]
#'pf <- tempfile()
#'
#'## write in suitable format for IMPUTE
#'nsnps <- ncol(A.small)
#'write.impute(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), bp=1:nsnps, pedfile=pf)
#'unlink(pf)
#'
write.impute <- function(X,a1,a2,bp,pedfile,snp.id=NULL) {
  check.args(X,a1,a2)
  if(length(bp)!=ncol(X))
    stop("require ncol(X) == length(bp)")
  if(is.null(snp.id))
    snp.id <- sprintf("SNP%s",1:ncol(X))
  res <- .C("write_impute", X@.Data, as.character(a1), as.character(a2),
            as.integer(bp), as.character(pedfile), as.integer(nrow(X)),
            as.integer(ncol(X)), rownames(X),
            colnames(X), as.character(snp.id), get.eol(),PACKAGE="snpStatsWriter")
  return(c(nrow(X), ncol(X)))
}

################################################################################

## beagle

##' Write a snpStats object in beagle format
##'
##' see \code{\link{write.simple}} for general information
##' 
#' @export
##' @inheritParams write.simple
#'@param trait disease trait (0=missing, 1=control, 2=case)
#'@param gfile,mfile \code{gfile}=genotype file, \code{pedfile}=pedigree file
##' 
##' @return No return value, but has the side effect of writing specified output
##' files.
##' @author Chris Wallace
#'@keywords manip
#'@examples
#'
#'data(testdata,package="snpStats")
#'A.small <- Autosomes[1:6,1:10]
#'gf <- tempfile()
#'mf <- tempfile()
#'
#'## write in suitable format for beagle
#'nsnps <- ncol(A.small)
#'write.beagle(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), bp=1:nsnps, gfile=gf, mfile=mf)
#'unlink(gf)
#'unlink(mf)
#'
write.beagle <- function(X,a1,a2,bp,trait=NULL,gfile,mfile) {
  check.args(X,a1,a2)
  if(length(bp)!=ncol(X))
    stop("require ncol(X) == length(bp)")
  if(!is.null(trait) & length(trait)!=nrow(X))
    stop("require trait (if given) to have length equal to ncol(X)")
  if(!is.null(trait)) {
    r <- range(trait)
    if(r[1]==0 & r[2]==2) {
      wh <- which(trait==0)
      cat("dropping",length(wh),"individuals with trait==0 (presumed missing)\n")
      trait <- trait[-wh]
      X <- X[-wh,]
    }
    if(r[1]==0 & r[2]==1) {
      trait <- trait+1 # require 1,2 labels
    }
  } else {
    trait <- numeric(0)
  }
  res <- .C("write_beagle", X@.Data, as.character(a1), as.character(a2),
            as.integer(bp), as.integer(trait),
            as.character(gfile), as.character(mfile),
            as.integer(nrow(X)),
            as.integer(ncol(X)), as.integer(length(trait)),
            rownames(X),
            colnames(X), get.eol(),PACKAGE="snpStatsWriter")
  return(c(nrow(X), ncol(X)))
}

################################################################################

## phase/fastphase
## beagle

##' Write a snpStats object in PHASE/FastPHASE format
##'
##' see \code{\link{write.simple}} for general information
##' 
#' @export
##' @inheritParams write.simple
##' @param file Output file name. 
##' 
##' @return No return value, but has the side effect of writing specified output
##' files.
##' @author Chris Wallace
#'@keywords manip
#'@examples
#'
#'data(testdata,package="snpStats")
#'A.small <- Autosomes[1:6,1:10]
#'f <- tempfile()
#'
#'## write in suitable format for PHASE
#'nsnps <- ncol(A.small)
#'write.phase(A.small, file=f)
#'unlink(f)
#'

write.phase <- function(X,a1=rep(1,ncol(X)),a2=rep(2,ncol(X)),bp=NULL,file) {
  check.args(X,a1,a2)
  if(!is.null(bp) && length(bp)!=ncol(X))
    stop("require ncol(X) == length(bp)")
  res <- .C("write_phase", X@.Data, as.character(a1), as.character(a2),
            as.character(file), as.integer(nrow(X)),
            as.integer(ncol(X)), rownames(X), colnames(X),
            as.integer(!is.null(bp)), as.integer(bp), get.eol(),PACKAGE="snpStatsWriter")
  return(c(nrow(X), ncol(X)))
}

################################################################################

##' write an sbams format file
##'
##' sbams is software from Xiaoquan Wen at https://github.com/xqwen/sbams
##' @title write.sbams
##' @inheritParams write.simple
##' @param response vector or matrix of response variables. rows index subjects, columns index variables
##' @return  No return value, but has the side effect of writing specified output
##' file.
##' @export
##' @author Chris Wallace
#'@keywords manip
#'@examples
#'
#'data(testdata,package="snpStats")
#'A.small <- Autosomes[1:6,1:10]
#'R <- matrix(rnorm(12),ncol=2)
#' colnames(R) <- c("var1","var2")
#'f <- tempfile()
#'
#'## write in suitable format for sbams
#'write.sbams(X=A.small, response=R, file=f)
#'unlink(f)
#'
write.sbams <- function(X,response,file) {
  if(is.data.frame(response)) ## data.frame response
    response <- as.matrix(response)
  if(!is.matrix(response)) ## vector response
    response <- as.matrix(response,ncol=1)
  if(nrow(response)!=nrow(X))
    stop("reponse matrix must have equal nrow() to X")

  ## genotypes need to be numeric
  N <- as(X,"numeric")

  ## all variables need to be zero centred, then transposed to be write a variable on each line
  N <- t(scale(N, center=TRUE, scale=FALSE))
  response <- t(scale(response, center=TRUE, scale=FALSE))

  ## add "response" or "covariate" to each rowname
  rownames(N) <- paste("covariate",rownames(N),sep=" ")
  rownames(response) <- paste("response",rownames(response),sep=" ")
    
  ## write the file
  write.table(response,file=file,row.names=TRUE,col.names=FALSE,quote=FALSE)
  write.table(N,file=file,row.names=TRUE,col.names=FALSE,quote=FALSE,append=TRUE)  
}
