#' Efficiently Read Sequencing Data (VCF format, METAL format) into R
#'
#' SeqMiner provides functions to easily load Variant Call Format (VCF) or METAL format into R
#'
#' The aim of this package is to save your time parsing large text file.
#' That means data processing time can be saved for other researches.
#' This packages requires Bgzip compressed and Tabix indexed files as input.
#' If input files contaians annotation by TabAnno (),
#' it is possible to extract information at the unit of genes.
#'
#' @docType package
#' @name SeqMiner
#' @useDynLib seqminer
#' @importFrom utils packageVersion
NULL

#' Check if the input is url e.g. http:// or ftp://
#' @param fileName character vector
#' @keywords internal
isURL <- function(fileName) {
  if (grepl(pattern = "^http://", fileName) ||
      grepl(pattern = "^ftp://", fileName) ) {
    return(TRUE)
  }
  return(FALSE)
}

local.file.exists <- function(fileName) {
  if (isURL(fileName)) {
    return(TRUE)
  }
  return(file.exists(fileName))
}

#' Check if the inputs are valid tabix range such as chr1:2-300
#' @param range characer vector
#' @export
#' @examples
#' valid <- isTabixRange(c("chr1:1-200", "X:1", "1:100-100", "chr1", "1:1-20,1:30-40"))
#' stopifnot(all(valid))
#' invalid <- isTabixRange(c(":1", "chr1::", ":-"))
#' stopifnot(all(!invalid))
isTabixRange <- function(range) {
  isValid <- function(x) {
    ret = strsplit(x = x, split = ":")[[1]]
    if (length(ret) == 1) { # e.g. "chr1"
      return(TRUE)
    }
    if (length(ret) != 2) {
      return(FALSE)
    }
    chrom = ret[1]
    if (nchar(chrom) == 0) {
      return (FALSE)
    }
    ret = strsplit(x = ret[2], split = "-")[[1]]
    if (length(ret) == 2) {
      beg = suppressWarnings(as.integer(ret[1]))
      end = suppressWarnings(as.integer(ret[2]))
      if (is.na(beg) || is.na(end) || beg > end) {
        return(FALSE)
      }
    } else if (length(ret) == 1) {
      beg = suppressWarnings(as.integer(ret[1]))
      if (is.na(beg)) {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
    return(TRUE)
  }
  ranges <- strsplit(x = range, split = ",")[[1]]
  sapply(ranges, isValid)
}

#' Check input file has tabix index
#'
#' @param fileName character vector, representing file names an input file name
#' @return TRUE if an index file with ".tbi" or "bci" exists, and created later than VCF/BCF file
#' @keywords internal
hasIndex <- function(fileName) {
  if (isURL(fileName)) {
    return(TRUE)
  }
  endsWith <- function(i, pattern = NULL) {
    if (is.null(pattern)) {
      return(FALSE)
    }
    if (length(fileName) != 1) {
      stop("hasIndex() only take vector of length 1.")
      return (FALSE)
    }
    p <- paste(pattern, "$", sep = "")
    if (length(grep(x=i, pattern=p)) != 1) {
      return (FALSE)
    }
    return (TRUE)
  }
  ## handle VCF
  fIndex <- NULL
  if (endsWith(fileName, ".vcf")) {
    stop(gettextf("Input file is not bgzip compressed. Please use: [ bgzip %s ]", fileName))
    return (FALSE)
  }
  if (endsWith(fileName, ".vcf.gz")) {
    fIndex <- paste(fileName, ".tbi", sep = "")
    if (!file.exists(fIndex)) {
      stop(gettextf("Cannot find index file (you can create it using: [ tabix -p vcf %s ])" , fileName))
      return (FALSE)
    }
  }

  ## handle BCF
  if (endsWith(fileName, ".bcf")) {
    stop(gettextf("Input file is not bgzip compressed. Please use: [ bgzip %s ]", fileName))
    return (FALSE)
  }
  if (endsWith(fileName, ".bcf.gz")) {
    fIndex <- paste(fileName, ".bci", sep = "")
    if (!file.exists(fIndex)) {
      stop(gettextf("Cannot find index file (you can create it using: [ bcftools index %s ])" , fileName))
      return (FALSE)
    }
  }

  ## Check index file
  if (is.null(fIndex)) {
    stop(gettextf("Cannot guess a valid tabix index, did your file have suffix: .vcf, .vcf.gz, .bcf, or .bcf.gz ?"))
    return (FALSE)
  }
  if (!file.exists(fIndex)) {
    stop(gettextf("The index file [ %s ] does not exists.", fIndex))
    return (FALSE)
  }

  ## Check index file timestamp
  time.data <- file.info(fileName)$mtime
  time.index <- file.info(fIndex)$mtime
  if (time.data > time.index) {
    stop(gettextf("Index file [ %s ] is older than data file [ %s ].", fIndex, fileName))
    return(FALSE)
  }
  return (TRUE)
}

#' Read a gene from VCF file and return a genotype matrix
#'
#' @param fileName character, represents an input VCF file (Bgzipped, with Tabix index)
#' @param range character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @param annoType character, annotated types you would like to extract, such as "Nonsynonymous", "Synonymous". This can be left empty.
#' @return genotype matrix
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' cfh <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Nonsynonymous")
readVCFToMatrixByRange <- function(fileName, range, annoType) {
  stopifnot(local.file.exists(fileName), length(fileName) == 1)
  stopifnot(hasIndex(fileName))
  stopifnot(all(isTabixRange(range)))
  fileName <- path.expand(fileName)

  storage.mode(fileName) <- "character"
  storage.mode(range)    <- "character"
  storage.mode(annoType) <- "character"
  .Call("readVCFToMatrixByRange", fileName, range, annoType, PACKAGE="seqminer");
}

#' Read a gene from VCF file and return a genotype matrix
#'
#' @param fileName character, represents an input VCF file (Bgzipped, with Tabix index)
#' @param geneFile character, a text file listing all genes in refFlat format
#' @param geneName character vector, which gene(s) to be extracted
#' @param annoType character, annotated types you would like to extract, such as "Nonsynonymous", "Synonymous". This can be left empty.
#' @return genotype matrix
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- readVCFToMatrixByGene(fileName, geneFile, "CFH", "Synonymous")
readVCFToMatrixByGene <- function(fileName, geneFile, geneName, annoType) {
  stopifnot(local.file.exists(fileName), length(fileName) == 1)
  stopifnot(local.file.exists(geneFile), length(geneFile) == 1)
  stopifnot(hasIndex(fileName))
  fileName <- path.expand(fileName)

  storage.mode(fileName) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(geneName) <- "character"
  storage.mode(annoType) <- "character"
  .Call("readVCFToMatrixByGene", fileName, geneFile, geneName, annoType, PACKAGE="seqminer");
}

#' Read information from VCF file in a given range and return a list
#'
#' @param fileName character, represents an input VCF file (Bgzipped, with Tabix index)
#' @param range character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @param annoType character, annotated types you would like to extract, such as "Nonsynonymous", "Synonymous". This can be left empty.
#' @param vcfColumn character vector, which vcf columns to extract. It can be chosen from CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT and etc.
#' @param vcfInfo character vector, which should be tags in the INFO columns to extarct. Common choices include: DP, AC, AF, NS
#' @param vcfIndv character vector, which values to extract at individual level. Common choices are: GT, GQ, GD
#' @return a list of genes, and each elements has specified vcfColumn, vcfinfo, vcfIndv
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' cfh <- readVCFToListByRange(fileName, "1:196621007-196716634", "Nonsynonymous",
#'                             c("CHROM", "POS"), c("AF", "AC"), c("GT") )
readVCFToListByRange <- function(fileName, range, annoType, vcfColumn, vcfInfo, vcfIndv) {
  stopifnot(local.file.exists(fileName), length(fileName) == 1)
  stopifnot(hasIndex(fileName))
  stopifnot(all(isTabixRange(range)))
  fileName <- path.expand(fileName)

  storage.mode(fileName) <- "character"
  storage.mode(range)    <- "character"
  storage.mode(annoType) <- "character"
  storage.mode(vcfColumn)<- "character"
  storage.mode(vcfInfo)  <- "character"
  storage.mode(vcfIndv)  <- "character"
  .Call("readVCFToListByRange", fileName, range, annoType, vcfColumn, vcfInfo, vcfIndv, PACKAGE="seqminer");
}

#' Read information from VCF file in a given range and return a list
#'
#' @param fileName character, represents an input VCF file (Bgzipped, with Tabix index)
#' @param geneFile character, a text file listing all genes in refFlat format
#' @param geneName character vector, which gene(s) to be extracted
#' @param annoType character, annotated types you would like to extract, such as "Nonsynonymous", "Synonymous". This can be left empty.
#' @param vcfColumn character vector, which vcf columns to extract. It can be chosen from CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT and etc.
#' @param vcfInfo character vector, which should be tags in the INFO columns to extarct. Common choices include: DP, AC, AF, NS
#' @param vcfIndv character vector, which values to extract at individual level. Common choices are: GT, GQ, GD
#' @return a list of genes, and each elements has specified vcfColumn, vcfinfo, vcfIndv
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous",
#'                            c("CHROM", "POS"), c("AF", "AC"), c("GT") )
readVCFToListByGene <- function(fileName, geneFile, geneName, annoType, vcfColumn, vcfInfo, vcfIndv) {
  stopifnot(local.file.exists(fileName), length(fileName) == 1)
  stopifnot(file.exists(geneFile), length(geneFile) == 1)
  stopifnot(hasIndex(fileName))
  fileName <- path.expand(fileName)

  storage.mode(fileName) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(geneName) <- "character"
  storage.mode(annoType) <- "character"
  storage.mode(vcfColumn)<- "character"
  storage.mode(vcfInfo)  <- "character"
  storage.mode(vcfIndv)  <- "character"
  .Call("readVCFToListByGene", fileName, geneFile, geneName, annoType, vcfColumn, vcfInfo, vcfIndv, PACKAGE="seqminer");
}

#' Read association statistics by gene from METAL-format files. Both score statistics and covariance statistics will be extracted.
#'
#' @param scoreTestFiles character vector, score test output files (rvtests outputs using --meta score)
#' @param covFiles character vector, covaraite files (rvtests outputs using --meta cov)
#' @param geneFile character, a text file listing all genes in refFlat format
#' @param geneName character vector, which gene(s) to be extracted
#' @return a list of statistics including chromosome, position, allele frequency, score statistics, covariance and annotation(if input files are annotated).
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- rvmeta.readDataByGene(scoreFileName, covFileName, geneFile, "CFH")
rvmeta.readDataByGene <- function(scoreTestFiles, covFiles, geneFile, geneName) {
  stopifnot(local.file.exists(scoreTestFiles))
  stopifnot(is.null(covFiles) || (local.file.exists(covFiles) && length(covFiles) == length(scoreTestFiles)))
  stopifnot(file.exists(geneFile), length(geneFile) == 1)
  scoreTestFiles <- path.expand(scoreTestFiles)
  if (!is.null(covFiles)) {
    covFiles <- path.expand(covFiles)
  }

  storage.mode(scoreTestFiles) <- "character"
  storage.mode(covFiles) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(geneName) <- "character"
  if (is.null(covFiles)) {
    .Call("rvMetaReadDataByGene", scoreTestFiles, "", geneFile, geneName, PACKAGE="seqminer");
  } else {
    .Call("rvMetaReadDataByGene", scoreTestFiles, covFiles, geneFile, geneName, PACKAGE="seqminer");
  }
}

#' Read association statistics by range from METAL-format files. Both score statistics and covariance statistics will be extracted.
#'
#' @param scoreTestFiles character vector, score test output files (rvtests outputs using --meta score)
#' @param covFiles character vector, covaraite files (rvtests outputs using --meta cov)
#' @param ranges character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return a list of statistics including chromosome, position, allele frequency, score statistics, covariance and annotation(if input files are annotated).
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- rvmeta.readDataByRange(scoreFileName, covFileName, "1:196621007-196716634")
rvmeta.readDataByRange <- function (scoreTestFiles, covFiles, ranges) {
  stopifnot(local.file.exists(scoreTestFiles))
  stopifnot(is.null(covFiles) || (local.file.exists(covFiles) && length(covFiles) == length(scoreTestFiles)))
  stopifnot(all(isTabixRange(ranges)))
  scoreTestFiles <- path.expand(scoreTestFiles)
  if (!is.null(covFiles)) {
    covFiles <- path.expand(covFiles)
  }

  storage.mode(scoreTestFiles) <- "character"
  storage.mode(covFiles) <- "character"
  storage.mode(ranges) <- "character"
  if (is.null(covFiles)) {
    .Call("rvMetaReadDataByRange", scoreTestFiles, "",
          ranges, PACKAGE = "seqminer")
  } else {
    .Call("rvMetaReadDataByRange", scoreTestFiles, covFiles,
          ranges, PACKAGE = "seqminer")
  }
}

#' Read covariance by range from METAL-format files.
#'
#' @param covFile character, a covariance file (rvtests outputs using --meta cov)
#' @param tabixRange  character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return a matrix of covariance within given range
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' cfh <- rvmeta.readCovByRange(covFileName, "1:196621007-196716634")
rvmeta.readCovByRange <- function(covFile, tabixRange) {
  stopifnot(local.file.exists(covFile))
  stopifnot(all(isTabixRange(tabixRange)))
  covFile <- path.expand(covFile)

  storage.mode(covFile) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readCovByRange", covFile, tabixRange, PACKAGE="seqminer");
}

#' Read score test statistics by range from METAL-format files.
#'
#' @param scoreTestFiles character vector, score test output files (rvtests outputs using --meta score)
#' @param tabixRange  character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return score test statistics within given range
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' cfh <- rvmeta.readScoreByRange(scoreFileName, "1:196621007-196716634")
rvmeta.readScoreByRange <- function(scoreTestFiles, tabixRange) {
  stopifnot(local.file.exists(scoreTestFiles))
  stopifnot(all(isTabixRange(tabixRange)))
  scoreTestFiles <- path.expand(scoreTestFiles)

  storage.mode(scoreTestFiles) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readScoreByRange", scoreTestFiles, tabixRange, PACKAGE="seqminer");
}

#' Read skew by range from METAL-format files.
#'
#' @param skewFile character, a skew file (rvtests outputs using --meta skew)
#' @param tabixRange  character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return an 3-dimensional array of skewness within given range
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' skewFileName = system.file("rvtests/rvtest.MetaSkew.assoc.gz", package = "seqminer")
#' cfh <- rvmeta.readSkewByRange(skewFileName, "1:196621007-196716634")
rvmeta.readSkewByRange <- function(skewFile, tabixRange) {
  stopifnot(local.file.exists(skewFile))
  stopifnot(all(isTabixRange(tabixRange)))
  skewFile <- path.expand(skewFile)

  storage.mode(skewFile) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readSkewByRange", skewFile, tabixRange, PACKAGE="seqminer");
}

#' Read tabix file, similar to running tabix in command line.
#'
#' @param tabixFile character, an tabix indexed file
#' @param tabixRange  character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return character vector, each elements is an individual line
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' snp <- tabix.read(fileName, "1:196623337-196632470")
tabix.read <- function(tabixFile, tabixRange) {
  stopifnot(local.file.exists(tabixFile))
  stopifnot(all(isTabixRange(tabixRange)))
  tabixFile <- path.expand(tabixFile)

  storage.mode(tabixFile) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readTabixByRange", tabixFile, tabixRange, PACKAGE="seqminer");
}

#' Read tabix file, similar to running tabix in command line.
#'
#' @param tabixFile character, an tabix indexed file
#' @param skippedLine  logical, whether to read tabix skipped lines (when used 'tabix -S NUM')
#' @return a list
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' snp <- tabix.read.header(fileName)
tabix.read.header <- function(tabixFile, skippedLine = FALSE) {
  stopifnot(local.file.exists(tabixFile))
  tabixFile <- path.expand(tabixFile)

  storage.mode(tabixFile) <- "character"
  header <- .Call("readTabixHeader", tabixFile, PACKAGE="seqminer");
  ret <- list(header = header)
  if (skippedLine) {
    skipped <- .Call("readTabixSkippedLine", tabixFile, PACKAGE="seqminer");
    ret[[length(ret) + 1]] <- skipped
    names(ret)[2] <- "skippedLine"
  }
  ret
}

#' Read tabix file, similar to running tabix in command line.
#'
#' @param tabixFile character, an tabix indexed file
#' @param tabixRange  character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @param col.names logical, use tabix file header as result headers (default: TRUE)
#' @param stringsAsFactors logical, store loaded data as factors (default: FALSE)
#' @return data frame, each elements is an individual line
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' snp <- tabix.read.table(fileName, "1:196623337-196632470")
tabix.read.table <- function(tabixFile, tabixRange, col.names = TRUE, stringsAsFactors = FALSE) {
  stopifnot(local.file.exists(tabixFile))
  stopifnot(all(isTabixRange(tabixRange)))
  tabixFile <- path.expand(tabixFile)

  storage.mode(tabixFile) <- "character"
  storage.mode(tabixRange) <- "character"
  header <- .Call("readTabixHeader", tabixFile, PACKAGE = "seqminer")
  body <- .Call("readTabixByRange", tabixFile, tabixRange, PACKAGE="seqminer");

  ## parse body to a table
  body <- do.call(rbind, strsplit(body, "\t"))
  body <- as.data.frame(body, stringsAsFactors = FALSE)
  if (ncol(body) > 0) {
    for (i in 1:ncol(body)) {
      body[,i] <- utils::type.convert(body[,i], as.is = !stringsAsFactors)
    }

    num.col <- ncol(body)
    header <- header[nchar(header) > 0]
    if (length(header) == 0 || !col.names) {
      colNames <- paste0("V", 1L:num.col)
    } else {
      hdrLine <- header[length(header)]
      hdrLine <- sub("^#", "", hdrLine)
      colNames <- make.names(strsplit(hdrLine, "\t")[[1]])
      if (length(colNames) > ncol(body)) {
        colNames <- colNames[1:ncol(body)]
      } else if (length(colNames) < ncol(body)) {
        tmpNames <- paste0("V", 1L:num.col)
        tmpNames[1:length(colNames )] <- colNames
        colNames <- tmpNames
      }
    }
    colnames(body) <- colNames
  }

  body
}

.onAttach <- function(libname, pkgname){
  # check new version
  newVersionLink = "http://zhanxw.com/seqminer/version"
  timeout <- options("timeout")$timeout
  warn <- options("warn")$warn
  options(timeout = 3)
  options(warn = -1)
  conn <- url(newVersionLink)
  ret <- tryCatch(readLines(conn, n = 2), error = function(e) {NULL})
  close(conn)
  options(warn = warn)
  options(timeout = timeout)

  if (!is.null(ret) && length(ret) == 2) {
    version <- ret[1]
    if (utils::packageVersion("seqminer") < version) {
      if (length(ret) > 1) {
        packageStartupMessage(ret[2])
      } else {
        packageStartupMessage("Found new version of seqminer: ", ret)
      }
    }
  }
  ## packageStartupMessage("seqminer loaded ...")
}

#' Read null model statistics
#'
#' @param scoreTestFiles character vector, score test output files (rvtests outputs using --meta score)
#' @return a list of statistics fitted under the null mode (without genetic effects)
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
rvmeta.readNullModel <- function(scoreTestFiles) {
  stopifnot(local.file.exists(scoreTestFiles))
  read.null.model <- function(fn) {
    ## cat("gzFile, fn = ", fn, "\n")
    f <- gzfile(fn, "rb")
    ret = list()
    beginStore = FALSE
    while(TRUE) {
      suppressWarnings(p <- readLines(con = f, n = 1, warn = FALSE))
      if (substr(p, 1, 1) != "#") {
        break
      }
      #print(p)
      if (!beginStore) {
        if ( p[1] == "##NullModelEstimates") {
          beginStore = TRUE
          next
        }
      } else {
        ret[length(ret) + 1] <- p
      }
    }
    close(f)
    ret
  }
  model.to.matrix <- function(ret) {
    if (is.null(ret) || length(ret) < 1) {
      return(numeric(0))
    }
    ## ret <- read.null.model(ret)
    ret <- lapply(ret, function(x) { strsplit(sub("## - ", "", x), "\t")[[1]]})
    cnames <- ret[[1]]
    ret[[1]] <- NULL
    rnames <- lapply(ret, function(x) {x[1]})
    rnames <- do.call(c, rnames)
    data <- lapply(ret, function(x) {x[-1]})
    mat <- do.call(rbind, data)
    colnames(mat) <- cnames[-1]
    rownames(mat) <- rnames
    suppressWarnings(class(mat) <- "numeric")
    mat
  }
  ret <- lapply(scoreTestFiles, function(x) { model.to.matrix(read.null.model(x))})
  ret
}

#' Write score-based association statistics files.
#'
#' @param rvmetaData a list vector. It's usually read by rvmeta.readDataByRange or rvmeta.readDataByGene function
#' @param outName character, a text indicating output file prefix
#' @param createIndex boolean, (default FALSE), whether or not to create the index
#' @return TRUE only if succeed
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- rvmeta.readDataByRange(scoreFileName, covFileName, "1:196621007-196716634")
#' rvmeta.writeScoreData(cfh, "cfh.MetaScore.assoc")
rvmeta.writeScoreData <- function (rvmetaData, outName, createIndex = FALSE) {
  outName <- path.expand(outName)
  storage.mode(outName) <- "character"
  .Call("rvMetaWriteScoreData", rvmetaData, outName, PACKAGE="seqminer");
  if (createIndex) {
    tabix.createIndex.meta(outName)
  }
  invisible(NULL)
}

#' Write covariance association statistics files.
#'
#' @param rvmetaData a list vector. It's usually read by rvmeta.readDataByRange or rvmeta.readDataByGene function
#' @param outName character, a text indicating output file prefix
#' @return TRUE only if succeed
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- rvmeta.readDataByRange(scoreFileName, covFileName, "1:196621007-196716634")
#' rvmeta.writeCovData(cfh, "cfh.MetaCov.assoc.gz")
rvmeta.writeCovData <- function (rvmetaData, outName) {
  outName <- path.expand(outName)
  storage.mode(outName) <- "character"
  .Call("rvMetaWriteCovData", rvmetaData, outName, PACKAGE="seqminer");
  tabix.createIndex.meta(outName)
  invisible(NULL)
}

#' Create tabix index file, similar to running tabix in command line.
#'
#' @param bgzipFile character, an tabix indexed file
#' @param sequenceColumn  integer, sequence name column
#' @param startColumn  integer, start column
#' @param endColumn  integer, end column
#' @param metaChar  character, symbol for comment/meta lines
#' @param skipLines  integer, first this number of lines will be skipped
#' @return NULL
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' tabix.createIndex(fileName, 1, 2, 0, '#', 0)
tabix.createIndex <- function(bgzipFile, sequenceColumn = 1, startColumn = 4, endColumn = 5, metaChar = "#", skipLines = 0) {
  stopifnot(file.exists(bgzipFile))
  bgzipFile <- path.expand(bgzipFile)
  storage.mode(bgzipFile) <- "character"
  storage.mode(sequenceColumn) <- "integer"
  storage.mode(startColumn) <- "integer"
  storage.mode(endColumn) <- "integer"
  storage.mode(metaChar) <- "character"
  storage.mode(skipLines) <- "integer"
  .Call("createTabixIndex", bgzipFile, sequenceColumn, startColumn, endColumn, metaChar, skipLines, PACKAGE="seqminer");
  invisible(NULL)
}

#' Create tabix index for bgzipped VCF file
#'
#' @param bgzipVcfFile character, input vcf file
#' @return NULL
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' tabix.createIndex.vcf(fileName)
tabix.createIndex.vcf <- function(bgzipVcfFile) {
  stopifnot(file.exists(bgzipVcfFile))
  tabix.createIndex(bgzipVcfFile, 1, 2, 0, '#', 0)
  invisible(NULL)
}

#' Create tabix index for bgzipped MetaScore/MetaCov file
#  MetaScore/MetaCov files can be generated by rvtests.
#'
#' @param bgzipFile character, input vcf file
#' @return NULL
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @seealso http://zhanxw.github.io/rvtests/ for rvtests
#' @examples
#' fileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' tabix.createIndex.meta(fileName)
tabix.createIndex.meta <- function(bgzipFile) {
  stopifnot(file.exists(bgzipFile))
  tabix.createIndex(bgzipFile, 1, 2, 2, '#', 0)
  invisible(NULL)
}

##################################################
## Annotations --------------------
##################################################

#' validate the inVcf can be created, and outVcf can be write to.
#' will stop if any error occurs
#' @param inVcf input file
#' @param outVcf output file
#' @return NULL
verifyFilename <- function(inVcf, outVcf) {
  if (!file.exists(inVcf)) {
    stop("Cannot open input file")
  }
  tryCatch(file.create(outVcf), warning = function(w) {stop(gettextf("Cannot create output file: %s", outVcf))})
}

#' Validate annotate parameter is valid
#' @param param a list of annotation elements
#' @param debug show extra debug information or not
#' @return list, first element is TRUE/FALSE if parameter is valid/invalid;
#                second element is a message string
#' @export
validateAnnotationParameter <- function(param, debug = FALSE) {
  status <- TRUE
  msg <- vector("character")
  if (debug) {
    print(param)
  }
  if (is.null(param$reference) || !file.exists(param$reference)) {
    msg <- c(msg, "Reference file does not exist")
  }
  if (is.null(param$geneFile) || !file.exists(param$geneFile)) {
    status <- FALSE
    msg <- c(msg, "Gene file does not exist")
  }
  if (! param$geneFileFormat %in% c("refFlat", "knownGene", "refGene")) {
    status <- FALSE
    msg <- c(msg, "Gene file format does not exist")
  }
  if (is.null(param$codonFile) || !file.exists(param$codonFile)) {
    status <- FALSE
    msg <- c(msg, "Codon file does not exist")
  }
  if (is.null(param$priorityFile) || !file.exists(param$priorityFile)) {
    status <- FALSE
    msg <- c(msg, "Priority file does not exist")
  }
  if (param$upstreamRange <= 0) {
    status <- FALSE
    msg <- c(msg, "Upstream range is less than zero")
  }
  if (param$downstreamRange <= 0) {
    status <- FALSE
    msg <- c(msg, "Downstream range is less than zero")
  }
  if (param$spliceIntoExon <= 0) {
    status <- FALSE
    msg <- c(msg, "SpliceIntoExon range is less than zero")
  }
  if (param$spliceIntoIntron <= 0) {
    status <- FALSE
    msg <- c(msg, "SpliceIntoIntron range is less than zero")
  }
  if (!is.null(param$bed)) {
    ##library(stringr)
    opt <- strsplit(param$bed, ",")
    n <- length(opt)
    parsed <- lapply(opt, function(x){unlist(strsplit(x, "="))})
    for (i in seq_along(n)) {
      if (length(parsed[[i]]) != 2) {
        status <- FALSE
        msg <- c(msg, paste("Cannot understand bed option: ", parsed[[i]]))
        next
      }
      if (!file.exists(parsed[[i]][2])) {
        status <- FALSE
        msg <- c(msg, paste("BED resource file not exists: ", parsed[[i]][2]))
        next
      }
    }
  }
  if (!is.null(param$genomeScore)) {
    ##library(stringr)
    opt <- strsplit(param$genomeScore, ",")
    n <- length(opt)
    parsed <- lapply(opt, function(x){unlist(strsplit(x, "="))})
    for (i in seq_along(n)) {
      if (length(parsed[[i]]) != 2) {
        status <- FALSE
        msg <- c(msg, paste("Cannot understand genomeScore option: ", parsed[[i]]))
        next
      }
      if (!file.exists(parsed[[i]][2])) {
        status <- FALSE
        msg <- c(msg, paste("Genomescore resource file not exists: ", parsed[[i]][2]))
        next
      }
    }
  }
  if (!is.null(param$tabix)) {
    ##library(stringr)
    opt <- strsplit(param$tabix, "),")
    n <- length(opt)
    parsed <- lapply(opt, function(x){unlist(strsplit(x, "\\("))})
    for (i in seq_along(n)) {
      if (length(parsed[[i]]) != 2) {
        status <- FALSE
        msg <- c(msg, paste("Cannot understand tabix option: ", parsed[[i]]))
        next
      }
      if (!file.exists(parsed[[i]][1])) {
        status <- FALSE
        msg <- c(msg, paste("Tabix resource file not exists: ", parsed[[i]][1]))
        next
      }
    }
  }

  res <- list(status, paste(msg, sep = " \n"))
  return(res)
}

#' Construct a usable set of annotation parameters
#' @param param a list of annotation elements
#' @return list, a complete list of supported parameters
#' @export
makeAnnotationParameter <- function(param = NULL) {
  ## default parameters
  defaultParam <- list(reference = NULL,
                       geneFile = NULL,
                       geneFileFormat = "refFlat",
                       codonFile = system.file("tabanno/codon.txt", package = "seqminer"),
                       priorityFile = system.file("tabanno/priority.txt", package = "seqminer"),
                       upstreamRange = 50,
                       downstreamRange = 50,
                       spliceIntoExon = 3,
                       spliceIntoIntron = 8,
                       bed = NULL,
                       genomeScore = NULL,
                       tabix = NULL,
                       indexOutput = FALSE,
                       inputFormat = "vcf",
                       checkReference = TRUE)
  ret <- defaultParam
  for (i in seq_along(param)) {
    if (names(param)[i] %in% names(ret)) {
      ret[[names(param)[i]]] <- param[[i]]
    }
  }
  ## expand tilde
  if(!is.null(ret$reference)) {
    ret$reference <- path.expand(ret$reference)
  }
  if(!is.null(ret$geneFile)) {
    ret$geneFile <- path.expand(ret$geneFile)
  }
  if(!is.null(ret$codonFile)) {
    ret$codonFile <- path.expand(ret$codonFile)
  }
  if(!is.null(ret$priorityFile)) {
    ret$priorityFile <- path.expand(ret$priorityFile)
  }
  if(!is.null(ret$bed)) {
    ret$bed <- gsub(pattern = "~", replacement = path.expand("~"), x = ret$bed)
  }
  if(!is.null(ret$tabix)) {
    ret$tabix <- gsub(pattern = "~", replacement = path.expand("~"), x = ret$tabix)
  }

  ## print(ret)
  ret
}

#' Annotate a test variant
#' @param param a list of annotation configuaration (e.g. reference file, gene definition)
#' @param chrom a vector of chromosome names
#' @param position a vector of chromosome positions
#' @param ref a vector of reference alleles
#' @param alt a vector of alternative alleles
#' @return annotated results in a data frame structure
#' @export
#' @seealso makeAnnotationParameter
#' @examples
#' param <- list(reference = system.file("tabanno/test.fa", package = "seqminer"),
#'               geneFile = system.file("tabanno/test.gene.txt", package = "seqminer"))
#' param <- makeAnnotationParameter(param)
#' print(param)
#' annotateGene(param, c("1", "1"), c(3, 5) , c("A", "C"), c("G", "C"))
annotateGene <- function(param, chrom, position, ref, alt) {
  param <- makeAnnotationParameter(param)
  res <- validateAnnotationParameter(param)
  if (!res[[1]])  {
    cat(paste(res[[2]], collapse = "\n"))
    stop("Stop due to critical error")
  }
  stopifnot(length(chrom) > 0)
  stopifnot(length(chrom) == length(position))
  stopifnot(length(chrom) == length(ref))
  stopifnot(length(chrom) == length(alt))

  storage.mode(chrom) <- "character"
  storage.mode(position) <- "integer"
  storage.mode(ref) <- "character"
  storage.mode(alt) <- "character"

  .Call("annotateGene", param, chrom, position, ref, alt, PACKAGE="seqminer")
}

#' Annotate a test variant
#' @param reference path to the reference genome file (.fa file)
#' @param chrom a vector of chromosome names
#' @param position a vector of chromosome positions
#' @param len a vector of length
#' @return based extracted from the reference genome
#' @export
getRefBase <- function(reference, chrom, position, len = NULL) {
  if (!file.exists(reference)) {
    stop("Reference file does not exists")
  }
  if (is.null(len)) {
    len <- rep(1, length(position))
  }
  reference <- path.expand(reference)

  storage.mode(reference) <- "character"
  storage.mode(chrom) <- "character"
  storage.mode(position) <- "integer"
  storage.mode(len) <- "integer"
  .Call("getRefBase", reference, chrom, position, len, PACKAGE="seqminer")
}

#' Annotate a VCF file
#' @param inVcf input VCF file name
#' @param outVcf output VCF file name
#' @param params parameters
#' @return 0 if succeed
#' @export
#' @examples
#' param <- list(reference = system.file("tabanno/test.fa", package = "seqminer"),
#'               geneFile = system.file("tabanno/test.gene.txt", package = "seqminer"))
#' param <- makeAnnotationParameter(param)
#' inVcf <- system.file("tabanno/input.test.vcf", package = "seqminer")
#' outVcf <- paste0(getwd(), "/", "out.vcf")
#' annotateVcf (inVcf, outVcf, param)
annotateVcf <- function(inVcf, outVcf, params) {
  params$inputFormat = "vcf"
  param <- makeAnnotationParameter(param)
  res <- validateAnnotationParameter(param)
  if (!res[[1]])  {
    cat(paste(res[[2]], collapse = "\n"))
    stop("Stop due to critical error")
  }

  verifyFilename(inVcf, outVcf)

  storage.mode(inVcf) <- "character"
  storage.mode(outVcf) <- "character"
  .Call("anno", inVcf, outVcf, params)
  invisible(NULL)
}

#' Annotate a plain text file
#' @param inFile input file name
#' @param outFile output file name
#' @param params parameters
#' @return 0 if succeed
#' @export
#' @examples
#' param <- list(reference = system.file("tabanno/test.fa", package = "seqminer"),
#'               geneFile = system.file("tabanno/test.gene.txt", package = "seqminer"),
#'               inputFormat = "plain")
#' param <- makeAnnotationParameter(param)
#' inFile <- system.file("tabanno/input.test.plain.txt", package = "seqminer")
#' outFile <- paste0(getwd(), "/", "out.annotated.txt")
#' annotatePlain(inFile, outFile, param)
annotatePlain <- function(inFile, outFile, params) {
  params$inputFormat = "plain"
  param <- makeAnnotationParameter(param)
  res <- validateAnnotationParameter(param)
  if (!res[[1]])  {
    cat(res[[2]])
    stop("Stop due to critical error")
  }
  verifyFilename(inFile, outFile)
  inFile <- path.expand(inFile)
  outFile <- path.expand(outFile)

  storage.mode(inFile) <- "character"
  storage.mode(outFile) <- "character"
  .Call("anno", inFile, outFile, params)
  invisible(NULL)
}

#' Test whether a vector of positions are inside given ranges
#' @param positions characters, positions. e.g. c("1:2-3", "1:4")
#' @param rangeList character, ranges, e.g. "1:1-3,1:2-4"
#' @return logical vector, TRUE/FALSE/NA
#' @export
#' @examples
#' positions <- c("1:2-3", "1:4", "XX")
#' ranges <- "1:1-3,1:2-4,1:5-10"
#' isInRange(positions, ranges)
isInRange <- function(positions, rangeList) {
  storage.mode(positions) <- "character"
  storage.mode(rangeList) <- "character"
  .Call("isInRange", positions, rangeList)
}

#' Extract pair of positions by ranges
#' @param covData a covariance matrix with positions as dimnames
#' @param rangeList1 character specify a range
#' @param rangeList2 character specify a range
#' @return a covariance matrix
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' cfh <- rvmeta.readCovByRange(covFileName, "1:196621007-196716634")
#' rangeList1 <- "1:196621007-196700000"
#' rangeList2 <- "1:196700000-196716634"
#' getCovPair(cfh, rangeList1, rangeList2)
getCovPair <- function(covData, rangeList1, rangeList2) {
  pos <- dimnames(covData)[1][[1]]
  idx1 <- isInRange(pos, rangeList1)
  idx2 <- isInRange(pos, rangeList2)
  ret <- covData[idx1, idx2]
  dimnames(ret) <- list(pos[idx1], pos[idx2])
  ret
}

#' Test whether directory is writable
#' @param outDir the name of the directory
#' @return TRUE if the file is writable
#' isDirWritable("~")
isDirWritable <- function(outDir) {
  name <- tempfile(tmpdir = outDir, fileext = "tmp")
  ret <- file.create(name, showWarnings = FALSE)
  if (ret) {
    unlink(name)
  }
  return (ret)
}

#' Download annotation resources to a directory
#' @param outputDirectory the directory to store annotation resources
#' @return will not return anything
#' @export
#' @examples
#' \dontrun{
#' download.annotation.resource("/tmp")
#' }
download.annotation.resource <- function(outputDirectory) {
  outDir = outputDirectory
  ## prepare a writable dir
  if (!dir.exists(outDir)) {
    message(gettextf("Create output directory: %s", outDir))
    dir.create(outDir, recursive = TRUE)
  }
  if (!isDirWritable(outDir)) {
    stop(gettextf("Unable to write to directory: %s", outDir))
  }

  ## download function
  download <- function(url) {
    fn <- basename(url)
    destfile <- file.path(outDir, fn)
    if (file.exists(destfile)) {
      warning(gettextf("Overwriting %s", fn))
    }
    utils::download.file(url, destfile)
  }
  ## download resources
  message("Begin download TabAnno resource files (human hg19)...")


  message("Download reference file and its index:")
  download("http://qbrc.swmed.edu/zhanxw/software/anno/resources/hs37d5.fa")
  download("http://qbrc.swmed.edu/zhanxw/software/anno/resources/hs37d5.fa.fai")

  message("Download gene definition:")
  download("http://qbrc.swmed.edu/zhanxw/software/anno/resources/refFlat_hg19.txt.gz")

  message("Download TabAnno codon definition and annotation priority files:")
  download("http://qbrc.swmed.edu/zhanxw/software/anno/codon.txt")
  download("http://qbrc.swmed.edu/zhanxw/software/anno/priority.txt")

  message("Download completed")
  message(gettextf("You can begin to use it:"))
  message(gettextf(" param <- makeAnnotationParameter(list(reference = \"%s\", geneFile = \"%s\", codonFile = \"%s\", priorityFile = \"%s\" ))",
                   file.path(outDir, "hs37d5.fa"),
                   file.path(outDir, "refFlat_hg19.txt.gz"),
                   file.path(outDir, "codon.txt"),
                   file.path(outDir, "priority.txt")))
  invisible(NULL)
}

