#' @name fastsimcoal
#' @title Run fastsimcoal
#' @description Run a fastsimcoal simulation and load results into a 
#'   \linkS4class{gtypes} object.
#'
#' @param pop.info matrix of population sampling information created by the 
#'   \code{\link{fscPopInfo}} function.
#' @param locus.params data.frame specifying loci to simulate created by the 
#'   \code{\link{fscLocusParams}} function.
#' @param mig.rates a matrix or list of matrices giving the migration rates 
#'   between pairs of populations.
#' @param hist.ev matrix of historical events created by the 
#'   \code{\link{fscHistEv}} function.
#' @param label character string to label files with.
#' @param quiet logical. Run quietly?
#' @param delete.files logical. Delete files when done?
#' @param exec name of fastsimcoal executable.
#' @param num.cores number of cores to use.
#' 
#' @note fastsimcoal is not included with \code{strataG} and must be downloaded 
#'   separately. Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs} for installation instructions.
#' 
#' @references Excoffier, L. and Foll, M (2011) fastsimcoal: a continuous-time 
#'   coalescent simulator of genomic diversity under arbitrarily complex 
#'   evolutionary scenarios Bioinformatics 27: 1332-1334.\cr
#'   \url{http://cmpg.unibe.ch/software/fastsimcoal2/}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
NULL


.fscWrite <- function(pop.info, locus.params, mig.rates = NULL, hist.ev = NULL, label = NULL) {
  opt <- options(scipen = 999)
  
  ploidy <- attr(locus.params, "ploidy")
  pop.info[, c("pop.size", "sample.size")] <- pop.info[, c("pop.size", "sample.size")] * ploidy
  
  if(is.null(label)) label <- "fastsimcoal.output"
  file <- paste(label, ".par", sep = "")
  mig.rates <- if(!is.null(mig.rates)) if(is.list(mig.rates)) mig.rates else list(mig.rates)
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)
  
  # Write input file
  write(paste("//  <<", label, ">>  (input from 'fastsimcoal.skeleSim.run')"), file)
  write(paste(nrow(pop.info), "populations to sample"), file, append = T)
  
  write("//Population effective sizes", file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, "pop.size"], file, append = T)
  
  write("//Sample sizes", file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, c("sample.size", "sample.times")], file, append = T)
  
  write("//Growth rates", file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, "growth.rate"], file, append = T)
  
  write("//Number of migration matrices", file, append = T)
  write(length(mig.rates), file, append = T)
  if(!is.null(mig.rates)) {
    for(i in 1:length(mig.rates)) {
      write("//migration matrix", file, append = T)
      for(r in 1:nrow(mig.rates[[i]])) write(mig.rates[[i]][r, ], file, append = T)
    }
  }
  
  write("//Historical events: time, source, sink, migrants, new size, growth rate, migr. matrix", file, append = T)
  write(ifelse(is.null(hist.ev), 0, nrow(hist.ev)), file, append = T)
  if(!is.null(hist.ev)) {
    for(i in 1:nrow(hist.ev)) {
      write(paste(hist.ev[i, ], collapse = " "), file, append = T)
    }
  }
  
  num.chrom <- attr(locus.params, "num.chrom")
  if(!is.null(num.chrom)) locus.params$chromosome <- 1
  locus.params <- split(locus.params, locus.params$chromosome)
  num.independent <- if(is.null(num.chrom)) length(locus.params) else num.chrom
  chrom.struct <- if(num.independent == 1 | !is.null(num.chrom)) 0 else 1
  write("//Number of independent loci [chromosomes]", file, append = T)
  write(paste(num.independent, chrom.struct), file, append = T)
  for(block in locus.params) {
    block$chromosome <- NULL
    write("//Per chromosome: Number of linkage blocks", file, append = T)
    write(nrow(block), file, append = T)
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file, append = T)
    for(i in 1:nrow(block)) {
      line <- paste(block[i, ], collapse = " ")
      write(gsub(" NA", "", line), file, append = T)
    }
  }
  
  options(opt)
  invisible(file)
}


.fscRead <- function(file, locus.params) {
  formatGenotypes <- function(x, ploidy) {
    # reformat matrix to have alleles side-by-side
    nloci <- ncol(x) - 2
    loc.end <- seq(ploidy, nrow(x), by = ploidy)
    gen.data <- do.call(rbind, lapply(loc.end, function(i) {
      allele.i <- (i - ploidy + 1):i
      loci <- as.vector(x[allele.i, -(1:2)])
      id <- paste(x[allele.i, 2], collapse = ".")
      pop <- x[allele.i[1], 1]
      c(id, pop, loci)
    }))
    # rename loci
    locus_names <- paste("Locus", zero.pad(1:nloci), sep = "_")
    locus_names <- paste(rep(locus_names, each = ploidy), 1:ploidy, sep = ".")
    colnames(gen.data) <- c("id", "pop", locus_names)
    gen.data
  }
  
  formatDNA <- function(dna.seq, pop, locus.params) {
    # create multidna object splitting chromosomes into loci
    num.chrom <- attr(locus.params, "num.chrom")
    chrom.pos <- if(is.null(num.chrom)) {
      tapply(locus.params$num.markers, locus.params$chromosome, sum)
    } else {
      rep(sum(locus.params$num.markers), num.chrom)
    }
    chrom.pos <- cumsum(chrom.pos)
    chrom.pos <- cbind(start = c(1, chrom.pos[-length(chrom.pos)] + 1), end = chrom.pos)
    
    rownames(dna.seq) <- pop
    dna.seq <- tolower(dna.seq)
    new("multidna", lapply(1:nrow(chrom.pos), function(i) {
      as.matrix(dna.seq)[, chrom.pos[i, "start"]:chrom.pos[i, "end"]]
    }))
  }
  
  f <- readLines(file)
  
  # get start and end points of data blocks
  start <- grep("SampleData=", f) + 1
  end <- which(f == "}") - 2
  pos <- cbind(start, end)
  
  # compile data into 3 column character matrix
  data.mat <- do.call(rbind, lapply(1:nrow(pos), function(i) {
    f.line <- f[pos[i, 1]:pos[i, 2]]
    f.line <- gsub("[[:space:]]+", "--", f.line)
    result <- do.call(rbind, strsplit(f.line, "--"))[, -2]
    cbind(rep(paste("Sample", i), nrow(result)), result)
  }))
  
  ploidy <- attr(locus.params, "ploidy")
  
  # get data type
  data.type <- f[grep("DataType=", f)]
  data.type <- gsub("\tDataType=", "", data.type)
  switch(
    data.type,
    DNA = { # diploid SNPs
      dna.seq <- do.call(rbind, strsplit(data.mat[, 3], ""))
      if(attr(locus.params, "ploidy") == 2) {
        gen.data <- formatGenotypes(cbind(data.mat[, 1:2], dna.seq), ploidy)
        df2gtypes(gen.data, ploidy, description = file)
      } else { # haploid DNA sequences
        dna.seq <- formatDNA(dna.seq, data.mat[, 2], locus.params)
        g <- sequence2gtypes(dna.seq, strata = data.mat[, 1], description = file)
        labelHaplotypes(g)$gtype
      }
    },
    MICROSAT = {
      gen.data <- formatGenotypes(data.mat, ploidy)
      df2gtypes(gen.data, ploidy, description = file)
    },
    NULL
  )
}


#' @rdname fastsimcoal
#' @export
#' 
fastsimcoal <- function(pop.info, locus.params, mig.rates = NULL, 
                        hist.ev = NULL, label = NULL, quiet = TRUE, 
                        delete.files = TRUE, exec = "fsc252", num.cores = NULL) {
  
  if(is.null(label)) label <- "fsc.run"
  if(file.exists(label)) for(f in dir(label, full.names = T)) file.remove(f)
  
  # Write fastsimcoal input file
  infile <- .fscWrite(
    pop.info = pop.info, locus.params = locus.params,
    mig.rates = mig.rates, hist.ev = hist.ev, label = label
  )
  
  # Run fastsimcoal
  cores.spec <- if(!is.null(num.cores)) {
    num.cores <- max(1, num.cores)
    num.cores <- min(num.cores, min(detectCores(), 12))
    paste(c("-c", "-B"), num.cores, collapse = " ")
  } else ""
  cmd.line <- paste(
    exec, "-i", infile, "-n 1",
    ifelse(quiet, "-q", ""), "-S", cores.spec,
    attr(locus.params, "opts")
  )
  err <- if(.Platform$OS.type == "unix") {
    system(cmd.line, intern = F)
  } else {
    shell(cmd.line, intern = F)
  }
  
  if(err == 0) {
    if(!quiet) cat("fastsimcoal exited normally\n")
  } else {
    stop("fastsimcoal exited with error ", err, "\n")
  }
  
  # Read and parse output
  arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
  if(!file.exists(arp.file)) stop("fastsimcoal did not generate output")
  g <- .fscRead(arp.file, locus.params)
  
  # Cleanup
  if(delete.files) {
    unlink(label, recursive = TRUE, force = TRUE)
    file.remove(infile)
  }
  
  return(g)
}