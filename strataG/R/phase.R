#' @name phase 
#' @title PHASE
#' @description Run PHASE to estimate the phase of loci in diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param loci vector or data.frame of loci in 'g' that are to be phased. If a 
#'   data.frame, it should have columns named \code{locus} (name of locus in 'g'),  
#'   \code{group} (number identifying loci in same linkage group), and
#'   \code{position} (integer identifying location of each locus in a linkage group).
#' @param positions position along chromosome of each locus.
#' @param type type of each locus.
#' @param num.iter number of PHASE MCMC iterations.
#' @param thinning number of PHASE MCMC iterations to thin by.
#' @param burnin number of PHASE MCMC iterations for burnin.
#' @param model PHASE model type.
#' @param ran.seed PHASE random number seed.
#' @param final.run.factor optional.
#' @param save.posterior logical. Save posterior sample in output list?
#' @param in.file name to use for PHASE input file.
#' @param out.file name to use for PHASE output files.
#' @param delete.files logical. Delete PHASE input and output files when done?
#' @param ph.res result from \code{phase.run}.
#' @param thresh minimum probability for a genotype to be selected (0.5 - 1).
#' @param keep.missing logical. T = keep missing data from original data set. F = Use estimated genotypes from PHASE.
#'  
#' @note PHASE is not included with \code{strataG} and must be downloaded 
#'   separately. Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs} for installation instructions. 
#' 
#' @details
#' \tabular{ll}{
#'   \code{phase} \tab runs PHASE assuming that the executable is installed properly and available on the command line.\cr
#'   \code{phaseWrite} \tab writes a PHASE formatted file.\cr
#'   \code{phaseReadPair} \tab reads the '_pair' output file.\cr
#'   \code{phaseReadSample} \tab reads the '_sample' output file.\cr
#'   \code{phaseFilter} \tab filters the result from \code{phase.run} to extract one genotype for each sample.\cr
#'   \code{phasePosterior} \tab create a data.frame all genotypes for each posterior sample.\cr
#' }
#'  
#' @return
#' \describe{
#'  \item{phase}{a list containing:
#'    \tabular{ll}{
#'      \code{locus.name} \tab new locus name, which is a combination of loci in group.\cr
#'      \code{gtype.probs} \tab a data.frame listing the estimated genotype for every sample along with probability.\cr
#'      \code{orig.gtypes} \tab the original gtypes object for the composite loci.\cr
#'      \code{posterior} \tab a list of \code{num.iter} data.frames representing posterior sample of genotypes for each sample.\cr
#'    }}
#'  \item{phaseWrite}{a list with the input filename and the \linkS4class{gtypes} object used.}
#'  \item{phaseReadPair}{a data.frame of genotype probabilities.}
#'  \item{phaseReadSample}{a list of data.frames representing the posterior sample of genotypes for one set of loci for each sample.}
#'  \item{phaseFilter}{a matrix of genotypes for each sample.}
#'  \item{phasePosterior}{a list of data.frames representing the posterior sample of all genotypes for each sample.}
#' }
#' 
#' @references Stephens, M., and Donnelly, P. (2003). A comparison of Bayesian methods for haplotype reconstruction from 
#'   population genotype data. American Journal of Human Genetics 73:1162-1169.
#'   Available at: \url{http://stephenslab.uchicago.edu/software.html#phase}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples \dontrun{
#' data(bowhead.snps)
#' data(bowhead.snp.position)
#' snps <- df2gtypes(bowhead.snps, ploidy = 2, description = "Bowhead SNPS")
#' summary(snps)
#' 
#' # Run PHASE on all data
#' phase.results <- phase(snps, bowhead.snp.position, num.iter = 100, 
#'   save.posterior = FALSE)
#' 
#' # Filter phase results
#' filtered.results <- phaseFilter(phase.results, thresh = 0.5)
#' 
#' # Convert phased genotypes to gtypes
#' ids <- rownames(filtered.results)
#' strata <- bowhead.snps$Stock[match(ids, bowhead.snps$LABID)]
#' filtered.df <- cbind(id = ids, strata = strata, filtered.results)
#' phased.snps <- df2gtypes(filtered.df, ploidy = 2, description = "Bowhead phased SNPs")
#' summary(phased.snps)
#' }
#' 
#' @export
#' 
phase <- function(g, loci, positions = NULL, type = NULL,
  num.iter = 100000, thinning = 100, burnin = 100000, model = "new", 
  ran.seed = NULL, final.run.factor = NULL, save.posterior = FALSE, 
  in.file = "phase_in", out.file = "phase_out", delete.files = TRUE) { 
  
  if(ploidy(g) != 2) stop("'g' must be diploid")
  
  # check loci format
  if(!is.data.frame(loci)) {
    if(!(is.character(loci) & is.vector(loci))) {
      stop("'loci' must be a data.frame or character vector")
    }
    if(is.null(positions)) positions <- rep(1, nLoc(g))
    if(length(positions) != length(loci)) {
      stop("'positions' must be same length as 'loci'")
    }
    loci <- data.frame(locus = loci, position = positions, group = 1)
  }
  loci$group <- as.character(loci$group)
  loci$position <- as.numeric(loci$position)
  
  if(is.null(type)) type <- rep("S", length(unique(loci$group)))
  if(length(type) != length(unique(loci$group))) {
    stop("'type' must be same length as number of locus groups")
  }
  names(type) <- unique(loci$group)
  
  result <- lapply(unique(loci$group), function(grp) {
    lets <- paste(sample(c(0:9, letters), 10, replace = TRUE), collapse = "")
    in.file <- paste("phase_in_", lets, sep = "")
    out.file <- paste("phase_out_", lets, sep = "")
    
    # Write input file
    group.df <- loci[loci$group == grp, ]  
    locus.type <- rep(type[grp], nrow(group.df))
    in.file.data <- phaseWrite(g, loci = group.df$locus, 
                               positions = group.df$position, 
                               type = locus.type, in.file)

    # Set parameters
    M.opt <- switch(model, new = "-MR", old = "-MS", hybrid = "-MQ", "")           
    S.opt <- ifelse(is.null(ran.seed), "", paste("-S", ran.seed, sep = ""))
    X.opt <- ifelse(is.null(final.run.factor), "", 
                    paste("-X", final.run.factor, sep = ""))
    s.opt <- ifelse(save.posterior, "-s", "")
    in.file.opt <- paste("\"", in.file, "\"", sep = "")
    out.file.opt <- paste("\"", out.file, "\"", sep = "")
    iter.params <- paste(trunc(num.iter), trunc(thinning), trunc(burnin))
    phase.cmd <- paste("PHASE", M.opt, S.opt, X.opt, s.opt, in.file.opt, 
                       out.file.opt, iter.params)
    
    # Run Phase
    err.code <- system(phase.cmd)  
    if(err.code == 127) {
      stop("You do not have PHASE installed.") 
    } else if(!err.code == 0) {
      stop(paste("Error running PHASE. Error code", err.code, "returned."))
      cat("\n")
    }
    
    # Read output
    gtype.probs <- phaseReadPair(paste(out.file, "_pairs", sep = ""))
    if(is.null(gtype.probs)) {
      alleles <- rep(NA, nrow(g$genotypes))
      gtype.probs <- data.frame(
        id = indNames(g), a1 = alleles, a2 = alleles, pr = rep(1, nInd(g))
      )
    }
    new.locus.name <- paste(group.df$locus, collapse = "_")
    alleles <- paste(new.locus.name, 1:2, sep = ".")
    colnames(gtype.probs)[1:3] <- c("id", alleles) 
    rownames(gtype.probs) <- NULL
    
    locus.result <- list(
      locus.name = new.locus.name, gtype.probs = gtype.probs, 
      orig.gtypes = in.file.data$gtypes
    )  
    
    if(save.posterior) {
      file <- paste(out.file, "_sample", sep = "")
      l.type <- paste(locus.type, collapse = "")
      locus.result$posterior <- phaseReadSample(file, l.type)
      for(i in 1:length(locus.result$posterior)) {
        colnames(locus.result$posterior[[i]]) <- c("id", alleles)
      }
    }
    
    if(delete.files) {
      file.remove(c(dir(pattern = in.file), dir(pattern = out.file)))
    }
    
    locus.result  
  })
  
  names(result) <- lapply(result, function(x) x$locus.name)
  class(result) <- c("phase.result", class(result))
  result
}


#' @rdname phase
#' @export
#' 
phaseReadSample <- function(out.file, type) {
  if(!file.exists(out.file)) return(NULL)
  post.file <- scan(file = out.file, what = "character", 
                    sep = "\n", quiet = TRUE)
  iter.start <- grep(type, post.file) + 1
  lapply(iter.start, function(start) {
    num.samples <- as.integer(post.file[start - 3])
    end <- start + (num.samples * 3) - 3
    as.data.frame(t(sapply(seq(start, end, by = 3), function(i) {
      id <- strsplit(post.file[i], " ")[[1]][2]    
      hap1 <- gsub(" ", "", post.file[i + 1])
      hap2 <- gsub(" ", "", post.file[i + 2])
      c(id, hap1, hap2)
    })), stringsAsFactors = FALSE)
  })
}


#' @rdname phase
#' @export
#' 
phaseReadPair <- function(out.file) {   
  if(!file.exists(out.file)) return(NULL)
  pair.file <- scan(file = out.file, what = "character", 
                    sep = "\n", quiet = TRUE)
  
  id.start <- grep("IND:", pair.file)
  gtype.probs <- lapply(1:length(id.start), function(i) {
    id.end <- ifelse(i == length(id.start), 
                     length(pair.file), 
                     id.start[i + 1] - 1)
    id <- sub("IND: ", "", pair.file[id.start[i]])
    t(sapply((id.start[i] + 1):id.end, function(j) {
      line.split <- unlist(strsplit(pair.file[j], " , "))
      names(line.split) <- c("hap1", "hap2", "pr")
      c(id = id, line.split)
    }))
  })             
  gtype.probs <- as.data.frame(do.call(rbind, gtype.probs), 
                               stringsAsFactors = FALSE)  
  gtype.probs$pr <- as.numeric(as.character(gtype.probs$pr))
  gtype.probs
}


#' @rdname phase
#' @export
#' 
phaseWrite <- function(g, loci, positions = NULL, 
                       type = rep("S", length(loci)), in.file = "phase_in") {
  
  if(ploidy(g) != 2) stop("'g' must be diploid")
  
  # Make sure locus.names and locus.positions are sorted properly
  if(is.null(positions)) positions <- rep(1, length(nLoc(g)))
  asc.order <- order(positions)
  loci <- loci[asc.order]
  positions <- positions[asc.order]
  
  sub.g <- g[, loci, ]
  write(c(
    nInd(sub.g), 
    length(loci),
    paste("P", paste(positions, collapse = " ")),
    paste(type, collapse = ""), ""
  ), file = in.file)
  
  g.mat <- as.matrix(sub.g)
  g.mat[is.na(g.mat)] <- "?"
  for(i in 1:nrow(g.mat)) {
    write(c(
      rownames(g.mat)[i],
      paste(g.mat[i, seq(1, ncol(g.mat) - 1, 2)], collapse = " "),
      paste(g.mat[i, seq(2, ncol(g.mat), 2)], collapse = " ")
    ), file = in.file, append = TRUE)
  }
  
  invisible(list(filename = in.file, gtypes = sub.g))
}


#' @rdname phase
#' @export
#' 
phasePosterior <- function(ph.res, keep.missing = TRUE) {    
  if(!"phase.result" %in% class(ph.res)) {
    stop("'ph.res' is not a result from 'phase.run'.")
  }
  
  num.iter <- length(ph.res[[1]]$posterior)
  lapply(1:num.iter, function(iter) {
    do.call(cbind, lapply(1:length(ph.res), function(locus) {
      ph.res <- ph.res[[locus]]
      post.df <- ph.res$posterior[[iter]]
      
      if(keep.missing) {
        for(i in 1:nrow(post.df)) {
          ids <- which(indNames(ph.res$orig.gtypes) == post.df[i, 1])
          if(any(is.na(loci(ph.res$orig.gtypes[ids, , ])))) {
            post.df[i, 2:3] <- NA
          }
        }
      }
      
      cols <- if(locus == 1) {1:3} else {2:3}
      post.df[, cols]
    })) 
  })
}


#' @rdname phase
#' @export
#' 
phaseFilter <- function(ph.res, thresh = 0.5, keep.missing = TRUE) {
  if(!"phase.result" %in% class(ph.res)) {
    stop("'ph.res' is not a result from 'phase.run'.")
  }
  
  filtered <- lapply(ph.res, function(x) {
    gtype.probs <- x$gtype.probs
    pr.vec <- unique(gtype.probs[, 1])
    locus.filtered <- do.call(rbind, lapply(pr.vec, function(i) {
      this.id <- gtype.probs[gtype.probs[, 1] == i, ]
      max.index <- which.max(this.id$pr)
      if(length(max.index) == 0) return(this.id[1, ])
      kept.line <- this.id[max.index, ]
      if(as.numeric(kept.line$pr) < thresh) kept.line[, 2:3] <- c(NA, NA)
      kept.line
    })) 
    rownames(locus.filtered) <- NULL
    
    if(keep.missing) {
      for(i in 1:nrow(locus.filtered)) {
        ids <- which(indNames(x$orig.gtypes) == locus.filtered[i, 1])
        if(any(is.na(loci(x$orig.gtypes[ids, , ])))) {
          locus.filtered[i, 2:3] <- NA
        }
      }
    }
    
    locus.filtered
  })
  
  id <- filtered[[1]][, 1]
  filtered <- as.matrix(do.call(cbind, lapply(filtered, function(x) x[, 2:3])))
  rownames(filtered) <- id
  colnames(filtered) <- paste(rep(names(ph.res), each = 2), ".", 1:2, sep = "")
  filtered
}