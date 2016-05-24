#' @name structure
#' @title STRUCTURE
#' @description Run STRUCTURE to assess group membership of samples.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param k.range vector of values to for \code{maxpop} in multiple runs. 
#'   If set to \code{NULL}, a single STRUCTURE run is conducted with 
#'   \code{maxpops} groups. If specified, do not also specify \code{maxpops}.
#' @param num.k.rep number of replicates for each value in \code{k.range}.
#' @param label label to use for input and output files
#' @param delete.files logical. Delete all files when STRUCTURE is finished?
#' @param exec name of executable for STRUCTURE. Defaults to "structure".
#' @param ... arguments to be passed to \code{structure.write}.
#' @param maxpops number of groups.
#' @param burnin number of iterations for MCMC burnin.
#' @param numreps number of MCMC replicates.
#' @param noadmix logical. No admixture?
#' @param freqscorr logical. Correlated frequencies?
#' @param randomize randomize.
#' @param seed set random seed.
#' @param pop.prior a character specifying which population prior model to 
#'   use: "locprior" or "usepopinfo".
#' @param locpriorinit parameterizes locprior parameter \emph{r} - how 
#'   informative the populations are. Only used when 
#'   \code{pop.prior} = "locprior".
#' @param maxlocprior specifies range of locprior parameter \emph{r}. Only used 
#'   when \code{pop.prior} = "locprior".
#' @param gensback integer defining the number of generations back to test 
#'   for immigrant ancestry. Only used when \code{pop.prior} = "usepopinfo".
#' @param migrprior numeric between 0 and 1 listing migration prior. Only used 
#'   when \code{pop.prior} = "usepopinfo".
#' @param pfrompopflagonly logical. update allele frequencies from individuals 
#'   specified by \code{popflag}. Only used when \code{pop.prior} = 
#'   "usepopinfo".
#' @param popflag a vector of integers (0, 1) or logicals identifiying whether 
#'   or not to use strata information. Only used when \code{pop.prior} 
#'   = "usepopinfo".
#' @param file name of the output file from STRUCTURE.
#' @param pops vector of population labels to be used in place of numbers in 
#'   STRUCTURE file.
#'    
#' @return
#' \describe{
#'  \item{structure.run}{a list where each element is a list with results 
#'    from \code{structure.read} and a vector of the filenames used.}
#'  \item{structure.write}{a vector of the filenames used by STRUCTURE.}
#'  \item{structure.read}{a list containing:
#'    \tabular{ll}{
#'      \code{summary} \tab new locus name, which is a combination of loci 
#'        in group.\cr
#'      \code{q.mat} \tab data.frame of assignment probabilities for 
#'        each id.\cr
#'      \code{prior.anc} \tab list of prior ancestry estimates for each 
#'        individual where population priors were used.\cr
#'      \code{files} \tab vector of input and output files used by STRUCTURE.\cr
#'      \code{label} \tab label for the run.\cr
#'    }
#'  }
#' }
#' 
#' @note STRUCTURE is not included with \code{strataG} and must be downloaded 
#'   separately. Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs} for installation instructions.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Pritchard, J.K., M. Stephens, P. Donnelly. 2000. Inference of 
#'   population structure using multilocus genotype data. Genetics 155:945-959.\cr 
#'   \url{http://pritchardlab.stanford.edu/structure.html}
#' 
#' @seealso \code{\link{structurePlot}}, \code{\link{evanno}}, 
#'   \code{\link{clumpp}} 
#' 
#' @examples
#' \dontrun{
#' data(msats.g)
#' 
#' # Run STRUCTURE
#' sr <- structureRun(msats, k.range = 1:4, num.k.rep = 10)
#' 
#' # Calculate Evanno metrics
#' evno <- evanno(sr)
#' evno
#' 
#' # Run CLUMPP to combine runs for K = 2
#' q.mat <- clumpp(sr, k = 3)
#' q.mat
#' 
#' # Plot CLUMPP results
#' structurePlot(q.mat)
#' }
#' 
#' @importFrom utils file_test
#' @export
#' 
structureRun <- function(g, k.range = NULL, num.k.rep = 1, label = NULL, 
                         delete.files = TRUE, exec = "structure", ...) {
  
  if(!ploidy(g) > 1) stop("'g' must have a ploidy > 1")
  
  # setup folder
  if(is.null(label)) label <- paste(description(g), "structureRun", sep = ".")
  label <- gsub("[[:punct:]|[:space:]]", ".", label)
  unlink(label, recursive = TRUE, force = TRUE)
  dir.create(label)
  if(!file_test("-d", label)) {
    stop(paste("'", label, "' is not a valid folder.", sep = ""))
  }
  label <- file.path(label, label)
  
  # setup k and replicate data.frame to cycle through
  if(is.null(k.range)) k.range <- 1:nlevels(strata(g))
  rep.df <- expand.grid(rep = 1:num.k.rep, k = k.range)
  
  rownames(rep.df) <- paste(label, ".k", rep.df$k, ".r", rep.df$rep, sep = "")
  out.files <- lapply(rownames(rep.df), function(x) {
    sw.out <- structureWrite(g, label = x, maxpops = rep.df[x, "k"], ...)
    files <- sw.out$files
    cmd <- paste(exec, " -m ", files["mainparams"], 
                 " -e ", files["extraparams"], 
                 " -i ", files["data"], 
                 " -o ", files["out"], 
                 sep = ""
    )
    
    err.code <- system(cmd)
    if(err.code == 127) {
      stop("You do not have STRUCTURE installed.")
    } else if(!err.code == 0) {
      stop(paste("Error running STRUCTURE. Error code", err.code, "returned."))
    }
    
    files["out"] <- paste(files["out"], "_f", sep = "")
    result <- structureRead(files["out"], sw.out$pops)
    
    if(file.exists("seed.txt")) file.remove("seed.txt")
    files <- if(delete.files) NULL else files
    
    result <- c(result, list(files = files, label = basename(x)))
    fname <- paste(x, ".ws.rdata", sep = "")
    save(result, file = fname)
    fname
  })
  
  run.result <- lapply(out.files, function(f) {
    result <- NULL
    load(f)
    result
  })
  names(run.result) <- sapply(run.result, function(x) x$label)
  class(run.result) <- c("structure.result", class(run.result))
  
  if(delete.files) unlink(dirname(label), recursive = TRUE, force = TRUE)
  run.result
}


#' @rdname structure
#' @export
#' 
structureWrite <- function(g, label = NULL, maxpops = nlevels(strata(g)), 
                           burnin = 1000, numreps = 1000, noadmix = TRUE, 
                           freqscorr = FALSE, randomize = TRUE, seed = 0, 
                           pop.prior = NULL, locpriorinit = 1, maxlocprior = 20, 
                           gensback = 2, migrprior = 0.05,
                           pfrompopflagonly = TRUE, popflag = NULL, ...) {
  
  if(ploidy(g) != 2) stop("'g' must be diploid")
  
  # check parameters
  if(!is.null(pop.prior)) {
    if(!pop.prior %in% c("locprior", "usepopinfo")) {
      stop("'pop.prior' must be 'locprior' or 'usepopinfo'.")
    }
  }
  if(is.null(popflag)) popflag <- rep(1, nInd(g))
  popflag <- as.numeric(popflag)
  if(length(popflag) != nInd(g)) {
    stop("'popflag' should be the same length as the number of individuals in 'g'.")
  }
  if(!all(popflag %in% c(0, 1))) {
    stop("all values in 'popflag' must be 0 or 1.")
  }
  
  in.file <- ifelse(is.null(label), "data", paste(label, "data", sep = "_"))
  out.file <- ifelse(is.null(label), "out", paste(label, "out", sep = "_"))
  main.file <- ifelse(is.null(label), "mainparams", 
                      paste(label, "mainparams", sep = "_"))
  extra.file <- ifelse(is.null(label), "extraparams", 
                       paste(label, "extraparams", sep = "_"))
  
  # write data
  write(paste(locNames(g), collapse = " "), file = in.file)
  popdata <- as.numeric(strata(g))
  
  for(i in 1:nInd(g)) {
    id <- indNames(g)[i]
    loci <- loci(g, id, locNames(g))
    loci <- c(as.matrix(sapply(loci, as.numeric)))
    loci[is.na(loci)] <- -9
    loci <- paste(loci, collapse = " ")
    write(paste(id, popdata[i], popflag[i], loci), 
          file = in.file, append = TRUE
    )
  }
  
  # write mainparams
  main.params <- c(
    paste("MAXPOPS", as.integer(maxpops)),
    paste("BURNIN", as.integer(burnin)),
    paste("NUMREPS", as.integer(numreps)),
    paste("INFILE", in.file),
    paste("OUTFILE", out.file),
    paste("NUMINDS", nInd(g)),
    paste("NUMLOCI", nLoc(g)),
    "MISSING -9",
    "ONEROWPERIND 1",
    "LABEL 1",
    "POPDATA 1",
    "POPFLAG 1",
    "LOCDATA 0",
    "PHENOTYPE 0",
    "EXTRACOLS 0",
    "MARKERNAMES 1"
  )
  main.params <- paste("#define", main.params)
  write(main.params, file = main.file)
  
  # write extraparams
  extra.params <- c(
    paste("NOADMIX", as.integer(noadmix)),
    paste("FREQSCORR", as.integer(freqscorr)),
    "INFERALPHA 1",
    "ALPHA 1.0",
    "FPRIORMEAN 0.01",
    "FPRIORSD 0.05",
    "LAMBDA 1.0",
    "UNIFPRIORALPHA 1", 
    "ALPHAMAX 20.0",
    "ALPHAPRIORA 0.05",
    "ALPHAPRIORB 0.001",
    "COMPUTEPROB 1",
    paste("ADMBURNIN", max(0, as.integer(burnin / 2))),
    "ALPHAPROPSD 0.025",
    "STARTATPOPINFO 0",
    paste("RANDOMIZE", as.integer(randomize)),
    paste("SEED", as.integer(seed)),
    "METROFREQ 10",
    "REPORTHITRATE 0" 
  )
  
  if(!is.null(pop.prior)) {
    pop.prior <- tolower(pop.prior)
    prior.params <- if(pop.prior == "locprior") {
      c("LOCPRIOR 1",
        "LOCISPOP 1",
        paste("LOCPRIORINIT", locpriorinit),
        paste("MAXLOCPRIOR", maxlocprior)
      )        
    } else if(pop.prior == "usepopinfo") {
      c("USEPOPINFO 1",
        paste("GENSBACK", trunc(gensback)),
        paste("MIGRPRIOR", migrprior),
        paste("PFROMPOPFLAGONLY", as.integer(pfrompopflagonly))
      )
    }
    extra.params <- c(extra.params, prior.params)
  }
  
  extra.params <- extra.params[!is.na(extra.params)]
  extra.params <- paste("#define", extra.params)
  write(extra.params, file = extra.file)
  
  invisible(list(
    files = c(data = in.file, mainparams = main.file, 
              extraparams = extra.file, out = out.file),
    pops = levels(strata(g))
  ))
}


#' @rdname structure
#' @export
#' 
structureRead <- function(file, pops = NULL) {
  if(!file.exists(file)) {
    stop(paste("the file '", file, "' can't be found.", sep = ""))
  }
  
  # Read file to get single results and parameter values
  result <- scan(file, "character", quiet = TRUE)
  
  loc <- grep("Estimated", result, ignore.case = FALSE, value = FALSE)
  est.ln.prob <- as.numeric(result[loc[1] + 6])
  
  loc <- grep("likelihood", result, ignore.case = FALSE, value = FALSE)
  mean.lnL <- as.numeric(result[loc[1] + 2])
  var.lnL <- as.numeric(result[loc[2] + 2])
  
  loc <- grep("MAXPOPS", result, value = F)
  maxpops <- result[loc]
  maxpops <- sub("MAXPOPS=", "", maxpops)
  maxpops <- as.integer(sub(",", "", maxpops))
  
  loc <- grep("GENSBACK", result, value = F)
  gensback <- result[loc]
  gensback <- sub("GENSBACK=", "", gensback)
  gensback <- as.integer(sub(",", "", gensback))
  
  smry <- c(k = maxpops, est.ln.prob = est.ln.prob, mean.lnL = mean.lnL, var.lnL = var.lnL)
  
  # Read file to get population assignment probability table
  result <- scan(file, "character", sep = "\n", quiet = TRUE)
  first <- grep("(%Miss)", result, value = FALSE) + 1
  last <- grep("Estimated Allele", result, value = FALSE) - 1
  tbl.txt <- result[first:last]
  
  # Remove special characters from table and whitespace from end of lines
  tbl.txt <- sub("[*]+", "", tbl.txt)
  tbl.txt <- sub("[(]", "", tbl.txt)
  tbl.txt <- sub("[)]", "", tbl.txt)
  tbl.txt <- sub("[|][ ]+$", "", tbl.txt)
  
  # Find which lines have had population priors
  prior.lines <- grep("[|]", tbl.txt)
  
  # Create table of population assignments for lines without population priors
  no.prior <- if(length(prior.lines) < length(tbl.txt)) {
    no.prior.q.txt <- if(length(prior.lines) == 0) tbl.txt else tbl.txt[-prior.lines]
    .structureParseQmat(no.prior.q.txt, pops)
  } else NULL
  
  # Return just this base table if MAXPOPS = 1
  if(maxpops == 1) {
    no.prior$row <- NULL
    return(list(summary = smry, q.mat = no.prior, prior.anc = NULL))
  }
  
  # Create table of population assignments for lines with population priors
  has.prior <- if(length(prior.lines) > 0) {
    prior.txt <- strsplit(tbl.txt[prior.lines], "[|]")
    # Get base table
    prior.q.txt <- unlist(lapply(prior.txt, function(x) x[1]))
    df <- .structureParseQmat(prior.q.txt, pops)
    # Parse ancestry assignments into matrix
    prior.anc <- lapply(prior.txt, function(x) {
      anc.mat <- matrix(NA, nrow = maxpops, ncol = gensback + 1)
      rownames(anc.mat) <- paste("Pop", 1:nrow(anc.mat), sep = ".")
      colnames(anc.mat) <- paste("Gen", 0:gensback, sep = ".")
      # Split on whitespace and colons
      x <- sapply(strsplit(x[-1], "\\s|[:]"), function(y) {
        y <- y[y != ""] # remove empty strings
        y[-1] # return vector with first element ("Pop") removed - vector has population # and ancestry assignments
      })
      # Populate ancestry matrix with probabilities
      for(i in 1:ncol(x)) {
        pop <- as.numeric(x[1, i])
        anc.mat[pop, ] <- as.numeric(x[-1, i])
      }
      anc.mat
    })
    names(prior.anc) <- df$id
    
    # Create population probability matrix for samples with priors and add to end of base table
    prob.mat <- t(sapply(1:nrow(df), function(i) {
      pop.probs <- rowSums(prior.anc[[i]])
      pop.probs[is.na(pop.probs)] <- df$prob.1[i]
      pop.probs
    }))
    colnames(prob.mat) <- paste("prob", 1:ncol(prob.mat), sep = ".")
    df$prob.1 <- NULL
    df <- cbind(df, prob.mat)
    
    list(df = df, prior.anc = prior.anc)
  } else NULL
  
  # Combine assignment probability matrices
  has.prior.df <- if(is.null(has.prior)) NULL else has.prior$df
  q.mat <- rbind(no.prior, has.prior.df)
  q.mat <- q.mat[order(q.mat$row), ]
  q.mat$row <- NULL
  rownames(q.mat) <- NULL
  # Make sure all probs sum to 1
  q.mat[, -(1:3)] <- t(apply(q.mat[, -(1:3)], 1, function(i) i / sum(i)))
  
  prior.anc <- if(is.null(has.prior)) NULL else has.prior$prior.anc
  
  list(summary = smry, q.mat = q.mat, prior.anc = prior.anc)
}


# Internal file used by 'stucture' and 'clumpp' to parse output files
# q.mat.txt: character vector of Q matrix from STRUCTURE output file.
# pops: vector of population labels to be used in place of numbers in STRUCTURE file.
#   
.structureParseQmat <- function(q.mat.txt, pops) {
  q.mat.txt <- sub("[*]+", "", q.mat.txt)
  q.mat.txt <- sub("[(]", "", q.mat.txt)
  q.mat.txt <- sub("[)]", "", q.mat.txt)
  q.mat.txt <- sub("[|][ ]+$", "", q.mat.txt)
  
  # Parse population assignment portion of table to create single-line 
  #   data.frame
  do.call(rbind, lapply(q.mat.txt, function(x) {
    # Split on spaces and remove empty spaces and colons
    x <- strsplit(x, " ")[[1]]
    x <- x[!x %in% c("", ":")]
    p <- as.numeric(x[4])
    
    df <- data.frame(row = as.numeric(x[1]), id = x[2], 
                     pct.miss = as.numeric(x[3]), 
                     orig.pop = if(is.null(pops)) p else pops[p], 
                     stringsAsFactors = FALSE
    )
    pop.prob <- as.data.frame(rbind(as.numeric(x[-(1:4)])))
    colnames(pop.prob) <- paste("Group", 1:ncol(pop.prob), sep = ".")
    cbind(df, pop.prob) 
  }))
}