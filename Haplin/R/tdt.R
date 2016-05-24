######################################################################################
#################### PUBLIC FUNCTIONS:  ##############################################
######################################################################################

#### haplinTDT
#### Transmission disequilibrium tests (tdt, hhrr, trimm)
## Input parameters
##  - filename: input file (haplin format) containing trio data
##              Assume furthermore the columns: 'orig.lines' (and 'sex' for an x chromsome)
##  - nsim.perm: Number of simulations used in permutation test
##  - select.gender: Do the analysis for a gender subset. Values: 1, 2, or NULL.
##                   1: Male, 2: Female, NULL: All
##  - method: character vector containing which of the methods tdt, hhrr, trimm to be used.
##  - names.marker: Names given to the markers. If NULL, markers are named 1,2,...
##  - use.haplotypes: If FALSE, markers are analyzed individually. If TRUE, haplotypes are
##                    reconstructed by Haplin and then analyzed as one multi-allelic marker.
##  - use.ambiguous: TRUE if ambiguous trios are to be kept
### - remaining arguments: arguments passed on to haplin
haplinTDT <- function(filename,
                      nsim.perm=0, select.gender=NULL, method=c("tdt", "hhrr", "trimm"),
                      names.marker=NULL, use.haplotypes=FALSE, use.ambiguous=TRUE,
                      design = "triad", markers="ALL", n.vars=0, sep= " ", allele.sep= ";",
                      na.strings="NA", use.missing = FALSE,
                      xchrom = FALSE, sex = NULL, threshold = 0.01, verbose = TRUE, printout = TRUE)
{
  .mcall <- lapply(match.call()[-1], function(x) eval.parent(x, 3))
  .defaults.haplin <- formals(haplin)
  .info <- f.check.pars(.mcall, .defaults.haplin)

  if (!file.exists(filename)) {
    stop(paste("File", filename, "does not exist"))
  }
  
  if (!all(method %in% c("tdt", "hhrr", "trimm"))) {
    stop("method should be one of \"tdt\", \"hhrr\", \"trimm\"")
  }

  if (!is.null(select.gender)) {
    if (!all(select.gender %in% c(1, 2))) {
      stop("select.gender should be  1 or 2 if specified")
    }
  }
  
  if (.info$model$use.missing && !use.ambiguous) {
    cat("Set use.ambiguous = TRUE when use.missing = TRUE\n")
    use.ambiguous <- TRUE
  }

  if (nsim.perm < 0) {
    stop("nsim.perm must be >= 0")
  }

  if (.info$model$design != "triad") {
    stop("Only \"triad\" design is allowed\n")
  }
 
  .args <- .mcall

  ## Remove arguments that are not used in call to haplin
  .args$nsim.perm <- NULL
  .args$select.gender <- NULL
  .args$method <- NULL
  .args$names.marker <- NULL
  .args$use.haplotypes <- NULL
  .args$use.ambiguous <- NULL

  .args$reference <- 1 # To avoid message from haplin
  .args$data.out <- "prelim" # Return after data read

  if (use.haplotypes) {
    cat("Construct haplotypes... Remember: SNPs must be in correct order!\n")
    trios <- try(do.call("haplin", .args))
    trios <- cbind(trios, marker=rep(1, nrow(trios))) # Only 1 marker = the haplotype
  }
  else {
    ## First: Find number of markers in file
    if (.info$control$verbose) 
      cat("\nReading data from file...  ")
    if ((.info$model$design == "triad") | (.info$model$design == 
           "cc.triad")) {
      .fam <- "mfc"
    }
    if (.info$model$design == "cc") 
      .fam <- "c"
	.data.read <- f.read.data(info = .info)
    if (.info$control$verbose) 
      cat("Done\n")
	.markers <- attr(.data.read, "info")$filespecs$markers
    .nmarker <- length(.markers)

    ## Second: Read and transform data using haplin
    trios <- NULL
    for (i in seq(length.out = .nmarker)) {
      .args$markers <- .markers[i]
      trios0 <- try(do.call("haplin", .args))

      if (!inherits(trios0, "try-error"))
        trios <- rbind(trios, cbind(trios0, marker=rep(i, nrow(trios0))))
    }
  }
  
  object <- tdt(trios=trios, nsim.perm=nsim.perm, select.gender=select.gender, method=method,
                marker.names=names.marker, xchrom=.info$model$xchrom, use.ambiguous=use.ambiguous)
  object$select.gender <- select.gender
  object$use.haplotypes <- use.haplotypes
  object$xchrom <- .info$model$xchrom
  object$nsim.perm <- nsim.perm
  object$use.ambiguous <- use.ambiguous
  object$use.missing <- .info$model$use.missing
  object$orig.call <- sys.call()
  class(object) <- "haplinTDT"

  if (.info$control$printout) {
    plot(object)
    if (.info$control$verbose) 
      cat("\n#################################\n")
    print(summary(object))
  }

  object
}

print.haplinTDT <- function(x, ...) {
  cat("haplintTDT object containing methods: ")
 
  methods <- c("tdt", "hhrr", "trimm")
  ind <- methods %in% names(x)
  for (m in methods[ind]) {
    cat(m, "")
  }
  cat("\nSee summary() for detailed information.\n")
}

print.tdt <- function(x, ...) {
  ##  cat("TDT object----------\n")
  print.default(x)
}

print.summary.haplinTDT <- function(x, ...) {
  cat("haplinTDT object\n")
  if (x$use.haplotypes) {
    cat(" - Use haplotypes constructed by haplin\n")
  }
  else {
    cat(" - Single marker analyses\n")
  }

  if (x$xchrom) {
    cat(" - X chromosome data\n")
  }
  else {
    cat(" - Autosomal chromosome data\n")
  }

  if (x$use.ambiguous) {
    cat(" - Keep ambiguous trios\n")
  }
  else {
    cat(" - Remove ambiguous trios\n")
  }

  if (!x$use.missing) {
    cat(" - Remove missing values\n")
  }
  else {
    cat(" - Keep trios with missing values (reconstruct data)\n")
  }
  
  if (x$nsim.perm > 0) {
    cat(" - Number of simulations used in permutation tests: ", x$nsim.perm, "\n")
  }

  if (!is.null(x$select.gender)) {
    gender.txt <- c("male", "female")
    cat(" - Selected gender: ", gender.txt[x$select.gender], "\n")
  }


  methods <- c("tdt", "hhrr", "trimm")
  ind <- methods %in% names(x)
  for (m in methods[ind]) {
    cat("\n-------- Method ", m, "---------\n")
    print(summary(x[[m]]))
  }
}

print.summary.tdt <- function(x, ...) {
  p <- x$pval
  nmarker <- ncol(p)
  marker.names <- colnames(p)

  if (nmarker > 1 && !is.null(x$pval.max)) {
    cat("Global test:\n")
    cat("Maximum chi square: ", format(round(x$chisq.max,5), nsmall=5), "\n")
    cat("P-value: ", x$pval.max, " (permutation test)\n")
  }

  df <- rep(x$df, len=nmarker) ## x$df could either be a scalar or a vector of length nmarker
  
  cat("\nIndividual Marker scores:\n")
  for (i in 1:nmarker) {
    if (!is.null(marker.names)) {
      cat(paste("Marker ", marker.names[i], ":", sep=""), "\n")
    }
    else {
      cat(paste("Marker ", i, ":", sep=""), "\n")
    }
    if (!is.null(x$transmissionMatrix)) {
      cat("Transmission table:\n")
      print(round(x$transmissionMatrix[,,i], 2))
    }
    cat("Chi square: ", format(round(x$chisq[i],5), nsmall=5), "  Degrees of freedom:", df[i], "\n")
    cat("P-value:    ", format(round(p["chisq",i],5), nsmall=5), "(chisq test) ")
    if ("permutation" %in% rownames(p)) {
      cat(p["permutation",i], " (permutation test)")
    }
    cat("\n\n")
  }

}


summary.tdt <- function(object, ...) {
  class(object) <- "summary.tdt"
  return(object)
}

summary.haplinTDT <- function(object, ...) {
  class(object) <- "summary.haplinTDT"
  return(object)
}

plot.haplinTDT <- function(x, separate.plots = FALSE, filename, filetype="png", ask=TRUE, ...) {
  methods <- c("tdt", "hhrr", "trimm")
  ind <- methods %in% names(x)
  n <- sum(ind)

  savepar <- par(no.readonly=TRUE)
  on.exit(par(savepar))
  
  if (!missing(filename)) {
    .jpeg.size <- c(1024, 768)
    if (filetype == "png") {
      png(filename = filename, width = .jpeg.size[1], 
          height = .jpeg.size[2], pointsize = 12)
    }
    else if (filetype == "jpeg") {
      jpeg(filename = filename, width = .jpeg.size[1], 
           height = .jpeg.size[2], pointsize = 12)
    }
    else if (filetype == "pdf") {
      .pdf.size <- c(12,12)
      pdf(file = filename, width = .pdf.size[1], 
          height = .pdf.size[2]) #, pointsize = 9)
    }
  }
  
  if (separate.plots) {
    par(mfrow=c(1, 1), mar=c(5.1, 5.1, 3.1, 2.1))
  }
  else {
    par(mfrow=c(n, 1), mar=c(5.1, 5.1, 3.1, 2.1))
  }

  
  for (i in 1:n) {
    m <- methods[ind][i]
    plot(x[[m]], text=toupper(m), ...)
    if (separate.plots && i<n && ask && missing(filename)) {
      answer <- readline("Hit <Return> to see next plot:")
    }
  }
  if (!missing(filename)) {
    dev.off()
  }
}

plot.tdt <- function(x, text, ...) {
  p <- x$pval[,,drop=FALSE]
  n <- ncol(p)
  nmethod <- nrow(p)
  pval.methods <- rownames(p)

  gender.txt <- ""
  if (0) {
    if (!is.null(x$select.gender)) {
      gender.txt <- paste("- Trios with", c("male", "female")[x$select.gender], "children")
    }
    else {
      gender.txt <- "- All trios"
    }
  }
  
  main.txt <- paste("P-values from", text, gender.txt)
  
  ##  hist(p, main=main.txt)
  ##  points(p, rep(0,n))
  yrange <- range(0, -log(p, base=10))
  plot(1,1, ..., type="n", xlim=c(1, n), ylim=yrange, main=main.txt,xlab="Marker", ylab=expression(-log[10](p)), cex.lab=1.3)
  for (i in 1:nmethod) {
    points(-log(p[i,], base=10), pch=i, col=1)
    ##  p.exp <- seq(1/(n+1), n/(n+1), len=n)
    ##  lines(p.exp)
    ##  abline(a=0,b=0, col="red")
  }
  abline(a=-log(0.05, base=10),b=0, col="red", lty=2)
  text(n, -log(0.05, base=10)+0.025*diff(yrange), "P=0.05", col="blue", cex=1.2)
  abline(a=-log(0.01, base=10),b=0, col="red", lty=2)
  text(n, -log(0.01, base=10)+0.025*diff(yrange), "P=0.01", col="blue", cex=1.2)
  legend("bottomleft", pval.methods, pch=1:nmethod, cex=1.3, bty="n")
}


######################################################################################
#################### LOCAL FUNCTIONS                ##################################
######################################################################################

#### Transmission disequilibrium tests (tdt, hhrr, trimm)
## Input parameters
##  - trios: data containing allele or haplotype information.
##           Assume furthermore the columns: orig.lines, sex
##  - nsim.perm: Number of simulations used in permutation test
##  - select.gender: Do the analysis for a gender subset. Values: 1, 2, or NULL.
##                   1: Male, 2: Female, NULL: All
##  - method: character vector containing which of the methods tdt, hhrr, trimm to be used.
##  - marker.names: Names given to the markers. If NULL, markers are named 1,2,...
##  - xcrom: TRUE if xchromosome marker
##  - use.ambiguous: TRUE if ambiguous trios are to be kept
tdt <- function(trios, nsim.perm=1000, select.gender=NULL, method=c("tdt", "hhrr", "trimm"),
                xchrom=FALSE, marker.names=NULL, use.ambiguous=TRUE) {
  if (xchrom) {
    cols <- c("m1", "m2", "sex", "orig.lines", "freq", "marker") #orig.lines = line numbers in input haplin file
  }
  else {
    cols <- c("m1", "m2", "f1", "f2", "orig.lines", "freq", "marker")
  }

  ##  Select subset (gender, non-ambiguous, columns)
  trios <- SelectSubset(trios=trios, select.gender=select.gender, use.ambiguous=use.ambiguous,
                        col=cols)

  if (!xchrom) { # number of columns == 7
    ## Transform to 5 columns by appending the father alleles to the mother alleles
    ## and adding a 'par' (parent) variable
    n <- nrow(trios)
    m12 <- cbind(trios[,c("m1", "m2", "orig.lines", "freq", "marker")], par=rep("m", n))
    f12 <- cbind(trios[,c("f1", "f2", "orig.lines", "freq", "marker")], par=rep("f", n))
    colnames(f12) <- colnames(m12)
    trios <- rbind(m12, f12)
    colnames(trios)[1:2] <- c("h1", "h2")
    cols <- colnames(m12)
  }
  else { # xchrom. number of columns == 6
    ## Transform to 5 columns by adding a 'par' (parent) variable
    n <- nrow(trios)
    colnames(trios)[1:2] <- c("h1", "h2")
    trios <- cbind(trios, par=rep("m", n))
  }
  
  save.seed <- .Random.seed
  on.exit(set.seed(save.seed))
  set.seed(21) ## To assure that same results are obtained for given input data

  object <- tdt.basic(trios, nsim.perm=nsim.perm,
                      method=method, marker.names=marker.names)

  object
}

## Basic tdt test on 2+4 columns
## allele1 allele2 orig.lines m/f freq marker
tdt.basic <- function(trios, nsim.perm=1000,
                      method=c("tdt", "hhrr", "trimm"), marker.names=NULL) {
  if (missing(trios)) {
    stop("trios missing")
  }
##  nColumns <- ncol(trios)
  
##  if (nColumns != 6) {
##    stop("Number of columns in trios should be equal to 6")
##  }

  ## Reorder haplotypes (or alleles) such that they are numbered consecutively
  allele.cols <- c("h1", "h2")
  res <- ReorderAlleles(trios, allele.cols)
  trios <- res$trios
  orig.alleles <- res$orig.alleles

  ## Check that max(allele) = # alleles
  allele <- NULL
  for (i in 1:length(allele.cols)) {
    allele <- c(allele, trios[, allele.cols[i]])
  }
  nAllele <- length(table(allele))
  if (max(allele) != nAllele) {
    stop("Max allele count != No alleles")
  }

  ## Interrupt tdt test if number of alleles == 1
  if (nAllele == 1) {
    warning("Number of alleles == 1")
    obj <- list(nAllele=1, pvalues=NULL)
    return(obj)
  }

  markers <- as.integer(names(table(trios[,"marker"])))
  nmarker <- length(markers)
  if (max(markers) != nmarker) {
    warning("Max marker count != No of markers (markers: ", paste(markers, collapse=" "), ")\n")
  }

  ## Find the methods to be used (set doTdt, doHhrr, doTrimm)
  if ("tdt" %in% method) {
    doTdt <- TRUE
  }
  else {
    doTdt <- FALSE
  }
  if ("hhrr" %in% method) {
    doHhrr <- TRUE
  }
  else {
    doHhrr <- FALSE
  }
  if ("trimm" %in% method) {
    if (nAllele > 2) {
      warning("TRIMM not available when number of alleles > 2")
      doTrimm <- FALSE
    }
    else {
      doTrimm <- TRUE
    }
  }
  else {
    doTrimm <- FALSE
  }

  ## TDT and HHRR: Construct transmission table
  if (doTdt || doHhrr) {
    trans.table <- MakeTriosTable(trios, nAllele, markers)
    dimnames(trans.table) <- list(orig.alleles, orig.alleles, marker.names)
    n <- nrow(trios)
  }

  ## TRIMM: Construct difference vector (TDT may also be based on this when #alleles == 2)
  if (doTrimm || (doTdt && nAllele == 2)) {
    D <- DifferenceVector(trios)
    n.diffvec <- table(D[,"marker"])
    ##    n.denom.tdt <-  tapply(trios[,"h1"]!=trios[,"h2"], trios[,"marker"], sum)
  }
  
  method.names <- character()
  if (doHhrr) { ## haplotype-based haplotype relative risk
    hhrr.test <- apply(trans.table, 3, Hhrr.test)
    names(hhrr.test) <- marker.names
    df.hhrr <- nAllele - 1
    p.hhrr <- 1-pchisq(hhrr.test, df=df.hhrr)
    hhrr.sim <- matrix(NA, nrow=nmarker, ncol=nsim.perm)
    method.names <- c(method.names, "hhrr")
  }
  
  if (doTdt) { ## transmission disequilibrium test
    testObj <- matrix(unlist(apply(trans.table, 3, Tdt.test, return.df=TRUE)), nrow=2)
    tdt.test <- testObj[1,]
    names(tdt.test) <- marker.names
    ##    if (nAllele == 2) {
    ##      tdt.test.2allele <- Tdt.test.2allele(D[,1], D[,"marker"], n.denom.tdt)
    ##    }
    
    df.tdt <- testObj[2,]
    p.tdt <- 1-pchisq(tdt.test,df=df.tdt)
    tdt.sim <- matrix(NA, nrow=nmarker, ncol=nsim.perm)
    method.names <- c(method.names, "tdt")
  }

  if (doTrimm) { # TRiad Multi-Marker test
    trimm.test <- Trimm.test(D[,1], D[,"marker"], n.diffvec)
    names(trimm.test) <- marker.names
    df.trimm <- 1
    
    p.trimm <- 1-pchisq(trimm.test, df=df.trimm)
    trimm.sim <- matrix(NA, nrow=nmarker, ncol=nsim.perm)
    method.names <- c(method.names, "trimm")
  }


  ## ## PERMUTATATION TESTS ## ##

  if (doTrimm) {
    ## orig.lines = line numbers in input haplin file (= fam number)
    famlevel <- sort(unique(D[,"orig.lines"])) # family numbers records in D
    fam <- match(D[,"orig.lines"], famlevel) # family numbers reordered 1...#fam
    nfam <- length(famlevel) # number of families present in D
  }

  if (doTdt || doHhrr) {
    n <- nrow(trios)
    trios.sim <- data.frame(h1=rep(NA, n), h2=rep(NA, n), trios[,c("freq", "marker")])
    
    fampar <- paste(trios[,"orig.lines"], trios[,"par"],sep=":") # family:parent(m/f)
    famparlevel <- sort(unique(fampar))
    parent <- match(fampar, famparlevel) # reordering 1,2...#family:parent
    nparent <- length(famparlevel)
  }
  
  for (i in seq(len=nsim.perm)) {
    if (doTrimm) {
      swap <- sample(c(-1,1), size=nfam, replace=TRUE) # swap whole families (mother and father together)
      D.sim <- swap[fam] * D[,1]
      
      trimm.sim[,i] <- Trimm.test(D.sim, D[,"marker"], n.diffvec)
    }
    if (doTdt || doHhrr) {
      ## Swap parents (mother or father independantly of each other)
      ind.swap <- sample(c(FALSE, TRUE), size=nparent, replace=TRUE)
      ind.swap <- ind.swap[parent]
      
      trios.sim[ind.swap,1:2] <- trios[ind.swap,2:1]
      trios.sim[!ind.swap,1:2] <- trios[!ind.swap,1:2]
      table.sim <- MakeTriosTable(trios.sim, nAllele, nmarker)
      if (doTdt) {
        tdt.sim[,i] <- apply(table.sim, 3, Tdt.test)
      }
      if (doHhrr) {
        hhrr.sim[,i] <- apply(table.sim, 3, Hhrr.test)
      }
    }
  }
  
  res <- vector("list", length(method.names))
  index <- 1
  if (doHhrr) {
    if (nsim.perm > 0) {
      p.hhrr.perm <- (1/nsim.perm) * apply((hhrr.sim - hhrr.test >= 0), 1, sum)
      p.hhrr.perm <- pmax(1/nsim.perm, p.hhrr.perm)
      pval.hhrr <- rbind(p.hhrr, p.hhrr.perm)
      rownames(pval.hhrr) <- c("chisq", "permutation")
      
      hhrr.sim.max <- apply(hhrr.sim, 2, max)
      pval.hhrr.max <- mean(hhrr.sim.max >= max(hhrr.test))
      pval.hhrr.max <- pmax(1/nsim.perm, pval.hhrr.max)
    }
    else {
      pval.hhrr <- rbind(p.hhrr)
      rownames(pval.hhrr) <- "chisq"

      pval.hhrr.max <- NULL
    }
    
    colnames(pval.hhrr) <- marker.names

    object.hhrr <- list(transmissionMatrix=trans.table,
                        chisq=hhrr.test,
                        pval=pval.hhrr,
                        chisq.max=max(hhrr.test),
                        pval.max=pval.hhrr.max,
                        nAllele=nAllele,
                        df=df.hhrr)
    
    class(object.hhrr) <- "tdt"
    res[[index]] <- object.hhrr
    index <- index + 1
  }

  if (doTdt) {
    tdt.sim.max <- apply(tdt.sim, 2, max)
    if (nAllele == 2) {
      if (nsim.perm>0) {
        p.tdt.perm <- (1/nsim.perm) * apply((tdt.sim - tdt.test >= 0), 1, sum)
        p.tdt.perm <- pmax((1/nsim.perm), p.tdt.perm)
        
        pval.tdt.max <- mean(tdt.sim.max >= max(tdt.test))
        pval.tdt.max <- pmax((1/nsim.perm), pval.tdt.max)

        pval.tdt <- rbind(p.tdt, p.tdt.perm)
        rownames(pval.tdt) <- c("chisq", "permutation")
      }
      else {
        pval.tdt.max <- NULL

        pval.tdt <- rbind(p.tdt)
        rownames(pval.tdt) <- "chisq"
      }
      pval.max <- pval.tdt.max
    }
    else {
      if (nsim.perm > 0) {
        p.tdt.perm <- (1/nsim.perm) * apply((tdt.sim - tdt.test >= 0), 1, sum)
        p.tdt.perm <- pmax((1/nsim.perm), p.tdt.perm)
        
        pval.tdt.max <- mean(tdt.sim.max >= max(tdt.test))
        pval.tdt.max <- pmax((1/nsim.perm), pval.tdt.max)

        pval.tdt <- rbind(p.tdt, p.tdt.perm)
        rownames(pval.tdt) <- c("chisq", "permutation")
      }
      else {
        pval.tdt.max <- NULL
        
        pval.tdt <- rbind(p.tdt)
        rownames(pval.tdt) <- "chisq"
      }
      pval.max <- pval.tdt.max
      
    }

    colnames(pval.tdt) <- marker.names

    object.tdt <- list(transmissionMatrix=trans.table,
                       chisq=tdt.test,
                       pval=pval.tdt,
                       chisq.max=max(tdt.test),
                       pval.max=pval.max,
                       nAllele=nAllele,
                       df=df.tdt)    
    class(object.tdt) <- "tdt"
    res[[index]] <- object.tdt
    index <- index + 1
  }

  if (doTrimm) {
    if (nsim.perm > 0) {
      p.trimm.perm <- (1/nsim.perm) * apply((trimm.sim - trimm.test >= 0), 1, sum)
      p.trimm.perm <- pmax((1/nsim.perm), p.trimm.perm)
      
      pval.trimm <- rbind(p.trimm, p.trimm.perm)
      rownames(pval.trimm) <- c("chisq", "permutation")

      trimm.sim.max <- apply(trimm.sim, 2, max)
      pval.trimm.max <- mean(trimm.sim.max >= max(trimm.test))
      pval.trimm.max <- pmax((1/nsim.perm), pval.trimm.max)
    }
    else {
      pval.trimm <- rbind(p.trimm)
      rownames(pval.trimm) <- "chisq"

      pval.trimm.max <- NULL
    }
    colnames(pval.trimm) <- marker.names

    object.trimm <- list(chisq=trimm.test,
                         pval=pval.trimm,
                         chisq.max=max(trimm.test),
                         pval.max=pval.trimm.max,
                         nAllele=nAllele,
                         df=df.trimm)
    class(object.trimm) <- "tdt"
    res[[index]] <- object.trimm
  }
  
  names(res) <- method.names
  
  res
}

SelectSubset <- function(trios, use.ambiguous=TRUE, col=1:2, select.gender=NULL) {
  ## Exclude missing values (orig.lines that are repeated more than once)
  if (!use.ambiguous) {
    if ("marker" %in% colnames(trios)) {
      orig.lines <- paste(trios[,"marker"], trios[,"orig.lines"], sep="-")
    }
    else {
      orig.lines <- trios[,"orig.lines"]
    }
    orig.lines.table <- table(orig.lines)
    lines.notmissing <- names(orig.lines.table[orig.lines.table==1])
    ind <- orig.lines %in% lines.notmissing
    trios <- trios[ind,]
  }

  ## Select gender
  if (!is.null(select.gender)) {
    gender <- trios[,"sex"]
    ind.gender <- gender == select.gender
    trios <- trios[ind.gender,]
  }

  ## Extract columns of interest
  trios[,col]
}

ReorderAlleles <- function(trios, cols) {
  ## Reorder haplotypes (or alleles) such that they are numbered consecutively
  trios.cols <- trios[,cols]
  orig.alleles <- names(table(as.matrix(trios.cols)))
  for (i in 1:length(orig.alleles)) {
    trios.cols[trios.cols==orig.alleles[i]] <- i
  }
  trios[,cols] <- trios.cols

  list(trios=trios, orig.alleles=orig.alleles)
}


DifferenceVector.nofreq <- function(trios) {
  group <- paste(trios[,"orig.lines"], trios[,"marker"], sep=":")
  diffvec <- tapply(trios[,2], group, sum) - tapply(trios[,1], group, sum)
  diffvec <- diffvec[!is.na(diffvec)] # skip those orig.lines:marker combinations that are not in the data
  diffvec <- diffvec[diffvec!=0] # skip uninformative markers
  group.matrix <- matrix(unlist(strsplit(names(diffvec), ":")), ncol=2, byrow=TRUE)
  
  data.frame(diff=diffvec, orig.lines=as.integer(group.matrix[,1]), marker=group.matrix[,2])
}

DifferenceVector <- function(trios) {
  group <- paste(trios[,"orig.lines"], trios[,"marker"], sep=":")
  diffvec <- tapply((trios[,2]-trios[,1])*trios[,"freq"], group, sum)
  diffvec <- diffvec[!is.na(diffvec)] # skip those orig.lines:marker combinations that are not in the data
  diffvec <- diffvec[diffvec!=0] # skip uninformative markers
  group.matrix <- matrix(unlist(strsplit(names(diffvec), ":")), ncol=2, byrow=TRUE)
  
  data.frame(diff=diffvec, orig.lines=as.integer(group.matrix[,1]), marker=group.matrix[,2])
}


MakeTriosTable <- function(trios, nAllele, markers) {
  nmarker <- length(markers)
  trans.table <- array(NA, dim=c(nAllele, nAllele, nmarker))
  for (i in seq(len=nmarker)) {
    ind <- trios[,"marker"] == markers[i]
    trans.table[,,i] <- MakeTriosTable0(trios[ind,c("h1", "h2", "freq")], nAllele)
  }

  trans.table
}

MakeTriosTable.nofreq <- function(trios, nAllele, markers) {
  nmarker <- length(markers)
  trans.table <- array(NA, dim=c(nAllele, nAllele, nmarker))
  for (i in seq(len=nmarker)) {
    ind <- trios[,"marker"] == markers[i]
    trans.table[,,i] <- MakeTriosTable0(trios[ind,1:2], nAllele)
    diag(trans.table[,,i]) <- 0
  }

  trans.table
}

MakeTriosTable0 <- function(trios, nAllele) {
  h1h2 <- paste(trios[,1],trios[,2], sep=":")

  all <- expand.grid(1:nAllele, 1:nAllele)
  nm <- paste(as.character(all[,1]), as.character(all[,2]), sep=":")

  tab <- tapply(trios[, "freq"], h1h2, sum)
  
  h1h2.full <- rep(0, length(nm))
  names(h1h2.full) <- nm
  h1h2.full[names(tab)] <- tab
  
  trans.table <- matrix(h1h2.full, ncol=nAllele, byrow=TRUE)
  
  trans.table
}

MakeTriosTable0.nofreq <- function(trios, nAllele) {
  h1h2 <- table(paste(trios[,1],trios[,2], sep=":"))

  all <- expand.grid(1:nAllele, 1:nAllele)
  
  nm <- paste(as.character(all[,1]), as.character(all[,2]), sep=":")
  h1h2.full <- rep(0, length(nm))
  names(h1h2.full) <- nm
  h1h2.full[names(h1h2)] <- as.integer(h1h2)
  
  trans.table <- matrix(h1h2.full, ncol=nAllele, byrow=TRUE)
  
  trans.table
}

Hhrr.test <- function(trans.table) {
  row.sum <- apply(trans.table, 1, sum)
  col.sum <- apply(trans.table, 2, sum)
  numerator <- (row.sum-col.sum)^2
  denominator <- row.sum+col.sum
  ind <- denominator > 0
  if (any(!ind)) {
    warning("denominator == 0")
  }
  hhrr.test <- sum(numerator[ind]/denominator[ind])

  hhrr.test
}

Tdt.test <- function(trans.table, return.df=FALSE) {
  numerator <- (trans.table-t(trans.table))^2
  denominator <- trans.table+t(trans.table)
  ind <- (row(numerator) < col(numerator)) & (denominator > 0)
  tdt.test <- sum(numerator[ind]/denominator[ind])
  if (return.df) {
    df <- sum(ind)
    return(list(chisq=tdt.test, df=df))
  }
  tdt.test
}

Tdt.test.2allele <- function(D, marker, n) {
  tdt.test <- (tapply(D, marker, sum))^2/n
  
  as.vector(tdt.test)
}

Trimm.test <- function(D, marker, n) {
  D.mean <- tapply(D, marker, mean)
  D.mean.vec <- D.mean[marker]
  D.var <- tapply((D-D.mean.vec)^2, marker, sum)/(n*(n-1))
  
  trimm.test <- as.vector(D.mean^2/D.var)
  names(trimm.test) <- names(D.mean)
  
  trimm.test
}
