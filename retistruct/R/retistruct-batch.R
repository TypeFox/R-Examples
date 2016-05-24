##' List valid datasets underneath a directory. This reports all
##' directories that appear to be valid.
##'
##' @title List datasets underneath a directory
##' @param path Directory path to start searching from
##' @param verbose If \code{TRUE} report on progress
##' @return A vector of directories containing datasets
##' @author David Sterratt
##' @export
list.datasets <- function(path='.', verbose=FALSE) {
  ## We have to define a function to do directory listings as
  ## (a) list.files doesn't work recursively on all platforms and (b)
  ## in any case, it doesn't list directories when running recursively
  list.dirs <- function(path='.') {
    files <- list.files(path, recursive=FALSE, full.names=TRUE)
    fi <- file.info(files)
    dirs <- row.names(fi[fi$isdir,])
    for (d in dirs) {
      dirs <- c(dirs,list.dirs(d))
    }
    return(dirs)
  }

  ## Now get the directories
  dirs <- list.dirs(path)

  ## Go through directories determining if datasets are valid
  datasets <- c()
  for (d in dirs) {
    ## Determine if directory is a valid data directory. If it's a
    ## faulty one, we will let it pass to be picked up later on
    
    id.data.dir <- TRUE
    ## Case of faulty directory
    tryCatch({
      is.data.dir <- is.character(checkDatadir(d))
    }, error=function(e) {})

    if (!is.data.dir) {
      if (verbose) message(paste(d, "is not a data directory."))
      next
    }
    datasets <- c(datasets, d)
    if (verbose) message(paste(d, "is a data directory."))
  }
  return(datasets)
}

##' This function reconstructs a number of  datasets, using the R
##' \code{parallel} package to distribute the reconstruction of
##' multiple datasets across CPUs. If \code{datasets} is not specified
##' the function recurses through a directory tree starting at
##' \code{tldir}, determining whether the directory contains valid raw
##' data and markup, and performing the reconstruction if it does.
##'
##' @title Batch operation using the parallel package
##' @param tldir If datasets is not specified, the top level of the
##' directory tree through which to recurse in order to find datasets.
##' @param outputdir directory in which to dump a log file and images
##' @param datasets Vector of dataset directories to reconstruct
##' @param device string indicating what type of graphics output
##' required. Options are "pdf" and "png".
##' @param titrate Whether to "titrate" the reconstruction for
##' different values of \code{phi0}. See \code{titrate.reconstructedOutline}.
##' @param cpu.time.limit amount of CPU after which to terminate the
##' process
##' @param mc.cores The number of cores to use. Defaults to the total
##' number available.
##' @author David Sterratt
##' @export
retistruct.batch <- function(tldir='.', outputdir=tldir, datasets=NULL, 
                             device="pdf", titrate=FALSE,
                             cpu.time.limit=3600,
                             mc.cores=getOption("cores")) {
  ## Get datasets
  if (is.null(datasets)) {
    datasets <- list.datasets(tldir)
  }
  message(paste("About to reconstruct", length(datasets), "datasets."))

  ## Function to pass to mclapply
  call <- function(dataset) {
    logfile <- file.path(outputdir,
                         paste(retistruct.cli.basepath(dataset),
                               ".log", sep=""))
    message(paste("Reconstructing. Logging to", logfile))
    flog <- file(logfile, open="wt")
    sink(flog)
    sink(flog, type="message")
    return(retistruct.cli(dataset, cpu.time.limit, outputdir, device,
                          titrate=titrate))
  }

  ## Run the reconstructions
  ret <- parallel::mclapply(datasets, call, mc.preschedule=FALSE,
                  mc.cores=mc.cores)

  ## Extract data from the return structures
  ## Function to replace NULL with NA - needed for creating data frames
  n <- function(x) {
    return(ifelse(is.null(x), NA, x))
  }

  dat <- data.frame(cbind(dataset=datasets,
                          status=sapply(ret, function(x) {return(n(x$status))}),
                          mess=as.character(sapply(ret, function(x) {return(x$mess)})),
                          time=sapply(ret, function(x) {return(n(x$time))})))
  ## Defensive output, unless anything goes wrong at the next stage
  write.csv(dat, file.path(outputdir, "retistruct-batch-out.csv"))
  
  summ <- retistruct.batch.summary(tldir)
  write.csv(dat, file.path(outputdir, "retistruct-batch-summ.csv"))
  
  dat <- merge(summ, dat, by="dataset", all=TRUE)
  write.csv(dat, file.path(outputdir, "retistruct-batch.csv"))
  return(dat)
}

##' Recurse through a directory tree, determining whether the
##' directory contains valid derived data and extracting summary data
##' if it does.
##'
##' @title Extract summary data for a batch of reconstructions
##' @param tldir The top level directory of the tree through which to
##' recurse.
##' @param cache If \code{TRUE} use the cached statistics rather than
##' generate on the fly (which is slower).
##' @return Data frame containing summary data
##' @author David Sterratt
##' @export
retistruct.batch.summary <- function(tldir=".", cache=TRUE) {
  datasets <- list.datasets(tldir)
  logdat <- data.frame()

  ## Function to replace NULL with NA - needed for creating data frames
  n <- function(x) {
    return(ifelse(is.null(x), NA, x))
  }

  ## Go through datasets
  for (dataset in datasets) {
    message(paste("Reading", dataset))
    suppressMessages(r <- retistruct.read.recdata(list(dataset=dataset),
                                                  check=FALSE))
    if (!is.null(r)) {
      dat <- data.frame(dataset=dataset,
                        E=n(r$opt$value),
                        El=n(r$E.l),
                        nflip=n(r$nflip),
                        EOD=n(r$EOD),
                        sqrt.E=n(sqrt(r$E.l)),
                        mean.strain=n(r$mean.strain),
                        mean.logstrain=n(r$mean.logstrain),
                        OD.phi=n(r$Dss$OD[1,"phi"]),
                        OD.lambda=n(r$Dss$OD[1,"lambda"]),
                        mean.dtheta=n(r$titration$Dtheta.mean),
                        phi0d=n(r$phi0*180/pi),
                        phi0d.opt=n(r$titration$phi0d.opt),
                        L.rim=getFlatRimLength(r),
                        A.tot=r$A.tot)
      hullarea <- getDssHullarea(r)
      if (length(hullarea) > 0) {
        dat <- data.frame(dat, hullarea=hullarea)
      }
      message(paste("Getting KDE"))
      KDE <- getKDE(r)
      if (length(KDE) > 0) {
        ## Get out bandwidths by going through each component of the KDE
        KDEdat <- lapply(KDE, function(x) {x$h})
        names(KDEdat) <- paste("kde.h.", names(KDEdat), sep="")
        dat <- cbind(dat, KDEdat)
        ## Get out contour areas by going through each component of the KDE
        for (name in names(KDE)) {
          KDEdat <- as.list(KDE[[name]]$tot.contour.areas[,"contour.areas"])
          names(KDEdat) <- paste("kde.c", KDE[[name]]$tot.contour.areas[,"labels"], "." , name, sep="")
          dat <- cbind(dat, KDEdat)
        }
      }
      message(paste("Getting KR"))
      KR <- getKR(r)
      if (length(KR) > 0) {
        ## Get out bandwidths by going through each component of the KR
        KRdat <- lapply(KR, function(x) {x$h})
        names(KRdat) <- paste("kr.h.", names(KRdat), sep="")
        dat <- cbind(dat, KRdat)
        ## Get out contour areas by going through each component of the KR
        for (name in names(KR)) {
          KRdat <- as.list(KR[[name]]$tot.contour.areas[,"contour.areas"])
          names(KRdat) <- paste("kr.c", KR[[name]]$tot.contour.areas[,"labels"], "." , name, sep="")
          dat <- cbind(dat, KRdat)
        }
      }

      logdat <- merge(logdat, dat, all=TRUE)
    }
  }
  return(logdat)
}

##' Recurse through a directory tree, determining whether the
##' directory contains valid derived data and plotting graphs if it
##' does.
##'
##' @title Plot figures for a batch of reconstructions
##' @param tldir The top level directory of the tree through which to
##' recurse.
##' @param outputdir Directory in which to dump a log file and images
##' @param ... Parameters passed to plotting functions
##' @author David Sterratt
##' @export
retistruct.batch.figures <- function(tldir=".", outputdir=tldir, ...) {
  datasets <- list.datasets(tldir)
  for (dataset in datasets) {
    message(paste("Attempting to produce figures from", dataset))
    try(retistruct.cli.figure(dataset, outputdir, ...))
  }
}

##' @title Get titrations from a directory of reconstructions
##' @export
##' @param tldir The top level directory of the tree through which to
##' recurse. The files have to have been reconstructed with the
##' \code{titrate} option to \code{\link{retistruct.batch}}
retistruct.batch.get.titrations <- function(tldir=".") {
  datasets <- list.datasets(tldir)
  dat <- list()
  for (dataset in datasets) {
    r <- retistruct.read.recdata(list(dataset=dataset), check=FALSE)
    if (!is.null(r)) {
      ndat <- list(r$titration$dat)
      names(ndat) <- dataset
      dat <- c(dat, ndat)
    }
  }
  return(dat)
}

##' @title Plot titrations
##' @export
##' @param tdat Output of \code{\link{retistruct.batch.get.titrations}}
retistruct.batch.plot.titrations <- function(tdat) {
  dat <- array(NA, dim=c(length(tdat), nrow(tdat[[1]]), 2))
  summ <- NULL
  for (i in 1:length(tdat)) {
    d <- tdat[[i]]
    phi0d <- d["phi0", "phi0d"]
    if (!is.null(phi0d)) {
      ## Ignore retinae whose rim angle is 0 - this is the default and
      ## means that it probably hasn't been set properly
      if ((phi0d != 0) & (min(d[,"sqrt.E"]) < 0.2)) {
        dat[i,,2] <- d[,"sqrt.E"] - min(d[,"sqrt.E"])
        phi0d.opt <- d[which.min(d[,"sqrt.E"]), "phi0d"]
        dat[i,,1] <- d[,"phi0d"] - phi0d.opt

        sqrt.E.opt <- d[which.min(d[,"sqrt.E"]), "sqrt.E"]
        summ <- rbind(summ,
                      data.frame(dataset=names(tdat)[i],
                                 phi0d=d["phi0", "phi0d"] + 90,
                                 sqrt.E=d["phi0", "sqrt.E"],
                                 phi0d.opt=phi0d.opt + 90,
                                 sqrt.E.opt=min(d[, "sqrt.E"])))
      }
    }
  }

  ## Do the detailed plot of all the lines
  par(mar=c(2.4, 3.2, 0.7, 0.2))
  par(mgp=c(1.3, 0.3, 0), tcl=-0.3)
  par(mfcol=c(1, 1))
  par(cex=0.66)

  plot(NA, NA, xlim=range(dat[,,1], na.rm=TRUE), ylim=c(0, 0.2),
       xlab=expression(paste(phi[0] - hat(phi)[0])),
       ylab=expression(paste(italic(e)[L] - hat(italic(e))[L])))
  for (i in 1:length(tdat)) {
    lines(dat[i,,1], dat[i,,2], col="#00000040")
  }
  mtext("C", adj=-0.2, font=2, line=-0.7)
  dev.copy2pdf(file=file.path("retistruct-titration.pdf"), width=6.83/3, height=6.83/4)

  svg(filename=file.path("fig4-retistruct-titration.svg"), width=6.83/3, height=6.83/6)

  ## Plot of optimal rim lattitude versus the sub
  par(mar=c(2.2, 2.9, 0.5, 0.4))
  par(mgp=c(1.2, 0.3, 0), tcl=-0.3)
  par(mfcol=c(1, 2))
  par(cex=0.66)

  with(summ, plot(phi0d.opt ~ phi0d, col="white",
                  asp=1,
                  xlab=expression(italic(phi)[0]),
                  ylab=expression(hat(italic(phi))[0])))
  with(summ, boxplot(phi0d.opt ~ round(phi0d), at=unique(sort(round(phi0d))),
                      xaxt="n", add=TRUE))
  abline(0,1)
  abline( 10, 1, col="grey")
  abline( 20, 1, col="grey")
  abline(-10, 1, col="grey")
  abline(-20, 1, col="grey")
  panlabel("D")
  
  with(summ, plot(sqrt.E.opt ~ sqrt.E,
                  asp=1,
                  xlab=expression(italic(e)[L]),
                  ylab=expression(hat(italic(e))[L]),
                  pch=".", cex=3))
  abline(0,1)

  panlabel("E")
  dev.off()


  with(summ, print(summary(1 - sqrt.E.opt/sqrt.E)))
  
  return(list(dat=dat, summ=summ))
}

##' Recurse through a directory tree, determining whether the
##' directory contains valid derived data and converting r.rData files
##' to files in matlab format named r.mat
##'
##' @title Export data from reconstruction data files to matlab
##' @param tldir The top level of the directory tree through which to
##' recurse
##' @author David Sterratt
##' @export
retistruct.batch.export.matlab <- function(tldir=".") {
  datasets <- list.datasets(tldir)
  for (dataset in datasets) {
    r <- retistruct.read.recdata(list(dataset=dataset), check=FALSE)
    retistruct.export.matlab(r)
  }
}

##' Extract statistics from the retistruct-batch.csv summary file
##'
##' @title Extract statistics from the retistruct-batch.csv summary file
##' @param path The path to the retistruct-batch.csv
##' @return list of various statistics
##' @author David Sterratt
##' @export
retistruct.batch.analyse.summary <- function(path) {
  dat <- read.csv(file.path(path, "retistruct-batch.csv"))
  
  ## Detailed output codes
  etab <- sort(table(dat[,"mess"]), decreasing=TRUE)
  message("OUPUT CODES")
  message(rbind(format(etab, width=4), " ", paste(gsub("\n", "\n   ", gsub("\n$", "", names(etab))), "\n")))

  sqrt.E.fail <- 0.2
  
  ## Get number of failures due to lack of time
  N.outtime <- sum(na.omit(dat$status) == 1)

  ## Get successful reconstructions
  sdat <- subset(dat, dat$status==0)

  message("\nFAILURES")
  cpufail <- dat[grepl("CPU", dat[,"mess"]),]
  message(paste("\n", nrow(cpufail), "of", nrow(dat), "did not finish"))
  print(cpufail[,c("dataset")])
  
  ## Get number of failures due to phi0d not being set
  nophis <- subset(sdat, sdat$phi0d == 0 | is.na(sdat$phi0d))
  N.nophi <- nrow(nophis)
  message(paste("\n", nrow(nophis), "of", nrow(sdat), "retinae do not have the rim angle set:"))
  print(nophis[,c("dataset", "phi0d")])
  sdat <- subset(sdat, sdat$phi0d != 0 & !is.na(sdat$phi0d))
  
  ## Get number of failures due to sqrt.E being too large
  failures <- subset(sdat, sqrt.E >= sqrt.E.fail)
  N.fail <- nrow(failures)
  failures <- failures[order(failures[,"sqrt.E"], decreasing=TRUE),]
  message(paste(nrow(failures), "of", nrow(sdat), "retinae have e_L (sqrt.E) greater than", sqrt.E.fail, ":"))
  print(failures[,c("dataset", "sqrt.E")])
  sdat <- subset(sdat, sqrt.E < sqrt.E.fail)
  
  message("\nSTATISTICS")
  
  message("sqrt.E")
  sqrt.E <- summary(sdat[,"sqrt.E"])
  print(sqrt.E)
  message(paste("SD of sqrt.E is", sd(sdat[,"sqrt.E"]), "; sqrt.E mean + 2SDs is", mean(sdat[,"sqrt.E"]) + 2*sd(sdat[,"sqrt.E"])))
  message("mean.strain")
  mean.strain <- summary(sdat[,"mean.strain"])
  print(mean.strain)
  message("mean.logstrain")
  mean.logstrain <- summary(sdat[,"mean.logstrain"])
  print(mean.logstrain)
  message("time")
  time <- summary(sdat[,"time"])
  print(time)
  message("nflip")
  nflip <- summary(sdat[,"nflip"])
  print(nflip)
  message("%with.flips")
  with.flips <- mean(sdat[,"nflip"] > 0) * 100
  
  message("\nOUTLIERS")
  ## outliers <- subset(sdat, sqrt.E > (mean(sqrt.E) + 2*sd(sqrt.E)))
  outliers <- subset(sdat, sqrt.E >  0.1)
  outliers <- outliers[order(outliers[,"sqrt.E"], decreasing=TRUE),]
  message(paste(nrow(outliers), "of", nrow(sdat), "retinae have e_L (sqrt.E) greater than 0.1:"))
  print(outliers[,c("dataset", "sqrt.E")])

  dev.new(width=6.83/2, height=6.83/4)
  ## Plot of various things
  ## Figure 4 in PLoS paper
  par(mar=c(2.4, 2.6, 0.7, 0.2))
  par(mgp=c(1.4, 0.3, 0), tcl=-0.3)
  par(mfcol=c(1, 2))
  par(cex=0.66)

  ## Fig 4A: Histogram of goodness measure over all retinae
  hist(sdat[,"sqrt.E"], breaks=seq(0, max(sdat[,"sqrt.E"]), len=100),
       xlab=expression(italic(e)[L]), main="")
  panlabel("A")
  
  ## Fig 4B: Boxplot of age for retinae of different ages

  ## Find datasets containing ("Pxx")
  fage <- grepl("^.*(P\\d+).*$",  sdat$dataset)
  sdat$age <- sub("^.*P(\\d+).*$", "\\1", sdat$dataset)
  sdat$age <- sub(".*adult.*", "adult", sdat$age)
  sdat$age[grepl(".{6}", sdat$age)] <- NA
  sdat$age <- ordered(sdat$age, c(unique(sort(as.numeric(sdat$age))), "adult"))
  ## print(factor(sdat$age))
  ## print(sort(as.numeric(factor(sdat$age))))
  ## levels(sdat$age) <- sub("(\\d+)", "P\\1", levels(sdat$age))
  levels(sdat$age) <- sub("adult", "A", levels(sdat$age))
  ##  print(factor(sdat$age))
  with(sdat, boxplot(sqrt.E ~ age,
                     xaxt="n",
                     xlab="Postnatal day",
                     ylab=expression(italic(e)[L])))
  panlabel("B")
  axis(1, labels=NA, at=seq(1, len=length(levels(sdat$age))))
  mtext(levels(sdat$age), 1, at=seq(1, len=length(levels(sdat$age))), line=0.3, cex=0.66)
  dev.print(svg, file=file.path(path, "fig3-retistruct-deformation.svg"), width=6.83/2, height=6.83/4)

  ## Fig 4A-C: Locations of optic discs
  dev.new(width=6.83/2, height=6.83/6)
  par(mfcol=c(1, 3))
  par(mar=c(0.7,0.7,0.7,0.7))
  
  retistruct.batch.plot.ods(subset(sdat, sdat$age=="A"),
                            phi0d=subset(sdat, sdat$age=="A")[1, "phi0d"])
  panlabel("A")

  summ <- retistruct.batch.plot.ods(subset(sdat, sdat$age=="A"),
                            phi0d=-60)
  panlabel("B")

  par(mar=c(2.4, 2.6, 0.7, 0.5))
  par(mgp=c(1.3, 0.3, 0), tcl=-0.3)
  
  summlm <- lm(OD.res ~ sqrt.E, summ)
  with(summ, plot(sqrt.E, OD.res,
                  xlab=expression(italic(e)[L]),
                  ylab=expression(italic(epsilon)[OD]),
                  pch=20, col="blue"))
  abline(summlm)
  panlabel("C")
  dev.print(svg, file=file.path(path, "fig4-retistruct-ods.svg"), width=6.83/2, height=6.83/6)
  print(summary(summlm))
  
  ## with(sdat, table(sqrt.E ~ age))

  ## Print Fig. 4 to EPS file
  ## dev.print(postscript, file=file.path(path, "fig4-retistruct-goodness.eps"),
  ##           width=6.83, height=6.83/3,
  ##           onefile=FALSE, horizontal=TRUE)
  
  ## Plot of kernel density features
  par(mar=c(2.4, 2.3, 0.7, 0.2))
  par(mgp=c(1.4, 0.3, 0), tcl=-0.3)
  par(mfcol=c(1, 2))

  hist(na.omit(sdat[,"kde.h.red"]),
       breaks=seq(0, max(na.omit(sdat[,"kde.h.red"])), len=100),
       xlab=expression(italic(h)[red]), main="")
  mtext("A", adj=-0.15, font=2, line=-0.7)

  ## Plot of age versus goodness
  ## Find datasets containing ("Pxx")

  fage <- grepl("^.*(P\\d+).*$",  sdat$dataset)
  sdat$age <- sub("^.*P(\\d+).*$", "\\1", sdat$dataset)
  sdat$age <- sub(".*adult.*", "adult", sdat$age)
  sdat$age[grepl(".{6}", sdat$age)] <- NA
  sdat$age <- ordered(sdat$age, c(unique(sort(as.numeric(sdat$age))), "adult"))
  ## print(factor(sdat$age))
  ## print(sort(as.numeric(factor(sdat$age))))
  ## levels(sdat$age) <- sub("(\\d+)", "P\\1", levels(sdat$age))
  levels(sdat$age) <- sub("adult", "A", levels(sdat$age))
  boxplot(sdat$kde.h.red ~ sdat$age,
          xaxt="n",
          xlab="Postnatal day",
          ylab=expression(italic(h)[red]))
  mtext("B", adj=-0.15, font=2, line=-0.7)
  axis(1, labels=NA, at=seq(1, len=length(levels(sdat$age))))
  mtext(levels(sdat$age), 1, at=seq(1, len=length(levels(sdat$age))), line=0.3, cex=0.66)
  
  ## with(sdat, table(sqrt.E ~ age))
  dev.copy2pdf(file=file.path(path, "retistruct-bandwidth.pdf"), width=4.5, height=3)

  ## More KDE analysis -- first ignore outliers
  sdat <- cbind(sdat, genotype="W")
  levels(sdat$genotype) <- c("W", "B")
  sdat[grep("GMB", sdat$dataset), "genotype"] <- "B"
  kdat <- subset(sdat, sdat$kde.h.red < 4)
  
  ## Plot of kernel density features
  par(mar=c(2.4, 2.3, 0.7, 0.2))
  par(mgp=c(1.4, 0.3, 0), tcl=-0.3)
  par(mfcol=c(1, 1))
  
  with(kdat, boxplot(kde.h.red ~ genotype+age, col=c("white", "red"),
                     xlab="Postnatal day",
                     ylab=expression(italic(h)[red])))
  dev.copy2pdf(file=file.path(path, "retistruct-bandwidth-genotype.pdf"), width=12, height=8)


  
  return(invisible(list(N=nrow(sdat),
                        N.outtime=N.outtime,
                        N.fail=N.fail,
                        N.nophi=N.nophi,
                        sqrt.E=sqrt.E, mean.strain=mean.strain,
                        mean.logstrain=mean.logstrain,
                        time=time, nflip=nflip,
                        with.flips=with.flips,
                        outliers=outliers,
                        sdat=sdat)))
}

##' Extract statistics from a directory containing
##' reconstruction directories. 
##'
##' @title Extract statistics from a directory containing
##' reconstruction directories. 
##' @param path Directory containing recontstruction directories
##' @return Data frame containg various statistics 
##' @author David Sterratt
##' @export
retistruct.batch.analyse.summaries <- function(path) {
  files <- list.files(path, recursive=FALSE, full.names=TRUE)
  fi <- file.info(files)
  dirs <- row.names(fi[fi$isdir,])
  out <- data.frame()
  for (d in dirs) {
    file <- file.path(d, "retistruct-batch.csv")
    if (file.exists(file)) {
      print(file)
      summ <- try(retistruct.batch.analyse.summary(d))
      try(print(summ$sqrt.E["Median"]))
      try(out <- rbind(out, data.frame(file=file,
                                       N=summ$N,
                                       N.outtime=summ$N.outtime,
                                       N.fail=summ$N.fail,
                                       sqrt.E.Median=summ$sqrt.E["Median"],
                                       sqrt.E.Mean=summ$sqrt.E["Mean"],
                                       nflip.Median=summ$nflip["Median"],
                                       nflip.Mean=summ$nflip["Mean"],
                                       nflip.Max=summ$nflip["Max."],
                                       pc.flipped=summ$with.flips,
                                       time.Mean=summ$time["Mean"])))
    }
  }
  return(out)
}

##' Polar plot of ODs of a group of retinae.
##'
##' @title Superposed plot of ODs on polar axes
##' @param summ Summary object returned by
##' \code{\link{retistruct.batch.summary}}
##' @param phi0d The rim angle for the plot
##' @param ... Other parameters, passed to projection
##' @return A pseudo retina, in which the optic disks are treated as
##' datapoints 
##' @author David Sterratt
##' @export
retistruct.batch.plot.ods <- function(summ, phi0d, ...) {
  ## Make a dummy retina
  o <- list()
  class(o) <- "reconstructedOutline"
  o$phi0 <- phi0d*pi/180
  r <- RetinalDataset(o)
  r <- ReconstructedDataset(r)
  r <- RetinalReconstructedDataset(r)
  r <- RetinalReconstructedOutline(r)
  r$side <- "Right"
  summ <- subset(summ, summ$age=="A")
  r$Dss$OD <- na.omit(summ[,c("OD.phi","OD.lambda")])
  colnames(r$Dss$OD) <- c("phi", "lambda")
  r$cols["OD"] <- "blue"
  
  km <- karcher.mean.sphere(r$Dss$OD, na.rm=TRUE, var=TRUE)
  message(nrow(summ), " points")
  message("Mean: Lat ", format(km$mean["phi"]*180/pi, digits=3),
          " Long ", format(km$mean["lambda"]*180/pi, digits=3),
          " ; SD: ", format(sqrt(km$var)*180/pi, digits=3))
  message("Mean location is ", 180/pi*central.angle(km$mean["phi"],
                                      km$mean["lambda"],
                                      -pi/2,
                                      0),  " away from geometric centre")
  
  summ$OD.res <- 180/pi*central.angle(km$mean["phi"],
                                      km$mean["lambda"],
                                      r$Dss$OD[,"phi"],
                                      r$Dss$OD[,"lambda"])

  projection(r, datapoint.contours=FALSE, philim=c(-90, phi0d), ...)
  ## dev.new()
  ## with(summ, plot(sqrt.E, OD.res))
  ## abline(summlm)
  return(summ)
}
