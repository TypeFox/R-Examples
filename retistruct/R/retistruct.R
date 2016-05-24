## Global variables
##' @title Version of reconstruction file data format
##' @export
recfile.version <- 5      # Version of reconstruction file data format

## Report function, with similar arguments to print
retistruct.report <- function(message, title="",...) {
  cat(paste(message, "\n", sep=""))
}

##' @title Check the whether  directory contains valid data 
##' @param dir Diectory to check.
##' @return  \code{TRUE} if \code{dir} contains valid data;
##' \code{FALSE} otherwise.
##' @author David Sterratt
##' @export
checkDatadir <- function(dir=NULL) {
  if (idt.checkDatadir(dir))   { return("idt") }
  if (csv.checkDatadir(dir))   { return("csv") }
  if (ijroi.checkDatadir(dir)) { return("ijroi") }
  return(FALSE)
}

##' Read a retinal dataset in one of three formats; for information on
##' formats see see \code{\link{idt.read.dataset}},
##' \code{\link{csv.read.dataset}} and
##' \code{\link{ijroi.read.dataset}}. The format is autodetected from
##' the files in the directory.
##' 
##' @title Read a retinal dataset
##' @param dataset Path to directory containing as SYS and MAP file
##' @param ... Parameters passed to the format-specific functions. 
##' @return An object that of classes \code{\link{RetinalDataset}} and
##' \code{\link{RetinalDataset}}. There may be extra fields too,
##' depending on the format.
##' @author David Sterratt
##' @export
retistruct.read.dataset <- function(dataset, ...) {
  ## Check to see if dataset is valid
  type <- checkDatadir(dataset)
  
  if (type=="idt")   { return(idt.read.dataset(dataset, ...))}
  if (type=="csv")   { return(csv.read.dataset(dataset, ...))}
  if (type=="ijroi") { return(ijroi.read.dataset(dataset, ...))}

  stop("No valid dataset format detected.")
}

##' Test the oputline object \code{o} for the prescense of
##' potential optice disc. This is done by checking that the list of
##' landmark lines \code{Ss} exists.
##'
##' @title Test for a potential optic disc
##' @param o Outline object
##' @return \code{TRUE} if an optic disc may be present; \code{FALSE} otherwise
##' @author David Sterratt
##' @export
retistruct.potential.od <- function(o) {
  if (inherits(o, "dataset")) {
    return(with(o, exists("Ss")) && is.list(o$Ss) && (length(o$Ss) > 0))
  }
  return(FALSE)
}

##' Read the markup data contained in the files \file{markup.csv},
##' \file{P.csv} and \file{T.csv} in the directory \file{dataset},
##' which is specified in the reconstruction object \code{r}. 
##'
##' The tear information is contained in the files \file{P.csv} and
##' \file{T.csv}. The first file contains the locations of outline
##' points that the tears were marked up on. The second file contains
##' the indicies of the apicies and backward and forward verticies of
##' each tear. It is necessary to have the file of points just in case
##' the algorithm that determines \code{P} in
##' \code{\link{retistruct.read.dataset}} has changed since the markup
##' of the tears.
##'
##' The remaining information is contained  containted in the file
##' \file{markup.csv}.
##'
##' If \code{DVflip} is specified, the locations of points \code{P}
##' flipped in the \eqn{y}-direction. This operation also requires the
##' swapping of \code{gf}  and \code{gb} and \code{VF} and \code{VB}.
##' @title Read the markup data
##' @param a Dataset object, containing \code{dataset} path
##' @param error Function to run on error, by default \code{stop()}
##' @return o \code{RetinalDataset} object
##' \item{V0}{Indicies in \code{P} of apicies of tears}
##' \item{VB}{Indicies in \code{P} of backward verticies of tears}
##' \item{VF}{Indicies in \code{P} of backward verticies of tears}
##' \item{iN}{Index in \code{P} of nasal point, or \code{NA} if not marked}
##' \item{iD}{Index in \code{P} of dorsal point, or \code{NA} if not marked}
##' \item{iOD}{Index in \code{Ss} of optic disc }
##' \item{phi0}{Angle of rim in degrees}
##' \item{DVflip}{Boolean variable indicating if DV axis has been flipped}
##' @author David Sterratt
##' @export
retistruct.read.markup <- function(a, error=stop) {
  ## Return index in P of closest point to x
  closest <- function(P, x) {
    if (any(is.na(x))) {
      return(NA)
    }
    return(which.min(vecnorm(t(t(P) - x))))
  }

  ## Function to map old tears (M.old) and old points (P.old) onto new
  ## points (P). It returns a new matrix of indicies (M).
  convert.markup <- function(M.old, P.old, P) {
    M <- sapply(M.old, function(i) {
      ifelse(is.numeric(i), closest(P, P.old[i,]), i)
    })
    return(M)
  }

  ## Read in the old P data
  Pfile <- file.path(a$dataset, "P.csv")
  if (file.exists(Pfile)) {
    P.old <- as.matrix(read.csv(Pfile))
  } else {
    P.old <- a$P
  }
  
  ## Read in markup file
  markupfile <- file.path(a$dataset, "markup.csv")
  if (file.exists(markupfile)) {
    ## read.csv() returns a data frame
    M.df <- read.csv(markupfile)
    ## Get strings out before converting to vector
    if ("side" %in% names(M.df)) {
      a$side <- as.character(M.df[1,"side"])
      ## Make sure that side is a valid value
      if (!(a$side %in% c("Left", "Right"))) {
        a$side <- "Right"
        warning("The side was neither Left nor Right; setting to Right.")
      }
    }
    ## Convert to vector
    M <- sapply(M.df, function(x) x)
    if (!is.na(M["iD"])) {
      M["iD"] <- convert.markup(M["iD"], P.old, a$P)
      a <- setFixedPoint(a, M["iD"], "Dorsal")
    }
    if (!is.na(M["iN"])) {
      M["iN"] <- convert.markup(M["iN"], P.old, a$P)
      a <- setFixedPoint(a, M["iN"], "Nasal")
    }
    a$phi0 <- M["phi0"]*pi/180
    if ("iOD" %in% names(M)) {
      a <- nameLandmark(a, M["iOD"], "OD")
    }
    if ("DVflip" %in% names(M)) {
      a$DVflip <- M["DVflip"]
    }
  } else {
    error("Markup file M.csv doesn't exist.")
  }
  
  ## Read in tearfile
  tearfile <- file.path(a$dataset, "T.csv")
  if (file.exists(tearfile)) {
    T.old <- read.csv(tearfile)
    cn <- colnames(T.old)
    T <- matrix(convert.markup(as.matrix(T.old), P.old, a$P), ncol=3)
    colnames(T) <- cn
    for (i in 1:nrow(T)) {
      a <- addTear(a, T[i,])
    }
  } else {
    error("Tear file T.csv doesn't exist.")
  }
  return(a)
}

##' Check that markup such as tears and the nasal or dorsal points are present.
##' 
##' @title Retistruct check markup
##' @param o Outline object
##' @return If all markup is present, return \code{TRUE}. Otherwise
##' return \code{FALSE}.
##' @author David Sterratt
##' @export
retistruct.check.markup <- function(o) {
  return(!is.null(names(o$i0)))
}

##' Given an outline object with a \code{dataset} field,  read the
##' reconstruction data from the file \file{\var{dataset}/r.Rdata}.
##'
##' @title Read the reconstruction data from file
##' @param o Outline object containing \code{dataset} field
##' @param check If \code{TRUE} check that the base information in the
##' reconstruction object is the same as the base data in source
##' files.
##' @return If the reconstruction data exists, return a reconstruction
##' object, else return the outline object \code{o}.
##' @author David Sterratt
##' @export
retistruct.read.recdata <- function(o, check=TRUE) {
  recfile <- file.path(o$dataset, "r.Rdata")
  if (file.exists(recfile)) {
    load(recfile)                       # This puts r in the environment
    ## If the algorithm in the codebase is newer than in the recdata
    ## file, reject the recfile data
    if (is.null(r$version) || (r$version != recfile.version)) {
      unlink(recfile)
      warning("The algorithm has changed significantly since this retina was last reconstructed, so the cached reconstruction data has been deleted.")
      return(NULL)
    }
    ## If the base data doesn't match the recfile data, reject the
    ## recfile data
    if (check) {
      if (!isTRUE(all.equal(o, r))) {
        unlink(recfile)
        warning("The base data has changed since this retina was last reconstructed, so the cached reconstruction data has been deleted.")
        return(NULL)
      } 
    }
    ## Otherwise regenerate data derived from dataset, such as KDE and
    ## KR; this was not stored by retistruct.recdata.save()
    ## FIXME: This should be deleted when recfile.version is next incremented
    if (is.null(r$KDE) | is.null(r$KR)) {
      r <- ReconstructedDataset(r)
      r <- RetinalReconstructedDataset(r)
    }

    ## Make sure the dataset information isn't overwritten
    r$dataset <- o$dataset
    return(r)
  }
  return(NULL)
}

##' @title Reconstruct a retina
##' @param o \code{\link{AnnotatedOutline}} object
##' @param report Function to report progress
##' @param plot.3d If \code{TRUE} show progress in a 3D plot 
##' @param dev.flat The ID of the device to which to plot the flat
##' representation
##' @param dev.polar The ID of the device to which to plot the polar
##' representation
##' @param ... Parameters to be passed to \code{\link{ReconstructedOutline}}
##' @return Object of classes
##' \code{\link{RetinalReconstructedOutline}} and
##' \code{\link{RetinalReconstructedDataset}} that contains all the
##' reconstruction information
##' @author David Sterratt
##' @export
retistruct.reconstruct <- function(o, report=retistruct.report,
                                   plot.3d=FALSE, dev.flat=NA, dev.polar=NA,
                                   ...) {
  ## Check that markup is there
  if (!retistruct.check.markup(o)) {
    stop("Neither dorsal nor nasal pole specified")
  }

  ## Check tears are valid
  ct <- checkTears(o)
  if (length(ct)) {
    stop(paste("Invalid tears", toString(ct), "marked up. Fix using \"Move Point\"."))
  }

  ## Set up fixed point 
  o$lambda0 <- 0
  ## Case of dorsal point specified...
  if (names(o$i0)[1]=="Dorsal") {
    o$lambda0 <- 90 * pi/180
  }
  if (names(o$i0)[1]=="Nasal") {
    if (o$side=="Right") {
      o$lambda0 <- 0
    } else {
      o$lambda0 <- pi
    }
  }

  ## Now do folding itself
  r <- NULL
  r <- ReconstructedOutline(o,
                            report=report,
                            plot.3d=plot.3d, dev.flat=dev.flat,
                            dev.polar=dev.polar,
                            ...)
  if (!is.null(r)) {
    r <- ReconstructedDataset(r, report=report)
    if (!is.na(getLandmarkID(r, "OD"))) {
      SssMean <- getSssMean(r)
      r$EOD <- 90 + SssMean[["OD"]][1,"phi"] * 180/pi
    }
    report(paste("Mapping optimised. Deformation eL:", format(sqrt(r$E.l), 5),
                 ";", r$nflip, "flipped triangles. OD displacement:",
                 format(r$EOD, 2),
                 "degrees."))
        
    r <- RetinalReconstructedOutline(r, report=report)
    r <- RetinalReconstructedDataset(r, report=report)
    report("")
  }
  return(r)
}

##' Save the makrup in the \code{\link{RetinalDataset}} \code{a} to a
##' file called \code{markup.csv} in the directory \code{a$dataset}.
##'
##' @title Save markup
##' @param a \code{\link{RetinalDataset}} object
##' @author David Sterratt
##' @export
retistruct.save.markup <- function(a) {
  if (inherits(a, "retinalDataset")) {
    with(a, {
      ## Save the tear information and the outline
      write.csv(cbind(V0, VB, VF), file.path(dataset, "T.csv"),
                row.names=FALSE)
      write.csv(P, file.path(dataset, "P.csv"), row.names=FALSE)
      
      ## Save the dorsal and nasal locations and phi0 to markup.csv
      iD <- ifelse(names(i0) == "Dorsal", i0, NA)
      iN <- ifelse(names(i0) == "Nasal" , i0, NA)
      iOD <- which(names(Ss)=="OD")
      if (length(iOD)==0)
        iOD <- NA
      markup <- data.frame(iD=iD, iN=iN, phi0=phi0*180/pi, iOD=iOD, DVflip=DVflip, side=side)     
      write.csv(markup, file.path(dataset, "markup.csv"), row.names=FALSE)
    })
  }
}


##' Save the reconstruction data in an object \code{r}  that inherits
##' \code{\link{ReconstructedDataset}} and
##' \code{\link{ReconstructedOutline}}  to a file called
##' \code{r.Rdata} in the directory \code{r$dataset}.
##'
##' @title Save reconstruction data
##' @param r \code{\link{ReconstructedDataset}} object
##' @author David Sterratt
##' @export
retistruct.save.recdata <- function(r) {
  if (!is.null(r$dataset)) {
    ## Save the derived data
    r$version <- recfile.version        # Datafile version
    if (!is.null(r)) {
      save(r, file=file.path(r$dataset, "r.Rdata"))
    }
  }
}

##' Save as a MATLAB object certain fields  of an object \code{r}
##' that inherits \code{\link{ReconstructedDataset}} and
##' \code{\link{ReconstructedOutline}}  to a file called \code{r.mat}
##' in the directory \code{r$dataset}.
##'
##' @title Save reconstruction data in MATLAB format
##' @param r \code{\link{ReconstructedDataset}} object
##' @author David Sterratt
##' @export
retistruct.export.matlab <- function(r) {
  ## writeMat() doesn't seem to cope with hierarchical structures, so
  ## we have to unlist the KDE and KR objects using this function
  unlist.kernel.estimate <- function(KDE) {
    if (length(KDE) > 0) {
      for (i in 1:length(KDE)) {
        KDEi <- KDE[[i]]
        KDEi$flevels <- NULL
        KDEi$labels <- NULL
        KDEi$tot.contour.areas <- NULL
        KDEi <- unlist(KDEi, recursive=FALSE)
        KDE[[i]] <- list(flevels=KDE[[i]]$flevels, labels=KDE[[i]]$labels,
                         tot.contour.areas=KDE[[i]]$tot.contour.areas)
        KDE[[i]] <- c(KDE[[i]], KDEi)
        names(KDE[[i]]$contour.areas) <- c()
      }
      ## print(names(KDE))
      ## print(names(unlist(KDE, recursive=FALSE)))
      KDE <- unlist(KDE, recursive=FALSE)
      names(KDE) <- gsub('\\.', '_', names(KDE))
    }
    return(KDE)
  }
  
  if (!is.null(r$dataset)) {
    if (!is.null(r)) {
      f <- file.path(r$dataset, "r.mat")
      message(paste("Saving", f))
      KDE <- unlist.kernel.estimate(getKDE(r))
      KR <-  unlist.kernel.estimate(getKR(r))
      R.matlab::writeMat(f,
                         phi0=r$phi0*180/pi,
                         Dss=getDss(r),
                         DssMean=getDssMean(r),
                         DssHullarea=na.omit(getDssHullarea(r)),
                         Sss=name.list(getSss(r)),
                         Tss=name.list(getTss(r)),
                         KDE=KDE,
                         KR=KR,
                         side=as.character(r$side), DVflip=r$DVflip,
                         Dsw=lapply(r$Dsc, function(x) {sphere.cart.to.sphere.wedge(x, r$phi0 + pi/2, r$R)}),
                         Dsdw=lapply(r$Dsc, function(x) {sphere.cart.to.sphere.dualwedge(x, r$phi0 + pi/2, r$R)}))
    }
  }
}

