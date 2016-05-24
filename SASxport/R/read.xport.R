##
## Code originally from Frank Harrell's 'Hmisc' library:
##   http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/Hmisc
## Copied with permission on 2007-08-04
##

#' @importFrom Hmisc label label<- label<-.default label<-.data.frame
#' @importFrom utils download.file
#' @export

read.xport <- function(file,
                       force.integer=TRUE,
                       formats=NULL,
                       name.chars=NULL,
                       names.tolower=FALSE,
                       keep=NULL,
                       drop=NULL,
                       as.is=0.95, # Prevent factor conversion if 95% or more unique
                       verbose=FALSE,
                       as.list=FALSE,
                       include.formats=FALSE
                       )
  {
    sasdateform <-
      toupper(c("date","mmddyy","yymmdd","ddmmyy","yyq","monyy",
                "julian","qtr","weekdate","weekdatx","weekday","month"))
    sastimeform     <- toupper(c("hhmm","hour","mmss","time"))
    sasdatetimeform <- toupper(c("datetime","tod"))

    if(verbose)
      {
        oldOptions <- options("DEBUG")
        options(DEBUG=TRUE)
        on.exit(options(oldOptions))
      }

    if(length(grep('http://', file))>0 || length(grep('ftp://', file))>0 )
      {
        scat("Downloading file...")
        tf <- tempfile()
        download.file(file, tf, mode='wb', quiet=TRUE)
        file <- tf
      }

    scat("Checking if the specified file has the appropriate header")
    xport.file.header <- "HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!000000000000000000000000000000  "
    file.header <- readBin( file, what=character(0), n=1, size=nchar(xport.file.header) )
    file.header <- substr(file.header, start=1, stop=nchar(xport.file.header) )
    if( !identical(xport.file.header, file.header) )
      stop("The specified file does not start with a SAS xport file header!")

    scat("Extracting data file information...")
    dsinfo <- lookup.xport.inner(file)

    dsLabels <- sapply(dsinfo, label)
    dsTypes  <- sapply(dsinfo, SAStype)

    if(length(keep))
      whichds <- toupper(keep)
    else
      whichds <- setdiff(names(dsinfo), c(toupper(drop),'_CONTENTS_','_contents_'))

    scat("Reading the data file...")
    ds <- read.xport.inner(file, stringsAsFactors=FALSE)
    if(any(duplicated(names(dsinfo))))  # only true if file contains has more than one data set
       {
         warning("Duplicate data set names in file.  Data set names have been made unique.")
         names(dsinfo) <- make.unique(names(dsinfo))
         names(ds) <- make.unique(names(ds))
       }


    if( (length(keep)>0 || length(drop)>0) )
      ds <- ds[whichds]

    scat("Processing contents...")
    ## PROC FORMAT CNTLOUT= dataset present?
    fds <- NULL
    if(!length(formats)) {
      fds <- sapply(dsinfo, function(x)
                    all(c('FMTNAME','START','END','MIN','MAX','FUZZ')
                        %in% x$name))
      fds <- names(fds)[fds]
      if(length(fds) > 1) {
        warning('transport file contains more than one PROC FORMAT CNTLOUT= dataset; using only the first')
        fds <- fds[1]
      }
    }

    finfo <- NULL
    if(length(formats) || length(fds)) {
      if(length(formats))
        finfo <- process.formats(formats)
      else
        finfo <- process.formats(ds[[fds]])
    }

    ## Number of non-format datasets
    nods <- length(whichds)
    nds  <- nods - (length(formats) == 0 && length(finfo) > 0)
    which.regular <- setdiff(whichds, fds)

    ## Handle lowercase name conversions
    if(names.tolower)
      names.tolower <- tolower
    else
      names.tolower <- function(x) x

    dsn <- names.tolower(which.regular)

    res <- vector('list', nds)
    names(res) <- gsub('_','.',dsn)


    possiblyConvertChar <- (is.logical(as.is) && !as.is) ||
                           (is.numeric(as.is) && as.is > 0)
    j <- 0
    for(k in which.regular) {
      j   <- j + 1
      scat('Processing SAS dataset', k)
      w   <-
        if(nods==1)
          ds
        else ds[[k]]

      scat('.')

      label(w, self=TRUE)   <- dsLabels[k]
      names(label(w, self=TRUE)) <- NULL
      SAStype(w) <- dsTypes[k]
      names(SAStype(w)) <- NULL

      nam      <- names.tolower(makeNames(names(w), allow=name.chars))
      names(w) <- nam
      dinfo    <- dsinfo[[k]]

      fmt <- dinfo$format
      formats  <- fstr( fmt, dinfo$flength, dinfo$fdigits)

      ifmt <- dinfo$iformat
      iformats <- fstr( ifmt, dinfo$iflength, dinfo$ifdigits)

      lab      <- dinfo$label

      ndinfo   <- names.tolower(makeNames(dinfo$name, allow=name.chars))
      names(lab) <- names(fmt) <- names(formats) <- names(iformats) <- ndinfo
      if(length(w)>0)
        for(i in 1:length(w)) {
          changed <- FALSE
          x  <- w[[i]]
          fi <- fmt[nam[i]];
          names(fi) <- NULL
          if(fi != '' && length(finfo) && (fi %in% names(finfo))) {
            f <- finfo[[fi]]
            if(length(f)) {  ## may be NULL because had a range in format
              x <- factor(x, f$value, f$label)
              attr(x, 'SASformat') <- fi
              changed <- TRUE
            }
          }

          if(is.numeric(x)) {
            if(fi %in% sasdateform) {
              x <- importConvertDateTime(x, 'date', 'sas')
              changed <- TRUE
            } else if(fi %in% sastimeform) {
              x <- importConvertDateTime(x, 'time', 'sas')
              changed <- TRUE
            } else if(fi %in% sasdatetimeform) {
              x <- importConvertDateTime(x, 'datetime', 'sas')
              changed <- TRUE
            } else if(force.integer) {
              if(all(is.na(x))) {
                storage.mode(x) <- 'integer'
                changed <- TRUE
              } else if(max(abs(x),na.rm=TRUE) <= (2^31-1) &&
                        all(floor(x) == x, na.rm=TRUE)) {
                storage.mode(x) <- 'integer'
                changed <- TRUE
            }
            }
          } else if(possiblyConvertChar && is.character(x)) {
            if((is.logical(as.is) && !as.is) ||
               (is.numeric(as.is) && length(unique(x)) < as.is*length(x))) {
              x <- factor(x, exclude='')
              changed <- TRUE
            }
          }

          lz <- lab[nam[i]]
          if(!is.null(lz) && length(lz)>0 && !is.na(lz) && lz != '') {
            names(lz) <- NULL
            label(x)  <- lz
            changed   <- TRUE
          }

          if(nam[i] %in% names(formats)  && formats[nam[i]] > "" )
            {
              SASformat(x) <- formats[[nam[i]]]
              changed <- TRUE
            }

          if(nam[i] %in% names(iformats) && iformats[nam[i]] > "" )
            {
              SASformat(x) <- formats[[nam[i]]]
            changed <- TRUE
            }

          if(changed)
            w[[i]] <- x
        }

      scat('.')

      res[[j]] <- w
    }

    scat("Done")

    if( include.formats )
      {
        nds <- nds+1
        if( length(fds)>0 )
          res$"FORMATS" <- ds[[fds]]
        else
          res$FORMATS <- empty.format.table()
      }

    if(nds > 1 || as.list)
      res
    else
      if(class(w)=="list")
        w[[1]]
    else
      w
  }
