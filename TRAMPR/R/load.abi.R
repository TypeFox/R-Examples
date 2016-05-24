## load.abi.R (part of the TRAMPR package)

## Functions to load and process an .abi file.

## The basic procedure is to call load.abi.template() to create a
## template, then edit the template, then call load.abi() to create
## the TRAMPsamples object.

## write.TRAMPsamples can then be used to output the final object.

## "template" files: These are created by load.abi.template() and used
## by load.abi().  These contain all the unique values of
## "sample.file.name" from an .abi file, and require entries for
## "enzyme",  indicating which enzyme was used and "sample.index",
## which is used to group different enzymes together as a single
## sample.  If "sample.index" is an integer field, then this is used
## as "sample.pk".  Otherwise an integer key is generated that matches
## sample.index.

## This loads an .abi file and creates a "template" file.
load.abi.create.template <- function(file, file.template) {
  if ( missing(file.template) )
    file.template <- load.abi.create.template.name(file)
  d <- read.abi(file)
  res <- data.frame(sample.file.name=sort(unique(d$sample.file.name)),
                    enzyme=NA, sample.index=NA)
  write.csv(res, file.template, row.names=FALSE, na="")
  cat(sprintf("Saved template file in %s\n",
              file_path_as_absolute(file.template)))
  invisible(res)
}

## This generates an 'info' file from a 'template' file.
load.abi.create.info <- function(file, file.template, file.info) {
  if ( missing(file.template) )
    file.template <- load.abi.create.template.name(file)
  if ( missing(file.info) )
    file.info <- load.abi.create.info.name(file)
  d <- read.csv.safe(file.template)
  res <- data.frame(sample.index=sort(unique(d$sample.index)),
                    species=NA)
  write.csv(res, file.info, row.names=FALSE, na="")
  cat(sprintf("Saved info file in %s\n",
              file_path_as_absolute(file.info)))
  invisible(res)
}

load.abi.create.info.name <- function(file)
  sprintf("%s_info.csv", file_path_sans_ext(file))

load.abi.create.template.name <- function(file)
  sprintf("%s_template.csv", file_path_sans_ext(file))

## This loads the .abi file, template file and optional information
## file and creates the TRAMPsamples object.
## TODO: Be more noisy about how NAs are handled?
load.abi <- function(file, file.template, file.info,
                     primer.translate, ...) {
  ## (1) ABI data
  d <- read.abi(file)
  must.contain.cols(d, c("sample.file.name", "dye", "size", "height"))
  d$primer <- NA
  for ( primer in names(primer.translate) )
    d$primer[d$dye %in% primer.translate[[primer]]] <- primer
  if ( any(is.na(d$primer)) )
    warning("Unknown dyes not in primer.translate: ",
            paste(dQuote(sort(na.omit(unique(d$dye[is.na(d$primer)])))),
                  collapse=", "))

  ## (2) Template data (with sample.info and enzyme)
  if ( missing(file.template) )
    file.template <- load.abi.create.template.name(file)
  d.template <- read.csv.safe(file.template)
  must.contain.cols(d.template,
                    c("sample.file.name", "sample.index", "enzyme"))

  good.sample.index <- is.integer(d.template$sample.index)
  if ( good.sample.index )
    d.template$sample.fk <- d.template$sample.index
  else {
    warning("Found non-integer index: creating new sample.pk\n",
            "Please see ?load.abi about referencing your data")
    d.template$sample.fk <-
      as.integer(factor(d.template$sample.index))
  }

  i <- match(d$sample.file.name, d.template$sample.file.name)
  cols <- setdiff(names(d.template), "sample.file.name")
  data <- cbind(d.template[i,cols], d)

  ## Remove cases where sample.index is NA, since these have no way of
  ## grouping together, and orphan data is not allowed (see
  ## ?TRAMPsamples).
  data <- data[!is.na(data$sample.index),]
  d.template <- d.template[!is.na(d.template$sample.index),]
  rownames(data) <- seq(length=nrow(data))

  ## Order the columns in a meaningful way:
  data.cols <- c("sample.fk", "sample.index", "sample.file.name",
                 "primer", "enzyme", "size", "height")
  data <- data[c(data.cols, setdiff(names(data), data.cols))]

  ## (TOCHECK) Remove cases where any of the data.cols are missing.
  data <- data[complete.cases(data[data.cols]),]

  ## (3) Create "info" data.frame
  if ( missing(file.info) )
    file.info <- load.abi.create.info.name(file)

  info <- unique(d.template[c("sample.fk", "sample.index")])
  info <- info[do.call(order, info),]
  colnames(info)[1] <- "sample.pk"
  
  if ( file.exists(file.info) ) {
    cat(sprintf("Found info file at %s\n", file.info))
    info.extra <- read.csv.safe(file.info)
    if ( any(duplicated(info.extra$sample.index)) )
      stop("Duplicated sample.index in file.info")

    missing <- setdiff(info$sample.index, info.extra$sample.index)
    if ( length(missing) > 0 ) {
      warning("sample.index values missing from file.info: ",
              paste(missing, collapse=", "),
              "\nThese have been removed!")
      info <- subset(info, info$sample.index %in% info.extra$sample.index)
      data <- subset(data, data$sample.fk %in% info$sample.pk)
    }

    unknown <- setdiff(info.extra$sample.index, info$sample.index)
    if ( length(unknown) > 0 ) {
      warning("Unknown sample.index values in file.info: ",
              paste(unknown, collapse=", "),
              "\nThese have been ignored!")
      info.extra <- subset(info.extra,
                           info.extra$sample.index %in%
                           info$sample.index)
    }

    ## Now, do this:
    i <- match(info$sample.index, info.extra$sample.index)
    cols <- setdiff(names(info.extra), "sample.index")
    info <- cbind(info, info.extra[i,cols,drop=FALSE])
  } else {
    cat("No info file found, will create default.\n")
  }

  rownames(info) <- seq(length=nrow(info))
  if ( good.sample.index )
    info <- info[setdiff(names(info), "sample.index")]

  TRAMPsamples(data, info, ...)
}

## Read in the .abi file format.  This is fairly straightforward,
## except that the .abi files apparently have a trailing tab at the
## end of each line.
read.abi <- function(file) {
  header <- read.table(file, sep="\t", nrows=1, as.is=TRUE)
  d <- read.table(file, sep="\t", skip=1, as.is=TRUE)
  n <- ncol(header)

  cols <- c("Dye/Sample Peak", "Sample File Name", "Size",  "Height")
  if ( !all(cols %in% header) )
    stop(sprintf("Required columns missing in abi file: %s",
                 paste(dQuote(cols[!(cols %in% header)]),
                       collapse=", ")))

  if ( ncol(d) != n )
    stop("Invalid file (header/body mismatch)")

  if ( is.na(header[[n]]) )
    if ( all(is.na(d[[n]])) ) {
      d <- d[-n]
      header <- header[-n]
    } else stop("Inconsistent trailing whitespace")
  names(d) <- gsub("[^0-9a-z]", ".", tolower(as.character(header[1,])))

  code <- strsplit(d$dye.sample.peak, ",")

  if ( !all(sapply(code, length) == 2) )
    stop("Invalid 'Dye/Sample Peak' codes")

  code <- as.data.frame(t(matrix(unlist(code), 2)))
  code[] <- lapply(code, as.character)
  names(code) <- c("dye", "sample.peak")

  d <- cbind(d, code)
  is.data <- regexpr("^[0-9]+$", d$sample.peak) == 1
  d <- d[is.data,]
  d$sample.peak <- as.integer(d$sample.peak)

  d[setdiff(names(d), "dye.sample.peak")]
}

## Peakscanner apparently produces files that are subtly different
## from abi genemapper files.  see ?peakscanner.to.genemapper for more
## information (in man/load.abi.Rd)
peakscanner.to.genemapper <- function(filename, output) {
  if ( missing(output) )
    output <- sub("^(.*)\\.[^.]*$", "\\1.txt", filename)
  if ( filename == output )
    stop(sprintf("Can't read and write from same file (%s)",
                 dQuote(filename)))
  d <- read.csv(filename, check.names=FALSE, as.is=TRUE)
  d[["Dye/Sample Peak"]] <- gsub(" ", "", d[["Dye/Sample Peak"]])
  write.table(d, output, sep="\t", row.names=FALSE)
}
