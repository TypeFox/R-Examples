#' Convert string to a numeric.
#'
#' This function removes spaces from a string and converts it to a numeric.
#' This is required when reading EDF data, since the data fields in the header
#' are of fixed size, and extra whitespace must hence be trimmed.
#'
#' @param s The string to be converted
#' @return A numeric.
#' @examples
#' s <- "123   "
#' edf:::edf.char.to.num(s)
#' @keywords internal
edf.char.to.num <- function(s) {
    as.numeric(gsub(" ", "", s, fixed = TRUE))
}


#' Remove trailing whitespace from a string.
#'
#' This function removes trailing whitespace from a string.
#' This is required when reading EDF data, since the data fields in the header
#' are of fixed size, and extra whitespace must hence be trimmed.
#'
#' @param s The string to be converted
#' @return The string gives as input without trailing whitespace
#' @examples
#' s <- "abc   "
#' edf:::trim.end(s)
#' @keywords internal
trim.end <- function(s) {
    gsub("[[:space:]]*$", "", s)
}


#' Create variable name from string.
#'
#' This function creates a sensible variable name by replacing
#' all non-alphanumeric characters with underscores. This is useful
#' when the string is to be used as a the name of a list element.
#'
#' @param s The string to be used as a variable name
#' @return The string with non-alphanumeric characters replaced by underscores
#' @examples
#' s <- "a-b"
#' edf:::create.variable.name(s)
#' @keywords internal
create.variable.name <- function(s) {
    gsub("[^[:alnum:]]", "_", s)
}


#' #' Parse the global header of a European Data Format (EDF and EDF+) file.
#'
#' This function parses the header of an EDF or EDF+ file.
#'
#' @param data The 256 bytes of raw data that contain the EDF/EDF+ header.
#' @return A named list containing the information in the EDF/EDF+ header.
#' @keywords internal
parse.edf.global.header <- function(data) {
    header <- vector(mode = "list")

    blocks <- list("data.format.version"  = list(n.bytes = 8,  type = "numeric"),
                   "patient.id"           = list(n.bytes = 80, type = "character"),
                   "recording.id"         = list(n.bytes = 80, type = "character"),
                   "startdate"            = list(n.bytes = 8,  type = "character"),
                   "starttime"            = list(n.bytes = 8,  type = "character"),
                   "bytes.in.header"      = list(n.bytes = 8,  type = "numeric"),
                   "reserved"             = list(n.bytes = 44, type = "character"),
                   "n.data.records"       = list(n.bytes = 8,  type = "numeric"),
                   "data.record.duration" = list(n.bytes = 8,  type = "numeric"),
                   "n.signals"            = list(n.bytes = 4,  type = "numeric"))

    bytes.start <- 1

    for (b in names(blocks)) {
        block <- blocks[[b]]
        tmp <- readBin(data[bytes.start:(bytes.start+block$n.bytes - 1)], character(), size = block$n.bytes, n = 1)

        if (block$type == "character")
            header[[b]] <- trim.end(tmp)
        if (block$type == "numeric")
            header[[b]] <- edf.char.to.num(tmp)
        bytes.start <- bytes.start + block$n.bytes
    }

    ## Add timestamp information
    ts <- paste(header$startdate, header$starttime, sep = "T")
    ts <- as.POSIXct(strptime(as.character(ts), "%d.%m.%yT%H.%M.%S"))
    
    header$timestamp.start <- ts
    header$timestamp.stop <- ts + (header$n.data.records * header$data.record.duration)

    header
}


#' Parse the signal header of a European Data Format (EDF and EDF+) file.
#'
#' This function parses the signal header of an EDF or EDF+ file.
#'
#' @param data The n*256 bytes of raw data that contain the EDF/EDF+ signal header,
#' @param header.global The global EDF/EDF+ header.
#' #' where n is the number of signals in the EDF/EDF+ file.
#' @return A named list containing the information in the EDF/EDF+ signal header.
#' @keywords internal
parse.edf.signal.header <- function(data, header.global) {
    n.signals <- header.global$n.signals
    header <- vector(mode = "list", length = n.signals)

    for (i in seq.int(n.signals)) {
        header[[i]] <- list()
    }


    blocks <- list("label"              = list(n.bytes =  16, type = "character"),
                   "transducer"         = list(n.bytes =  80, type = "character"),
                   "physical.dimension" = list(n.bytes =  8,  type = "character"),
                   "physical.minimum"   = list(n.bytes =  8,  type = "numeric"),
                   "physical.maximum"   = list(n.bytes =  8,  type = "numeric"),
                   "digital.minimum"    = list(n.bytes =  8,  type = "numeric"),
                   "digital.maximum"    = list(n.bytes =  8,  type = "numeric"),
                   "prefiltering"       = list(n.bytes =  80, type = "character"),
                   "n.samples"          = list(n.bytes =  8,  type = "numeric"),
                   "reserved"           = list(n.bytes =  32, type = "character"))

    bytes.start <- 1

    for (b in names(blocks)) {
        block <- blocks[[b]]
        for (i in seq.int(n.signals)) {
            tmp <- readBin(data[bytes.start:(bytes.start+block$n.bytes - 1)], character(), size = block$n.bytes, n = 1)
            if (block$type == "character")
                header[[i]][[b]] <- trim.end(tmp)
            if (block$type == "numeric")
                header[[i]][[b]] <- edf.char.to.num(tmp)
            bytes.start <- bytes.start + block$n.bytes
        }
    }

    for (i in seq.int(n.signals)) {
        header[[i]]$label        <- create.variable.name(header[[i]]$label)
        names(header)[i]         <- header[[i]]$label
        header[[i]]$samplingrate <- header[[i]]$n.samples / header.global$data.record.duration
        header[[i]]$gain         <- (header[[i]]$physical.maximum - header[[i]]$physical.minimum) / (header[[i]]$digital.maximum - header[[i]]$digital.minimum)
    }

    header
}


#' Parse one EDF+ Time-stamped Annotation List (TAL).
#'
#' This function parses one EDF+ Time-stamped Annotation List (TAL).
#'
#' @param data The bytes of raw data that contain the TAL.
#' @return A data frame containing the information in the TAL: onset, duration and annotation.
#' @keywords internal
parse.event <- function(data) {
    event <- readBin(as.raw(data), character(), size = 1, signed = FALSE, endian = "little")
    ## Split into two parts: (1) onset, duration and (2) event annotation
    p1 <- unlist(strsplit(event, split = "\024"))
    p2 <- unlist(strsplit(p1[1], split = "\025"))
    onset      <- as.numeric(gsub("+", "", p2[1]))
    duration   <- as.numeric(p2[2])
    if (is.na(duration))
        duration <- 0
    annotation <- p1[2]
    data.frame(record = NA, onset, duration, annotation, stringsAsFactors = FALSE)
}

#' Parse the annotation data in an EDF+ file.
#'
#' This function parses the annotations in an EDF+ file.
#'
#' @param data The data containining the annotations (as 8-bit unsigned integers).
#' @return A data frame containing the annotatations in the EDF+ file.
#' @keywords internal
parse.edf.annotations <- function(data) {
    ## Read the annotations
    annotations <- data.frame(record = numeric(), onset = numeric(), duration= numeric(), event = character(), stringsAsFactors = FALSE)

    ## --- process each TAL within the records
    n.records <- length(data)
    for (rc in seq.int(n.records)) {
        record <- data[[rc]]

        ## --- count the number of events
        last.event.stop <- rev(which(record == 20))[1] + 1
        event.start     <- which(record == 43)
        n.events        <- length(event.start)

        if (n.events == 1) {
            event.stop <- last.event.stop
        } else {
            event.stop <- c(event.start[2:length(event.start)] - 1, last.event.stop)
        }

        events          <- vector(mode = "list", length = n.events)

        ## --- separate the events
        for (i in seq.int(n.events)) {
            events[[i]] <- record[event.start[i]:event.stop[i]]
        }

        ## --- process the events
        events.df        <- do.call(rbind, lapply(events, parse.event))
        events.df$record <- rc
        annotations      <- rbind(annotations, events.df)
    }

    annotations
}
