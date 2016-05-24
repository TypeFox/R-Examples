#' Read a European Data Format (EDF and EDF+) file.
#'
#' This function reads the data stored in an EDF or EDF+ file. This data
#' consists of, e.g., physiologic signals and possibly also annotations (EDF+ files only).
#'
#' @param filename The full path to the EDF/EDF+ file to be read.
#' @param read.annotations Boolean denoting whether or not annotations should be read,
#'        if they are present. Defaults is TRUE.
#' @param header.only Boolean denoting whether to only read the headers in the EDF file. Default is FALSE.
#' #' @return A list containing
#' \enumerate{
#' \item the global header (recording start date, duration etc).
#' \item the signal headers (transducer types, sampling rates etc).
#' \item the signal data
#' \item the events (annotations). This is always NA for EDF files. This is NA for EDF+
#'       files without annotations.
#' }
#' @references Kemp B., V\"{a}rri, A., Rosa, A.C., Nielsen, K.D. and Gade, J. (1992).
#'             A simple format for exchange of digitized polygraphic recordings.
#'             Electroencephalogr Clin Neurophysiol. 1992 May;82(5):391-3.
#'             \url{http://www.ncbi.nlm.nih.gov/pubmed/1374708}
#'
#'             Kemp, B. and Olivan, J. (2003). European data format 'plus' (EDF+),
#'             an EDF alike standard format for the exchange of physiological data.
#'             Clin Neurophysiol. 2003 Sep;114(9):1755-61.
#'             \url{http://www.ncbi.nlm.nih.gov/pubmed/12948806}
#'
#' @export
read.edf <- function(filename, read.annotations = TRUE, header.only = FALSE) {
    ## Read headers
    edf.file      <- file(filename, "rb")
    header.global <- parse.edf.global.header(readBin(edf.file, "raw", size = 1, n = 256, signed = FALSE, endian = "little"))
    header.signal <- parse.edf.signal.header(readBin(edf.file, "raw", size = 1, n = header.global$n.signals*256, signed = FALSE, endian = "little"), header.global = header.global)

    ## Only return the header
    if (header.only) {
        edf <- list("header.global" = header.global, "header.signal" = header.signal)

        close(edf.file)
        return(edf)
    }

    ## Preallocate signal vectors
    signal.sample.count <- 1 + numeric(length = header.global$n.signals)
    signal <- vector(mode = "list", length = header.global$n.signals)

    for (i in seq.int(header.global$n.signals)) {
        signal[[i]]       <- list()
        signal[[i]]$data  <- vector(mode = "numeric", length = header.signal[[i]]$n.samples * header.global$n.data.records)
        signal[[i]]$t     <- (0:(length(signal[[i]]$data)-1)) / header.signal[[i]]$samplingrate
        names(signal)[i]  <- header.signal[[i]]$label
    }

    ## Does the measurement have annotations?
    if ((header.global$reserved %in% c("EDF+C", "EDF+D")) & ("EDF_Annotations" %in% names(header.signal))) {
        annotations_id  <- which(names(header.signal) == "EDF_Annotations")
        data <- vector(mode = "list", header.global$n.data.records)
    } else {
        annotations_id <- FALSE
    }

    ## Read signal data
    samples.max.1   <- length(signal[[1]]$data)
    signal.seq      <- seq.int(header.global$n.signals)
    annotation.data <- vector(mode = "list", length = header.global$n.data.records)
    
    j  <- 1
    ac <- 1

    while (signal.sample.count[1] <= samples.max.1) {
        for (i in signal.seq) {
            if (i == annotations_id) {
                tmp <- readBin(edf.file, integer(), size = 1, n = 2*header.signal[[i]]$n.samples, signed = FALSE, endian = "little")
                annotation.data[[ac]] <- tmp
                signal[[i]]$data[signal.sample.count[i]:(signal.sample.count[i] + 2*header.signal[[i]]$n.samples - 1)] <- tmp
                ac <- ac + 1
            } else {
                signal[[i]]$data[signal.sample.count[i]:(signal.sample.count[i] + header.signal[[i]]$n.samples - 1)] <- readBin(edf.file, integer(), size = 2, n = header.signal[[i]]$n.samples, signed = TRUE, endian = "little")
            }
            signal.sample.count[i] <- signal.sample.count[i] + header.signal[[i]]$n.samples
        }
    }

    close(edf.file)

    ## Scale the signals properly
    for (i in names(signal)) {
        signal[[i]]$data  <- (signal[[i]]$data  - header.signal[[i]]$digital.minimum ) * header.signal[[i]]$gain + header.signal[[i]]$physical.minimum
    }

    
    ## Create an EDF structure that will be returned
    edf <- list("header.global" = header.global, "header.signal" = header.signal, "signal" = signal, "events" = NA)

    ## Read annotations
    if (read.annotations) {
        if (annotations_id) {
            ## Clean up the EDF structure
            edf$events                        <- parse.edf.annotations(annotation.data)
            edf$header.signal$EDF_Annotations <- NULL
            edf$signal$EDF_Annotations        <- NULL
        }
    }

    edf
}
