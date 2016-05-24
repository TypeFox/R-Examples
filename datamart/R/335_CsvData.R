#' CsvData -- tabular data and its metadata
#' 
#' This class provides the infrastructure to
#' read and write tabular data, including
#' metadata on the dataset as a whole and
#' metadata on the columns of the dataset.
#'
#' Metadata on the dataset that are not interpreted are
#' title, description, license, publisher, keywords.
#' Metadata on the dataset that is interpreted are
#' update.scheme (always, never, and everything that 
#' can be passed to seq.POSIXt as by parameter),
#' update.method (lastmod, fullreplace),
#' encoding (an element of iconvlist()), type (csv, csv2),
#' modified (automatically set, YYYY-mm-dd HH:MM:SS UTC format).
#'
#' Metadata on a dataset column that are currently interpreted are name,
#' type (logical, integer, numeric, complex, character, raw, factor, Date, POSIXct), 
#' format (for Date and POSIXct, format string for strptime). 
#'
#' @seealso \code{\link{csvdata}}, \code{\link{read.csvdata}}
#'
#' @examples
#' getSlots("CsvData")
#' 
#' @exportClass CsvData
#' @name CsvData-class
#' @rdname CsvData-class
setClass(
    Class="CsvData", 
    representation=representation(
        resource="character",
        name="character",
        location="Location",
        dset.meta="list", 
        cols.meta="data.frame",
        update.fct="function"
    ),
    contains=c("Xdata", "Target")
)

#' Read CSV from local directory or Internet
#'
#' It is lazy, i.e. it does not read the csv in until requested. It does read in the metadata, though.
#' If according to the dataset meta data an update is due, an update process is performed.
#'
#' @param resource        the name of the resource. Required.
#' @param name            the physical name of the resource. Defaults to resource.
#' @param location        either a Location object, or a character pointing to a local directory or an url. See details.
#' @param update.fct      function for updating the data. the update method and update interval is specified as meta data.
#'                        Default is a function that returns an empty data.frame.
#' @param verbose         Print diagnostic messages, default is TRUE.
#' @param clss            class to construct. Defaults to CsvData.
#'
#' @export
read.csvdata <- function(
    resource, 
    location=NULL, 
    name=resource,
    update.fct=function(csv) return(data.frame()),
    verbose=TRUE,
    clss="CsvData"
) {
    if(is.null(location)) {
        location <- normalizePath(dirname(resource))
        resource <- basename(resource)
    }
    if(is.character(location)) {
        if(isTRUE(file.info(location)$isdir)) 
            location <- dirloc(location)
        else if(RCurl::url.exists(file.path(location, paste(name, ".csv", sep=""))))
            location <- webloc(location)
    }
    dcfreader <- function(p) {
        con <- file(p)
        res <- read.dcf(con)
        close(con)
        return(res)
    }
    
    # dataset metadata
    if(verbose) cat("reading dataset metadata...")
    dset.meta <- try(query(location, paste(name, ".dsm", sep=""), extract.fct=dcfreader), silent=TRUE)
    if(inherits(dset.meta, "try-error")) stop("no dataset metadata '", paste(name, ".dsm", sep=""), "' found at location '", as.character(location), "'")
    dset.meta <- as.list(dset.meta[1,])
    names(dset.meta) <- tolower(names(dset.meta))
    
    invalid.meta <- setdiff(names(dset.meta), c(
        "title", "description", "license", "publisher", "keywords",
        "update.scheme", "update.method", "encoding", "type", "modified"
    ))
    if(length(invalid.meta)>0) stop("invalid dataset metadata: unknown key(s): '", paste(invalid.meta, collapse="', '"), "'.", sep="")
    
    tmp <- dset.meta[["na.strings"]]
    if(!is.null(tmp)) {
        tmp <- strsplit(tmp, ",[\\s]*", perl=TRUE)[[1]]
        tmp <- gsub("\"", "", tmp)
        dset.meta[["na.strings"]] <- tmp
    }

    tmp <- dset.meta[["keywords"]]
    if(!is.null(tmp)) {
        tmp <- strsplit(dset.meta[["keywords"]], ",[\\s]*", perl=TRUE)[[1]]
        tmp <- gsub("\"", "", tmp)
        dset.meta[["keywords"]] <- tmp
    }
    
    tmp <- dset.meta[["update.scheme"]]
    if(!is.null(tmp)) {
        tmp <- strsplit(tmp, " ", fixed = TRUE)[[1L]]
        if(length(tmp) > 2L || length(tmp) < 1L) 
            stop("invalid 'update.scheme' metadata: '", dset.meta[["update.scheme"]], "'", sep="")
        
        valid <- pmatch(tmp[length(tmp)], c("secs", "mins", "hours", "days", "weeks", "months", "years", "DSTdays", "quarters", "always", "never"))
        if(is.na(valid)) 
            stop("invalid 'update.scheme' metadata: '", dset.meta[["update.scheme"]], "'", sep="")
        
        if((tmp %in% c("always", "never")) && length(tmp!=1)) 
            stop("invalid 'update.scheme' metadata: '", dset.meta[["update.scheme"]], "'", sep="") 
    }
    else dset.meta[["update.scheme"]] <- "never"

    tmp <- dset.meta[["update.method"]]
    if(is.null(tmp)) tmp <- "none"
    if(tmp=="none" && dset.meta[["update.scheme"]] != "never") stop("inconsistent metadata: update.method='none' while update.scheme='", dset.meta[["update.scheme"]], "'.", sep="")
    if(!(tmp %in% c("lastmod", "fullreplace", "none"))) 
        stop("invalid 'update.method' metadata: '", tmp, "'. Must be one of 'none', 'lastmod' and 'fullreplace'.", sep="")
    
    tmp <- dset.meta[["modified"]]
    if(!is.null(tmp)) {
        tmp <- strptime(tmp, "%Y-%m-%d %H:%M:%S")
        if(!is.na(tmp)) dset.meta[["modified"]] <- tmp
            else stop("invalid 'modified' metadata, not in '%Y-%m-%d %H:%M:%S' format.")
    }
    if(is.null(tmp) && inherits(location, "DirectoryLocation")) {
        warning("no 'modified' metadata found, trying to use file attribute.")
        file.mtime <- function(p) file.info(p)[["mtime"]]
        tmp <- try(query(location, paste(name, ".csv", sep=""), extract.fct=file.mtime), silent=TRUE)
        if(!inherits(tmp, "try-error") && !is.na(tmp)) dset.meta[["modified"]] <- tmp
    }
    if(is.null(tmp) && !(dset.meta[["update.scheme"]] %in% c("always", "never")))
        stop("no 'modified' metadata found, but is needed for 'update.scheme=",dset.meta[["update.scheme"]], "'", sep="")
        
    # TODO: test if encoding really works.
    tmp <- dset.meta[["encoding"]]
    if(!is.null(tmp)) {
        if(tmp!="latin1" && !(tmp %in% iconvlist())) stop("invalid/unknown 'encoding' metadata: '", tmp, "'", sep="")
    }
    
    #column metadata
    if(verbose) cat("OK\nreading column metadata...")
    cols.meta <- try(query(location, paste(name, ".dcm", sep=""), extract.fct=dcfreader), silent=TRUE)
    if(!inherits(cols.meta, "try-error")) {
        colnames(cols.meta) <- tolower(colnames(cols.meta))
        cols.meta <- as.data.frame(cols.meta, stringsAsFactors=FALSE)
    } else 
        stop("no column metadata '", paste(name, ".dcm", sep=""), "' found at location '", as.character(location), "'")

    col.classes <- unique(cols.meta[["type"]])
    if(!is.null(col.classes)) {
        unknown.types <- setdiff(
            col.classes,
            c("logical", "integer", "numeric", "complex", "character", "raw", "factor", "NULL", "Date", "POSIXct")
        )
        if(length(unknown.types)>0) stop("invalid column types defined: '", paste(unknown.types, collapse=", "), "'")
    }
    
    # instantiate object and update data if update.scheme says so
    res <- new(clss, resource=resource, name=name, location=location, dset.meta=dset.meta, cols.meta=cols.meta)
    upd.sch <- dset.meta[["update.scheme"]]
    if(upd.sch != "never" && (upd.sch=="always" || seq(dset.meta[["modified"]], length.out=2, by=upd.sch)[2]<Sys.time())) {
        if(verbose) cat("OK\nupdating dataset...")
        res <- revise.csvdata(res)
    }
    if(verbose) cat("OK\n")
    return(res)
}

#' Create CsvData object from data.frame
#'
#' This function takes a data.frame, its metadata and a location. The function writes the data
#' to the location and returns a CsvData object. The data.frame itself is not stored in this object.
#'
#' If metadata is missing, it is set with sensible default values. The default values for the dataset metadata are
#' title=<resource name>, modified=<current timestamp>, encoding='UTF-8', type='csv2', all other entries missing.
#' The default values for column metadata are name=<column name of the data.frame>, type=<class(column)>, 
#' format='%Y-%m-%d' for Date and '%Y-%m-%d %H:%M:%S' for POSIXt classes (otherwise missing).
#'
#' No updating takes place here.
#'
#' @param resource        the name of the resource. Required.
#' @param dat             data.frame to convert to CsvData object. Required.
#' @param name            the physical name of the resource. Defaults to resource.
#' @param location        either a Location object, or a character pointing to a local directory. See details.
#' @param update.fct      function for updating the data. the update method and update interval is specified as meta data.
#'                        Default is a function that returns an empty data.frame.
#' @param dset.meta       A list of metadata on the dataset.
#' @param cols.meta       A data.frame of metadata on the columns, rows are column names, columns are type, name, format
#' @param verbose         print diagnostic messages (default=TRUE)
#' @param clss            class to construct. Defaults to CsvData.
#'
#' @export
# TODO redirect calls to as.csvdata here
csvdata <- function(
    resource, dat, 
    location, name=resource, 
    dset.meta=NULL, cols.meta=NULL, 
    update.fct=function(csv) return(data.frame()),
    verbose=TRUE,
    clss="CsvData"
    ) {
    arg.lst <- list()
    arg.lst[[paste(name, ".csv", sep="")]] <- as.data.frame(dat)
    loc <- do.call(memloc, arg.lst)
    if(is.null(dset.meta)) dset.meta <- list(title=resource, type="csv2", encoding="UTF-8", modified=Sys.time())
    if(is.null(cols.meta)) {
        cols.meta <- data.frame(name=colnames(dat), type=sapply(dat, FUN=class), format=NA, title=NA)
    }
    tmp <- new(clss, resource=resource, name=name, location=loc, dset.meta=dset.meta, cols.meta=cols.meta)
    put(tmp, location)
    
    read.csvdata(resource=resource, name=name, location=location, update.fct=update.fct)
}

# internal function/private method to update CsvData object.
revise.csvdata <- function(csv, ...) {
    method <- csv@dset.meta[["update.method"]]

    if(method=="lastmod") {
        dat <- query(csv, resource=csv@resource, for.update=TRUE)
        newdat <- csv@update.fct(csv)
        miss_col <- setdiff(colnames(dat), colnames(newdat))
        for (cn in miss_col) newdat[, cn] <- NA
        dat <- unique(rbind(dat, newdat[,colnames(dat)]))
    }
    
    if(method=="fullreplace") dat <- csv@update.fct(csv)
    
    res <- csvdata(
        resource=csv@resource, dat=dat, 
        dset.meta=csv@dset.meta, cols.meta=csv@cols.meta, 
        location=csv@location, name=csv@name, 
        update.fct=csv@update.fct,
        clss=class(csv)
    )
    return(res)
}

#' @param for.update     (CsvData) update before loading (default FALSE)
#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,CsvData,character-method
setMethod(
    f="query",
    signature=c(self="CsvData", resource="character"),
    definition=function(self, resource, verbose=getOption("verbose"), for.update=FALSE, ...) {
        if(resource==self@resource) {
            if(verbose) cat("loading CSV...")
           
            used <- intersect(
                names(self@dset.meta), 
                c("sep", "dec", "encoding", "na.strings", "skip", "blank.lines.skip")
            )
            
            arg.lst <- self@dset.meta[used]
            date.cols <- NULL
            time.cols <- NULL
            if(nrow(self@cols.meta)>0) {
                col.classes <- self@cols.meta[, "type"] 
                date.cols <- which(col.classes == "Date")
                time.cols <- which(col.classes == "POSIXt")
                col.classes[c(date.cols, time.cols)] <- "character"
                arg.lst[["colClasses"]] <- col.classes
            }
            reader <- function(p) {
                the.args <- arg.lst
                the.args[["file"]] <- p
                do.call(
                    paste("read.", tolower(self@dset.meta[["type"]]), sep=""), 
                    the.args
                )
            }
            res <- query(self@location, paste(self@name, ".csv", sep=""), extract.fct=reader)
            for (i in date.cols) {
                fmt <- self@cols.meta[i, "format"]
                if(is.na(fmt)) {
                    warning("no date format specified for column '", rownames(self@cols.meta)[[i]], "', using '%Y-%m-%d'")
                    fmt <- "%Y-%m-%d"
                }
                res[[i]] <- as.Date(strptime(res[[i]], format=fmt))
            }
            for (i in time.cols) { 
                fmt <- self@cols.meta[i, "format"]
                if(is.na(fmt)) {
                    warning("no date format specified for column '", rownames(self@cols.meta)[[i]], "', using '%Y-%m-%d %H:%M:%S'")
                    fmt <- "%Y-%m-%d %H:%M:%S"
                }
                res[[i]] <- strptime(res[[i]], format=fmt)
            }
            if(verbose) cat("OK.\n")
            return(res)
        }
        
        if(verbose) cat("trying inherited method..\n")
        callNextMethod(self=self, resource=resource, verbose=verbose, ...)
    }
)

#' @rdname queries-methods
#' @name queries
#' @export
#' @docType methods
#' @aliases queries queries,CsvData-method
setMethod(
    f="queries",
    signature="CsvData",
    definition=function(self) self@resource
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,CsvData-method
setMethod(
    f="meta",
    signature="CsvData",
    definition=function(self, ...) { 
        res <- self@cols.meta
        row.names(res) <- res[, "name"]
        res[, "name"] <- NULL
        dsm <- self@dset.meta
        for (nm in names(dsm)) attr(res, nm) <- dsm[[nm]]
        return(res)
    }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,CsvData,DirectoryLocation-method
setMethod(
    f="put",
    signature=c(target="CsvData", where="DirectoryLocation"),
    definition=function(target, where, verbose=TRUE, ...) {
        tostr <- function(s) paste(s, collapse=", ")
        csv <- query(target, target@resource)
        
        # dataset metadata
        if(verbose) cat("writing dataset metadata..")
        f <- file.path(as.character(where), paste(target@name, ".dsm", sep=""))
        dsm <- list()
        for (nm in names(target@dset.meta)) dsm[[nm]] <- target@dset.meta[[nm]]
        dsm <- lapply(target@dset.meta, FUN=tostr)
        dsm[["modified"]] <- strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S UTC", tz="UTC")
        dsm <- do.call(data.frame, dsm)
        write.dcf(x=dsm, file=f)
        
        # column metadata
        if(verbose) cat("OK\nwriting column metadata..")
        f <- file.path(as.character(where), paste(target@name, ".dcm", sep=""))
        dcm <- target@cols.meta
        type <- sapply(csv, FUN=class)
        idx <- which(dcm[, "type"] != type)
        if(length(idx)>0) {
            err <- paste("col ", dcm[["name"]][idx], ": expected ", type[idx], ", got ", dcm[,"type"][idx], sep="")
            stop("invalid type(s) found: ", paste(err, collapse="\n"))
        }
        for (i in which(dcm[, "type"]  == "Date")) {
            fmt <- dcm[i, "format"]
            if(is.na(fmt)) {
                fmt <- "%Y-%m-%d"
                dcm[i, "format"] <- fmt
            }
            csv[,i] <- strftime(csv[,i], format=fmt)
        }
        for (i in which(dcm[, "type"] == "POSIXct")) {
            fmt <- dcm[i, "format"]
            if(is.na(fmt)) {
                fmt <- "%Y-%m-%d %H:%M:%S"
                dcm[i, "format"] <- fmt
            }
            csv[,i] <- strftime(csv[,i], format=fmt)
        }
        write.dcf(x=dcm, file=f)

        # the data itself
        if(verbose) cat("OK\nwriting the dataset..")
        arg.lst <- list(x=csv, row.names=FALSE)
        f <- file.path(as.character(where), paste(target@name, ".csv", sep=""))
        arg.lst[["file"]] <- f
        enc <- target@dset.meta[["encoding"]]
        if(!is.null(enc)) arg.lst[["fileEncoding"]] <- enc 
        na.s <- target@dset.meta[["na.strings"]][[1]]
        if(!is.null(na.s)) arg.lst[["na"]] <- na.s
        
        type <- tolower(target@dset.meta[["type"]])
        if(is.null(type)) type <- "csv2"
        do.call(paste("write.", type, sep=""), arg.lst)
        if(verbose) cat("OK\n")
        
        file.path(as.character(where), paste(target@name, c(".dsm", ".dcm", ".csv"), sep=""))
    }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,CsvData,SftpLocation-method
setMethod(
    f="put",
    signature=c(target="CsvData", where="SftpLocation"),
    definition=function(target, where, verbose=TRUE, ...) {
        d <- tempfile()
        dir.create(d); on.exit(unlink(d, recursive=TRUE), add=TRUE)
        dl <- dirloc(d)
        put(target, dl)
        for (f in rownames(meta(dl))) put(file.path(d,f), where)
    }
)


