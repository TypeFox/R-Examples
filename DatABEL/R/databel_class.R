#' DatABEL class
#'
#' DatABEL stores matrix-shape data in such a way that it can be
#' retrieved fast.
#'
#' @slot usedRowIndex (\code{"integer"})
#' @slot usedColIndex (\code{"integer"})
#' @slot uninames (\code{"list"})
#' @slot backingfilename Name of the (stem of the) file that contains
#'       the data stored in DatABEL format (\code{"character"})
#' @slot cachesizeMb Amount (in MB) of RAM to use for caching DatABEL
#'       data. (\code{"integer"})
#' @slot data (\code{"externalptr"})
#' @note Will extend description here
#' @name databel-class
#' @rdname databel-class
#' @aliases databel-class
#' @exportClass databel
#' @import methods
#' @author Yurii Aulchenko
#'

#
# databel  R class
# (C) 2009, 2010, Yurii Aulchenko, EMCR
#
# This class works with datable_cpp object using direct access
#
# SLOTS
#
# METHODS
# initilize, new
# [      keeping class
# [<-    keeping class
# dim
# length
# dimnames
# show
#
# !cbind
# !rbind
#
# is()
# as("databel","matrix")
# as("databel","vector")
# as("matrix","databel")
#
# backingfilename
# cachesizeMb
# cachesizeMb<-
# !get_data_type
# get_dimnames
# set_dimnames<-
# connect
# disconnect
# save_as
# dimnames<-
#
#

setClass(
    Class = "databel",
    representation = representation(
        usedRowIndex = "integer",
        usedColIndex = "integer",
        uninames = "list",
        backingfilename = "character",
        cachesizeMb = "integer",
        data = "externalptr"
        ),
    package = "DatABEL"
    );


setMethod(
    f = "initialize",
    signature = "databel",
    definition = function(.Object, baseobject, cachesizeMb=64, readonly=TRUE)
    {
#			cat("----- databel ini start -----\n");

        if (!(is(baseobject, "character") || is(baseobject, "databel"))) {
            stop("databel initialize: baseobject should be of character (filename) or databel class");
        }

        if (!is.numeric(cachesizeMb)) {
            stop("databel initialize: cache size must be numeric")
        }

        readonly <- as(readonly, "logical")
        cachesizeMb <- as(cachesizeMb, "integer")
        if (cachesizeMb < 0) {
            stop(paste("databel initialize: cache size must be positive integer; now",
                       cachesizeMb))
        }

        if (is(baseobject, "character")){
            address <- .Call("open_FilteredMatrix_R",
                             fname = baseobject,
                             csize = cachesizeMb,
                             rof = readonly,
                             PACKAGE="DatABEL");
            if (is(address, "null")) {
                stop("databel initialize: can not create databel object at step 1, NULL pointer returned")
            }

            .Object@data <- address

            nrows <- .Call("get_nobs_R", .Object@data, PACKAGE="DatABEL")
            ncols <- .Call("get_nvars_R", .Object@data, PACKAGE="DatABEL")
            .Object@usedRowIndex <- c(1:nrows)
            .Object@usedColIndex <- c(1:ncols)

            .Object@uninames <- uninames(.Object@data)

            .Object@backingfilename <- baseobject
            .Object@cachesizeMb <- cachesizeMb

        } else if (is(baseobject, "databel")) {

            for (sn in slotNames(baseobject))
                slot(.Object, sn) <- slot(baseobject, sn)

            address <- .Call("open_FilteredMatrix_R",
                             fname = backingfilename(baseobject),
                             csize = cachesizeMb,
                             rof = readonly,
                             PACKAGE="DatABEL");
            if (is(address, "null")) {
                stop("databel initialize: can not create databel object at step 2-1, NULL pointer returned")
            }

            .Object@data <- address
            .Object@data <- .Call("setFilteredArea_R",
                                  .Object@data,
                                  .Object@usedColIndex,
                                  .Object@usedRowIndex);

        } else {
            stop("databel initialize: unreachable statement -- baseobject should be of databel class of FV-file name");
        }

        #			cat("----- databel ini end -------\n");

        return(.Object)
    }
    );

# replace standard methods

#' @rdname  databel-class
#' @export
setMethod(
    f = "show",
    signature = "databel",
    definition = function(object)
    {
        connected <- databel_check(object, reconnect=TRUE)
        # SHOULD ACTUALLY SHOW ONLY NON-MASKED DATA

        cat("uninames$unique.names =", object@uninames$unique.names, "\n")
        cat("uninames$unique.rownames =", object@uninames$unique.rownames, "\n")
        cat("uninames$unique.colnames =", object@uninames$unique.colnames, "\n")
        cat("backingfilename =", object@backingfilename, "\n")
        cat("cachesizeMb =", object@cachesizeMb, "\n")
        cat("number of columns (variables) = ", ncol(object), "\n");
        cat("number of rows (observations) = ", nrow(object), "\n");
        toCol <- 10
        toRow <- 5
        if (ncol(object)<toCol) toCol <- ncol(object)
        if (nrow(object)<toRow) toRow <- nrow(object)
        cat("usedRowIndex: ")
        for (i in 1:toRow) cat(object@usedRowIndex[i], " ")
        if (toRow<dim(object)[1]) cat("...")
        cat("\n")
        cat("usedColIndex: ")
        for (i in 1:toCol) cat(object@usedColIndex[i], " ")
        if (toCol<dim(object)[2]) cat("...")
        cat("\n")

        if (!connected) {
            cat("databel show: object is not connected\n")
            return();
        }
        cat("Upper-left", toCol, "columns and ", toRow, "rows:\n")
        showout <- as(object[1:toRow, 1:toCol], "matrix")
        print(showout)
    }
    );


#' @aliases dim,databel-method
#' @rdname databel-class
setMethod(
    f = "dim",
    signature = "databel",
    definition = function(x)
    {
        connected <- databel_check(x);
        if (!connected) stop("databel dim: databel_check failed")
        return(c(length(x@usedRowIndex), length(x@usedColIndex)))
    }
    );


#' @aliases length,databel-method
#' @rdname databel-class
setMethod(
    f = "length",
    signature = "databel",
    definition = function(x)
    {
        connected <- databel_check(x)
        if (!connected) stop("databel length: databel_check failed")
        dm <- dim(x)
        return(dm[1]*dm[2])
    }
    );


#' @rdname databel-class
setMethod(
    f = "dimnames",
    signature = "databel",
    definition = function(x)
    {
        connected <- databel_check(x)
        if (!connected) stop("databel dimnames: databel_check failed")
        if (x@uninames$unique.names) {
            return(get_dimnames(x))
        } else if (x@uninames$unique.rownames) {
            return(list(get_dimnames(x)[[1]], NULL))
        } else if (x@uninames$unique.colnames) {
            return(list(NULL, get_dimnames(x)[[2]]))
        } else {
            return(NULL)
        }
    }
    );


#' @rdname databel-class
#' @param value Values to be replaced/inserted
setMethod(
    f = "dimnames<-",
    signature = "databel",
    definition = function(x, value)
    {
        connected <- databel_check(x)
        if (!connected) stop("databel dimnames<-: databel_check failed")
        if (anyDuplicated(value[[1]])) stop("non-unigue names in dim [[1]] (use set_dimnames?)")
        if (anyDuplicated(value[[2]])) stop("non-unigue names in dim [[2]] (use set_dimnames?)")
        set_dimnames(x) <- value
        x@uninames <- uninames(x@data)
        return(x)
    }
    );


#' @rdname databel-class
#' @param j Column index
#' @param drop Boolean (FALSE by default); UNUSED
setMethod(
    f = "[",
    signature = "databel",
    definition = function(x, i, j, drop)
    {
        #			print("[ started")
        connected <- databel_check(x)
        if (!connected) stop("databel [: object is not connected")
        if (missing(drop)) drop = FALSE;
        newi <- convert_intlogcha_index_to_int(i, x, 1)
        newj <- convert_intlogcha_index_to_int(j, x, 2)
        out <- databel(x)
        #			print(c("out dims orig are", dim(out)))
        out@usedRowIndex <- out@usedRowIndex[newi]
        out@usedColIndex <- out@usedColIndex[newj]
        out@data <- .Call("setFilteredArea_R", out@data,
                          out@usedColIndex, out@usedRowIndex);
        #			print(c("out dims after are", dim(out)))
        out@uninames <- uninames(out@data)
        #			print("[ ended")
        return(out);
    }
    );


#' @rdname databel-class
#' @param x A DatABEL object
#' @param i Row index
setMethod(
    f = "[<-",
    signature = "databel",
    definition = function(x, i, j, value)
    {
                                        #			print("started [<-")
        connected <- databel_check(x)
        if (!connected) stop("databel [<-: databel_check failed")

        newi <- convert_intlogcha_index_to_int(i, x, 1)
        newj <- convert_intlogcha_index_to_int(j, x, 2)

        if (length(value) != length(newi)*length(newj)) {
            stop("databel [<-: dimensions of i, j, value do not match")
        }

        #			value <- matrix(value, ncol=length(newj), nrow=length(newi))
        #            print(newi);
        #            print(newj);

        if(!.Call("assignDoubleMatrix", x@data, newi, newj,
                  as.double(value), as.integer(0))) {
            stop("databel [<-: can't write variable.");
        }
        #			print("finished [<-")
        return(x)
    }
    );


# IS

setGeneric('is.databel', function(x) standardGeneric('is.databel'))

setMethod('is.databel', signature(x='databel'),
          function(x) return(TRUE))

setMethod('is.databel', definition=function(x) return(FALSE))

                                        # setAs

as.matrix.databel <- function(x, ... )
{
    return(as(x, "matrix"))
}

as.vector.databel <- function(x, ... )
{
    return(as(x, "vector"))
}

#' @export
as.double.databel <- function(x, ... )
{
    to <- as(x, "matrix")
    return(as(to, "double"))
}



setAs("databel", "vector",
      function(from) {
          to <- as(from, "matrix")
          return(as(to, "vector"))
      }
      );

setAs("databel", "matrix",
      function(from) {
          connected <- databel_check(from)
          if (!connected) stop("setAs('databel', 'matrix'): check_connected failed")
          return(databel2matrix(from))
      }
      );


setAs("matrix", "databel",
      function(from) {
                                        #print("as matrix->databel begin");
          if (!is.numeric(from)) stop("from must be numeric (integer or double)")
          type <- "DOUBLE"
          tmpfilename <- get_temporary_file_name();
          to <- matrix2databel(from, filename=tmpfilename, cachesizeMb=64, type=type)
          cat("coersion from 'matrix' to 'databel' of type", type, "; object connected to file", tmpfilename, "\n")
                                        #print("as matrix->databel end");
          return(to)
      }
      );



### new generics

#' @aliases get_dimnames,databel-method
#' @rdname databel-class
#' @export
setGeneric(
    name = "get_dimnames",
    def = function(object) {standardGeneric("get_dimnames");}
    );

setMethod(
    f = "get_dimnames",
    signature = "databel",
    definition = function(object)
    {
        connected <- databel_check(object)
        if (!connected) stop("object is not connected", immediate.=TRUE)
        return(list(.Call("get_all_obsnames_R", object@data,
                          PACKAGE="DatABEL"),
                    .Call("get_all_varnames_R", object@data,
                          PACKAGE="DatABEL")))
    }
    );


#' @aliases set_dimnames<-,databel-method
#' @rdname databel-class
#' @export
setGeneric(
    name = "set_dimnames<-",
    def = function(x, value) {standardGeneric("set_dimnames<-");}
    );


setMethod(
    f = "set_dimnames<-",
    signature = "databel",
    definition = function(x, value)
    {

        connected <- databel_check(x);
        if (!connected) stop("set_dimnames<-: databel_check failed")

        if (!is.list(value)) stop("set_dimnames<-: value is not a list")
        if (length(value)!=2) {
            stop("set_dimnames<-: value should be a list with two vectors")
        }

        if (length(value[[1]]) != dim(x)[1]) {
            if (is.null(value[[1]])) {
                value[[1]] <- as.character(c(1:dim(x)[1]))
            } else {
                stop("set_dimnames<-: dimention 1 of x and lengthof list[[1]] do not match")
            }
        }

        if (length(value[[2]]) != dim(x)[2]) {
            if (is.null(value[[2]])) {
                value[[2]] <- as.character(c(1:dim(x)[2]))
            } else {
                stop("set_dimnames<-: dimention 2 of x and lengthof list[[2]] do not match")
            }
        }

        if (!is.character(value[[1]])) {
            stop("set_dimnames<-: colnames must be characters")
        }
        if (!is.character(value[[2]])) {
            stop("set_dimnames<-: rownames must be characters")
        }
        if (!is.null(value[[2]])) {
            r1 <- .Call("set_all_varnames_R", x@data,
                        as.character(value[[2]]))
        }
        if (!is.null(value[[1]])) {
            r2 <- .Call("set_all_obsnames_R", x@data,
                        as.character(value[[1]]))
        }

        x@uninames <- uninames(x@data)

        # if (length(unique(value[[1]]))==dim(x)[1])
        # x@uninames$unique.rownames <- TRUE else x@uninames$unique.rownames <- FALSE
        # if (length(unique(value[[2]]))==dim(x)[2])
        # x@uninames$unique.colnames <- TRUE else x@uninames$unique.colnames <- FALSE
        # if (x@uninames$unique.colnames && x@uninames$unique.rownames)
        # x@uninames$unique.names <- TRUE else x@uninames$unique.names <- FALSE

        return(x)
    }
    );


#' @aliases backingfilename,databel-method
#' @rdname databel-class
#' @param object A DatABEL object
#' @export
setGeneric(
    name = "backingfilename",
    def = function(object) {standardGeneric("backingfilename");}
    );

setMethod(
    f = "backingfilename",
    signature = "databel",
    definition = function(object)
    {
        return(object@backingfilename)
    }
    );


#' @aliases cachesizeMb,databel-method
#' @rdname databel-class
#' @export
setGeneric(
    name = "cachesizeMb",
    def = function(object) {standardGeneric("cachesizeMb");}
    );

setMethod(
    f = "cachesizeMb",
    signature = "databel",
    definition = function(object)
    {
        return(object@cachesizeMb)
    }
    );


#' @aliases cachesizeMb<-,databel-method
#' @rdname databel-class
#' @export
setGeneric(
    name = "cachesizeMb<-",
    def = function(x, value) {standardGeneric("cachesizeMb<-");}
    );

setMethod(
    f = "cachesizeMb<-",
    signature = "databel",
    definition = function(x, value)
    {
        # cat("set_cachesizeMb not implemented yet, leaving cachesizeMb unchanged\n")
        connected <- databel_check(x)
        if (!connected) stop("cachesizeMb<-: databel_check failed")
        if (!is.numeric(value)) stop("value must be numeric")
        value <- as.integer(value)
        if (value < 0) stop("can not set cachesizeMb to <0")
        .Call("set_cachesizeMb_R", x@data, value)
        x@cachesizeMb <- .Call("get_cachesizeMb_R", x@data)
        return(x)
    }
    );


#' @aliases save_as,databel-method
#' @rdname databel-class
#' @export
#' @param file Filename to save to
#' @param rows Index for the rows
#' @param cols Index for the columns
#' @param cachesizeMb Amount (in MB) of RAM to use for caching DatABEL
#' data.
#' @param readonly Boolean that specifies whether the file is to be
#' used in read-only mode or not
setGeneric(
    name = "save_as",
    def = function(x, rows, cols, file, cachesizeMb=64, readonly=TRUE)
    {
        standardGeneric("save_as");
    }
    );


setMethod(
    f = "save_as",
    signature = "databel",
    definition = function(x, rows, cols, file, cachesizeMb=64, readonly=TRUE)
    {
                                        #allowd_types <- c("databel", "text")
        if (!is.character(file)) {
            stop("databel save_as: file argument should be character")
        }
        if (!missing(rows)) {
            newi <- convert_intlogcha_index_to_int(rows, x, 1)
        } else {
            newi <- 1:dim(x)[1]
        }
        if (!missing(cols)) {
            newj <- convert_intlogcha_index_to_int(cols, x, 2)
        } else {
            newj <- 1:dim(x)[2]
        }

### check order!!!
        # print(newj)
        # print(newi)
        intpar <- as.integer(c(length(newj), length(newi), (newj-1), (newi-1)))
        # print(intpar)
        if (!.Call("save_R", file, intpar, x@data)) {
            stop("can not save_as(): save_R failed")
        }
        newobj <- databel(file, cachesizeMb=cachesizeMb, readonly=readonly)
        return(newobj)
    }
    );


#' @aliases connect,databel-method
#' @rdname databel-class
#' @export
setGeneric(
    name = "connect",
    def = function(object, readonly=TRUE) {standardGeneric("connect");}
    );

setMethod(
    f = "connect",
    signature = "databel",
    definition = function(object, readonly=TRUE)
    {
        # print("connect for databel started")
        connected <- databel_check(object, reconnect = FALSE,
                                   stop_on_error = FALSE, quiet =
                                   TRUE)
        if (connected) {
            # print("connect for databel finished")
            warning("object already connected, nothing done") #, immediate. = TRUE)
            return();
        } else {
            new_obj <- databel(backingfilename(object),
                               cachesizeMb=cachesizeMb(object),
                               readonly=readonly)
            new_obj@usedColIndex <- object@usedColIndex
            new_obj@usedRowIndex <- object@usedRowIndex
            res <- .Call("setFilteredArea_R", new_obj@data,
                         object@usedColIndex, object@usedRowIndex);
            # print("after setFilteredArea_R in connect")
            # print("before eval.parent(sub...")
            eval.parent(substitute(object <- new_obj));
            # print("after eval.parent(sub...")
        }
        # print("connect for databel finished")
    }
    );


#' @rdname databel-class
#' @aliases disconnect,databel-method
#' @export
setGeneric(
    name = "disconnect",
    def = function(object) {standardGeneric("disconnect");}
    );

setMethod(
    f = "disconnect",
    signature = "databel",
    definition = function(object)
    {
        connected <- databel_check(object, reconnect = FALSE, stop_on_error = FALSE)
        if (!connected) {
            warning("object is already disconnected")
            return()
        }
        tmp0 <- .Call("disconnectFilteredAndAbstract_R", object@data,
                      PACKAGE="DatABEL");
        tmp1 <- gc();
    }
    );

#' @aliases setReadOnly<-,databel-method
#' @rdname databel-class
#' @export
setGeneric(
    name = "setReadOnly<-",
    def = function(x, value) {standardGeneric("setReadOnly<-");}
    );

setMethod(
    f = "setReadOnly<-",
    signature = "databel",
    definition = function(x, value)
    {
        connected <- databel_check(x)
        if (!connected) stop("databel setReadOnly: object is not connected")
        res <- .Call("setReadOnly_R", x@data, as(value, "logical"),
                     PACKAGE="DatABEL");
        if (!res) {
            stop(paste("databel setReadOnly: can not set ReadOnly flag to",
                       as(value, "logical")))
            }
        return(x)
    }
    );
