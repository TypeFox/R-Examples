#  sybilLogClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.

# sybilLogClass


#------------------------------------------------------------------------------#
#                      definition of the class sybilLog                        #
#------------------------------------------------------------------------------#

setOldClass("file")

setClass(Class = "sybilLog",
         representation(
             fh        = "file",
             fname     = "character",
             fpath     = "character",
             fenc      = "character",
             loglevel  = "integer",
             verblevel = "integer",
             lastStep  = "list",
             lstname   = "character",
             didFoot   = "logical"
         )
)

# derivatives
#setClass(Class = "sybilLog_MODEL", contains = "sybilLog")


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

# sybilLog <- function(loglevel, verblevel) {
#     if ( (missing(fluxes)) || (missing(verblevel)) ) {
#         stop("creating instance of class sybilLog needs loglevel and verblevel!")
#     }
#
#     new("sybilLog", loglevel = loglevel, verblevel = verblevel)
# }


sybilLog <- function(filename = NA,
                     filepath = ".",
                     loglevel = -1,
                     verblevel = 0,
                     logfileEnc = NA,
                     ...) {

    new(Class = "sybilLog",
        filename   = as.character(filename),
        filepath   = as.character(filepath),
        loglevel   = as.integer(loglevel),
        verblevel  = as.integer(verblevel),
        logfileEnc = as.character(logfileEnc),
        ...
    )
}


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "sybilLog",
          definition = function(.Object,
                                filename,
                                filepath,
                                loglevel,
                                verblevel,
                                logfileEnc, ...) {

    if ( (missing(loglevel)) || (missing(verblevel)) ) {
        warning("creation of instances of class sybilLog need \
                'loglevel' and 'logverb'; set to -1 and 0 respectively")
        loglevel  <- -1
        verblevel <- 0
    }

    if (missing(filename)) {
        filename <- NA
    }

    if (missing(filepath)) {
        filepath <- "."
    }

    if (missing(logfileEnc)) {
        logFileEnc <- NA
    }

    tst <- strftime(Sys.time(), format = "%Y%m%d%H%M%OS6")
    if (is.na(filename)) {
        filename <- paste("sybil_", tst, ".log", sep = "")
    }

    if (loglevel > -1) {
        tofile <- file.path(filepath, filename)

        if (is.na(logfileEnc)) {
            fenc <- getOption("encoding")
        }
        else {
            fenc <- logfileEnc
            options(useFancyQuotes = FALSE)
        }

        fh <- try(file(description = tofile, open = "wt",
                       encoding = fenc, ...), silent = TRUE)

        if (is(fh, "try-error")) {
            warning("can not write logfile")
            fh <- NA
        }
    }
    else {
        fh   <- NA
        fenc <- NA
    }

    if (!is.na(fh)) {
        .Object@fh    = fh
    }
    .Object@fname     = as.character(filename)
    .Object@fpath     = as.character(filepath)
    .Object@fenc      = as.character(fenc)
    .Object@loglevel  = as.integer(loglevel)
    .Object@verblevel = as.integer(verblevel)
    .Object@lastStep  = stinit(tst)
    .Object@lstname   = as.character(tst)
    #validObject(.Object)
    return(.Object)
}
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# fh
setMethod(f = "fh",
          signature = "sybilLog",
          definition = function(object) {
              return(object@fh)
          }
)

setReplaceMethod(f = "fh",
                 signature = "sybilLog",
                 definition = function(object, value) {
                     object@fh <- value
                     return(object)
                 }
)


# fname
setMethod(f = "fname",
          signature = "sybilLog",
          definition = function(object) {
              return(object@fname)
          }
)

setReplaceMethod(f = "fname",
                 signature = "sybilLog",
                 definition = function(object, value) {
                     object@fname <- value
                     return(object)
                 }
)


# fpath
setMethod(f = "fpath",
          signature = "sybilLog",
          definition = function(object) {
              return(object@fpath)
          }
)

setReplaceMethod(f = "fpath",
                 signature = "sybilLog",
                 definition = function(object, value) {
                     object@fpath <- value
                     return(object)
                 }
)


# fenc
setMethod(f = "fenc",
          signature = "sybilLog",
          definition = function(object) {
              return(object@fenc)
          }
)

setReplaceMethod(f = "fenc",
                 signature = "sybilLog",
                 definition = function(object, value) {
                     object@fenc <- value
                     return(object)
                 }
)


# loglevel
setMethod(f = "loglevel",
                 signature = "sybilLog",
          definition = function(object) {
              return(object@loglevel)
          }
)

setReplaceMethod(f = "loglevel",
                 signature = "sybilLog",
                 definition = function(object, value) {
                     object@loglevel <- value
                     return(object)
                 }
)


# verblevel
setMethod(f = "verblevel",
          signature = "sybilLog",
          definition = function(object) {
              return(object@verblevel)
          }
)

setReplaceMethod(f = "verblevel",
                 signature = "sybilLog",
                 definition = function(object, value) {
                     object@verblevel <- value
                     return(object)
                 }
)


# lstname
setMethod(f = "lstname",
          signature = "sybilLog",
          definition = function(object) {
              return(object@lstname)
          }
)


# didFoot
setMethod(f = "didFoot",
          signature = "sybilLog",
          definition = function(object) {
              return(object@didFoot)
          }
)

setReplaceMethod(f = "didFoot",
                 signature = "sybilLog",
                 definition = function(object, value) {
                     object@didFoot <- value
                     return(object)
                 }
)


# ---------------------------------------------------------------------------- #
# start/stop logging methods
# ---------------------------------------------------------------------------- #

# close log file
setReplaceMethod(f = "logClose",
                 signature = "sybilLog",
                 definition = function(object, value) {

    if (is(object@fh, "file")) {

        if (!isTRUE(didFoot(object))) {
            lc <- .printLogComment("end unexpected:")
            cat("\n", lc, sep = "", file = object@fh, append = TRUE)
        }

        check <- try(isOpen(object@fh), silent = TRUE)
        if ( ! is(check, "try-error")) {
            close(object@fh)
        }
    }

    if (stexists(object@lstname)) {
        stclear(object@lstname)
    }

    return(object)

}
)


# log file header
setMethod(f = "logHead",
          signature = "sybilLog",
          definition = function(object) {

    if (!is(object@fh, "file")) {
        return(FALSE)
    }

    lc <- .printLogComment("start:")
    cat(lc, file = object@fh, append = TRUE)

    return(TRUE)
}
)


# log file footer
setReplaceMethod(f = "logFoot",
                 signature = "sybilLog",
                 definition = function(object, value) {

    if (object@verblevel > 1) {
        message("Done.")
    }

    if (!is(object@fh, "file")) {
        return(object)
    }

    lc <- .printLogComment("end:")
    cat("\n", lc, sep = "", file = object@fh, append = TRUE)

    object@didFoot <- value

    return(object)
}
)


# ---------------------------------------------------------------------------- #
# errors, warnings and messages
# ---------------------------------------------------------------------------- #

# error message
setMethod(f = "logError",
          signature = "sybilLog",
          definition = function(object, msg, num) {

    errmsg <- gettext(paste(msg, collapse = " "))

    if (is(object@fh, "file")) {
        cat("E\t", errmsg, "\t", date(), file = object@fh, append = TRUE)
    }

    err <- sybilError(errmsg = errmsg)

    return(err)
}
)

# error message
setMethod(f = "logError",
          signature(object = "sybilLog", num = "numeric"),
          definition = function(object, msg, num) {

    errmsg <- gettext(paste(msg, collapse = " "))

    if (is(object@fh, "file")) {
        cat("E\t", num, "\t", errmsg, "\t", date(),
            file = object@fh, append = TRUE)
    }
    else {}

    err <- sybilError(errmsg = errmsg, number = num)

    return(err)
}
)


# warning
setMethod(f = "logWarning",
          signature = "sybilLog",
          definition = function(object, ...) {

    msg <- gettext(paste(..., collapse = " "))

    if (object@verblevel > 0) {
        warning(msg, call. = FALSE)
    }
    else {}

    if (object@loglevel > 0) {
        if (is(object@fh, "file")) {
            cat("W\t", msg, "\n", file = object@fh, append = TRUE)
        }
        else {}
    }
    else {}

    return(TRUE)
}
)


# message
setMethod(f = "logMessage",
          signature = "sybilLog",
          definition = function(object, appendEllipsis = FALSE, ...) {

    msg <- gettext(paste(..., collapse = " "))

    if (object@verblevel > 1) {
        if (isTRUE(appendEllipsis)) {
            message(msg, " ... ", appendLF = FALSE)
        }
        else {
            message(msg, appendLF = FALSE)
        }
    }

    if (object@loglevel > 1) {
        if (is(object@fh, "file")) {
            cat("M\t", msg, "\n", file = object@fh, append = TRUE)
        }
    }

    return(TRUE)
}
)


# ---------------------------------------------------------------------------- #
# other logging methods
# ---------------------------------------------------------------------------- #

# return TRUE if object@fh is a connection (file)
setMethod(f = "logFH",
          signature = "sybilLog",
          definition = function(object) {

    out <- is(object@fh, "file")

    return(out)
}
)


# comments to the logfile
setMethod(f = "logComment",
          signature = "sybilLog",
          definition = function(object, cmt, cmtChar) {

    msg <- gettext(paste(cmt, collapse = " "))

    if (missing(cmtChar)) {
        cmtChar <- "# "
    }

    if (object@verblevel > 2) {
        cat(cmtChar, msg, "\n", sep = "")
    }

    if (object@loglevel > 2) {
        if (is(object@fh, "file")) {
            cat(cmtChar, msg, "\n",
                sep = "", file = object@fh, append = TRUE)
        }
    }

    return(TRUE)
}
)


# results of optimization
setMethod(f = "logOptimizationTH",
          signature = "sybilLog",
          definition = function(object) {

    th <- "opt no.  | ret | stat | obj value | dir | obj c | flux no. \n"
    
    if (object@verblevel > 2) {
        cat(th)
    }

    if (object@loglevel > 2) {
        if (is(object@fh, "file")) {
            cat(th, file = object@fh, append = TRUE)
        }
    }

    return(TRUE)
}
)


# results of optimization
setMethod(f = "logOptimization",
          signature = "sybilLog",
          definition = function(object, ok, stat, obj, dir, objc, del, i) {

    if ( (object@verblevel > 2) || (object@loglevel > 2) ) {
        fi    <- sprintf("  %-6s", paste("[", i, "]", sep = ""))
        fok   <- sprintf("%-3i", ok)
        fstat <- sprintf("%-4i", stat)
        fobj  <- sprintf("%9.3f", obj)
        if (is.null(dir)) {
            fdir  <- "   "
        }
        else {
            fdir  <- sprintf("%3s", dir)
        }
        if (is.null(objc)) {
            fobjc  <- "     "
        }
        else {
            fobjc  <- sprintf("%5.1f", objc)
        }
        fdel  <- paste(del, collapse = " ")
    
        prstr <- paste(fi, fok, fstat, fobj, fdir, fobjc, fdel, sep = " | ")

        if (object@verblevel > 2) {
            cat(prstr, "\n", sep = "")
        }

        if (object@loglevel > 2) {
            if (is(object@fh, "file")) {
                cat(prstr, "\n", sep = "", file = object@fh, append = TRUE)
            }
        }
    }

    return(invisible(TRUE))
}
)


# performing step foo
setReplaceMethod(f = "logStep",
                 signature = "sybilLog",
                 definition = function(object, value) {

    if (is.na(value)) {
        didstep <- stpop(object@lstname)
        if (object@verblevel > 1) {
            if (stlength(object@lstname) < 1) {
                message("OK")
            }
        }

        if (object@loglevel > 1) {
            if (is(object@fh, "file")) {
                write(paste("# done", didstep),
                      file = object@fh, append = TRUE)
            }
        }

    }
    else {

        msg <- gettext(paste(value, collapse = " "))
        stpush(object@lstname, msg)

        if (object@verblevel > 1) {
            if (pmatch(gettext("FAILED"), msg, nomatch = 0) == 1) {
                message(msg, appendLF = TRUE)
            }
            else {
                message(paste(msg, "... "), appendLF = FALSE)
            }
        }

        if (object@loglevel > 1) {
            if (is(object@fh, "file")) {
                write(paste("#", msg), file = object@fh, append = TRUE)
            }
        }
    }

    return(object)
}
)


# log function call
setMethod(f = "logCall",
          signature = "sybilLog",
          definition = function(object, nog) {

      if (missing(nog)) {
          nog <- 1
      }
      else {}

      if (object@loglevel > 2) {
          if (is(object@fh, "file")) {
              fc <- sys.call(sys.parent(n = nog))
              cat("# call to function", dQuote(fc[[1]]),
                  "with arguments:\n", file = object@fh, append = TRUE)
              .printNamedList(nList = as.list(fc)[-1],
                                      file = object@fh, append = TRUE)
          }
          else {}
      }
      else {}

      return(TRUE)
}
)


# setMethod("logCall", signature(object = "sybilLog"),
#           function(object, func, fargl, thdargs) {
#               if (object@loglevel > 2) {
#                   if (is(object@fh, "file")) {
#                       fc <- as.character(func)
#                       cat("# formal arguments to ", fc, "()\n",
#                           sep = "", file = object@fh, append = TRUE)
#                       .printNamedList(fargl, file = object@fh, append = TRUE)
#
#                       if ( (length(thdargs) > 0) && (!is.na(thdargs)) ) {
#                           cat("# further arguments to", fc, "(...)\n",
#                               file = object@fh, append = TRUE)
#                           if (length(thdargs) > 0) {
#                               .printNamedList(thdargs,
#                                                       file = object@fh, append = TRUE)
#                           }
#                           else {
#                               cat("none\n", file = object@fh, append = TRUE)
#                           }
#                       }
#                   }
#               }
#
#               return(TRUE)
#           }
# )

## !! required for the above method !! ##
#     # log function call
#     if (loglevel > 2) {
#         fargs <- formals()
#         for (na in names(fargs)) {
#             val <- try(eval(parse(text=na)), silent = TRUE)
#             if (!is(val, "try-error")) {
#                 fargs[na] <- val
#             }
#         }
#         print(match.call())
#         logCall(logObj, sys.call(), fargs, list(...))
#     }
