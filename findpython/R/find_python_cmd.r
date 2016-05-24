# Copyright (c) 2012-2014 Trevor L. Davis
# Copyright (c) 2014 Paul Gilbert
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#' Tests whether the python command is sufficient
#'
#' \code{is_python_sufficient} checks whether a given python binary has all the
#' desired features (minimum and/or maximum version number and/or access to
#' certain modules).
#'
#' @param path The path to a given python binary.  
#'      If binary is on system path just the binary name will work.
#' @param minimum_version The minimum version of python it should be.  
#'      Should be a string with major and minor number separated by a '.'.
#'      If left NULL won't impose such a restriction.  
#' @param maximum_version The maximum version of python it should be. 
#'      Should be a string with major and minor number separated by a '.'.
#'      If left NULL won't impose such a restriction.
#' @param required_modules Which modules should be required.  
#'      Can use a single "|" to represent a single either-or requirement like "json|simplejson".
#'      If left NULL won't impose such a restriction.
#' @return \code{TRUE} or \code{FALSE} depending on whether the python binary met all requirements
#' @export
is_python_sufficient <- function(path, minimum_version=NULL,
                                 maximum_version=NULL, required_modules=NULL) {
    python_code <- vector("character")
    if(!is.null(required_modules)) {
        import_code <- .create_import_code(required_modules)
        python_code <- append(python_code, import_code)
    }
    version_code <- .create_version_checking_code(minimum_version, maximum_version)
    python_code <- append(python_code, version_code)
    ok_message <- "Everything worked out okay"
    python_code <- append(python_code, paste("print(", shQuote(ok_message), ")", sep=""))
    tryCatch({
            output <- system(path, intern=TRUE, input=python_code, ignore.stderr=TRUE)
            any(grepl(ok_message, output))
        }, 
        warning = function(w) { 
            # qpath <- sQuote(path)
            # warning(qpath, 
            #     "Doesn't seem to have required modules or is of insufficient version")
            FALSE
        },
        error = function(e) {
            FALSE
        })
}

#' Find a suitable python cmd or give error if not possible
#'
#' \code{find_python_cmd} finds a suitable python cmd or raises an error if not possible
#'
#' @inheritParams is_python_sufficient
#' @param error_message What error message the user will see if couldn't find a sufficient python binary.
#'     If left NULL will print out a default message.
#' @return The path to an appropriate python binary.  If such a path wasn't found then it will throw an error.
#' @examples 
#'      \dontrun{
#'              find_python_cmd()
#'              find_python_cmd(minimum_version='2.6', maximum_version='2.7')
#'              find_python_cmd(required_modules = c('argparse', 'json | simplejson'))
#'      }
#' @seealso \code{\link{can_find_python_cmd}} for a wrapper which doesn't throw an error
#' @export
find_python_cmd <- function(minimum_version=NULL, maximum_version=NULL,
                            required_modules=NULL, error_message=NULL) {
    python_cmds <- c(getOption("python_cmd", ""), "python", Sys.getenv("PYTHON", ""), 
                    sprintf("python%.1f", c(seq(3.4, 3.0, by=-0.1), seq(2.7, 2.0, by=-0.1))), 
                    "python3", "python2", "pypy",
                    Sys.getenv("PYTHON3", ""), Sys.getenv("PYTHON2", ""), 
                    sprintf("C:/Python%s/python", c(34:30, 27:20)))
    python_cmds <- unique(python_cmds)
    python_cmds <- Sys.which(python_cmds)
    python_cmds <- python_cmds[which(python_cmds != "")]
    for(cmd in python_cmds) {
        if(is_python_sufficient(cmd, minimum_version, maximum_version, required_modules)) { 
            return(cmd)
        }
    }
    if(is.null(error_message)) { error_message <- paste("Couldn't find a sufficient Python binary.",
                    "If you haven't installed the Python dependency yet please do so.",
                    "If you have but it isn't on the system path (as is default on Windows) please add it to path",
                    "or set options('python_cmd'='/path/to/binary') ",
                    "or set the PYTHON, PYTHON2, or PYTHON3 environmental variables.",
                    if(!is.null(minimum_version)) paste('Python must be at least version', minimum_version),
                    if(!is.null(maximum_version)) paste('Python must be at most version', maximum_version),
                    if(!is.null(required_modules)) paste('Python must have access to the modules:', 
                                                         paste(required_modules, collapse=', ')))}
    stop(error_message)
}
#' Determins whether or not it can find a suitable python cmd 
#'
#' \code{can_find_python_cmd} runs \code{find_python_cmd} and returns whether it could find a suitable python cmd.  If it was successful its output also saves the found command as an attribute. 
#'
#' @inheritParams find_python_cmd
#' @param silent Passed to \code{try}, whether any error messages from \code{find_python_cmd} should be suppressed
#' @return \code{TRUE} or \code{FALSE} depending on whether \code{find_python_cmd} could find an appropriate python binary.
#'     If \code{TRUE} the path to an appropriate python binary is also set as an attribute.
#' @examples 
#'      did_find_cmd <- can_find_python_cmd()
#'      python_cmd <- attr(did_find_cmd, "python_cmd")
#' @seealso \code{\link{find_python_cmd}}
#' @export
can_find_python_cmd <- function(minimum_version = NULL,
          maximum_version = NULL, required_modules = NULL, 
          error_message = NULL, silent = FALSE) {
    python_cmd <- try(find_python_cmd(minimum_version = minimum_version,
                          maximum_version = maximum_version,
                          required_modules = required_modules,
                          error_message = error_message),
                    silent = silent)
    if(inherits(python_cmd, "try-error")) {
        r <- FALSE
    } else {
        r <- TRUE
        attr(r, 'python_cmd') <- python_cmd
    }
    r
}

# Create appropriate module import code
.create_import_code <- function(required_modules) {
    import_code <- vector("character")
    for(module in required_modules) {
        if(grepl("\\|", module)) {
            module <- gsub(" ", "", module)
            module <- strsplit(module, "\\|")[[1]]
            import_code <- append(import_code, 
                                  c("try:", 
                                   paste("    import", module[1]),
                                   "except ImportError:", 
                                   paste("    import", module[2])))
        } else {
            import_code <- append(import_code, paste("import", module))
        }
    }
    import_code
}

# Create appropriate version checking code
.create_version_checking_code <- function(minimum_version = NULL, maximum_version = NULL) {
    if(is.null(minimum_version) && is.null(maximum_version)) 
        return(c())
    else {
        version_code <- c("import sys")
    }
    if(!is.null(minimum_version)) {
        min_version <- strsplit(minimum_version, '\\.')[[1]]
        major = min_version[1]
        minor = min_version[2]
        version_code <- append(version_code, 
                           c(paste("if sys.version_info.major <", major, 
                                 ": raise Exception('Major version too low')"),
                           paste("if sys.version_info.major ==", major, 
                                 "and sys.version_info.minor <", minor, 
                                 ": raise Exception('Minor version too low')")))
    }
    if(!is.null(maximum_version)) {
        max_version <- strsplit(maximum_version, '\\.')[[1]]
        major = max_version[1]
        minor = max_version[2]
        version_code <- append(version_code, 
                           c(paste("if sys.version_info.major >", major, 
                                 ": raise Exception('Major version too high')"),
                           paste("if sys.version_info.major ==", major, 
                                 "and sys.version_info.minor >", minor, 
                                 ": raise Exception('Minor version too high')")))
    }
    version_code
}
