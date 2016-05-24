# Interpreted GRASS 7 interface functions
# Copyright (c) 2015 Roger S. Bivand
#

set.ignore.stderrOption <- function(value) {
	if (!is.logical(value)) stop ("logical argument required")
	res <- get("ignore.stderr", envir = .GRASS_CACHE)
	assign("ignore.stderr", value, envir = .GRASS_CACHE)
	res
}

get.ignore.stderrOption <- function() {
	get("ignore.stderr", envir = .GRASS_CACHE)
}

set.stop_on_no_flags_parasOption <- function(value) {
	if (!is.logical(value)) stop ("logical argument required")
	res <- get("stop_on_no_flags_paras", envir = .GRASS_CACHE)
	assign("stop_on_no_flags_paras", value, envir = .GRASS_CACHE)
	res
}

get.stop_on_no_flags_parasOption <- function() {
	get("stop_on_no_flags_paras", envir = .GRASS_CACHE)
}

set.useGDALOption <- function(value) {
	if (!is.logical(value)) stop ("logical argument required")
	res <- get("useGDAL", envir = .GRASS_CACHE)
	assign("useGDAL", value, envir = .GRASS_CACHE)
        if (value) {
          if (requireNamespace("rgdal", quietly = TRUE)) {
	    assign("useGDAL", value, envir = .GRASS_CACHE)
          } else {
            warning("rgdal not available")
          }
        }
# require(rgdal)
	res
}

get.useGDALOption <- function() {
	get("useGDAL", envir = .GRASS_CACHE)
}

set.pluginOption <- function(value) {
	res <- get("plugin", envir = .GRASS_CACHE)
        if (is.null(value)) {
	    assign("plugin", value, envir = .GRASS_CACHE)
        } else if (is.logical(value)) {
	    assign("plugin", value, envir = .GRASS_CACHE)
        } else stop ("logical or NULL argument required")
	res
}

get.pluginOption <- function() {
	get("plugin", envir = .GRASS_CACHE)
}

set.echoCmdOption <- function(value) {
	if (!is.logical(value)) stop ("logical argument required")
	res <- get("echoCmd", envir = .GRASS_CACHE)
	assign("echoCmd", value, envir = .GRASS_CACHE)
	res
}

get.echoCmdOption <- function() {
	get("echoCmd", envir = .GRASS_CACHE)
}

set.useInternOption <- function(value) {
	if (!is.logical(value)) stop ("logical argument required")
	res <- get("useIntern", envir = .GRASS_CACHE)
	assign("useIntern", value, envir = .GRASS_CACHE)
	res
}

get.useInternOption <- function() {
	get("useIntern", envir = .GRASS_CACHE)
}

set.legacyExecOption <- function(value) {
	if (!is.logical(value)) stop ("logical argument required")
	res <- get("legacyExec", envir = .GRASS_CACHE)
	assign("legacyExec", value, envir = .GRASS_CACHE)
	res
}

get.legacyExecOption <- function() {
	get("legacyExec", envir = .GRASS_CACHE)
}

set.defaultFlagsOption <- function(value) {
	res <- get("defaultFlags", envir = .GRASS_CACHE)
        if (is.null(value)) {
	    assign("defaultFlags", value, envir = .GRASS_CACHE)
        } else {
	    if (!is.character(value)) stop ("character argument required")
            stopifnot(length(value) > 0)
            stopifnot(all(value %in% c("quiet", "verbose")))
	    assign("defaultFlags", value, envir = .GRASS_CACHE)
        }
	res
}

get.defaultFlagsOption <- function() {
	get("defaultFlags", envir = .GRASS_CACHE)
}

set.suppressEchoCmdInFuncOption <- function(value) {
	if (!is.logical(value)) stop ("logical argument required")
	res <- get("suppressEchoCmdInFunc", envir = .GRASS_CACHE)
	assign("suppressEchoCmdInFunc", value, envir = .GRASS_CACHE)
	res
}

get.suppressEchoCmdInFuncOption <- function() {
	get("suppressEchoCmdInFunc", envir = .GRASS_CACHE)
}

