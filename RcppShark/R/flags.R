## RcppShark.package.skeleton.R: makes a skeleton for a package that wants to use RcppShark
##
## Copyright (C)  2014  Qiang Kou
##
## This file was part of RcppMLPACK.
##
## RcppMLPACK is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppMLPACK distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppMLPACK.  If not, see <http://www.gnu.org/licenses/>.

RcppShark.system.file <- function(...){
    tools::file_path_as_absolute( base::system.file( ..., package = "RcppShark" ) )
}

staticLinking <- function() {
    ! grepl( "^linux", R.version$os )
}

RcppSharkLdPath <- function() {
    if (nzchar(.Platform$r_arch)) {	## eg amd64, ia64, mips
        path <- RcppShark.system.file("lib",.Platform$r_arch)
    } else {
        path <- RcppShark.system.file("lib")
    }
    path
}

RcppSharkLdFlags <- function(static=staticLinking()) {
    RcppSharkdir <- RcppSharkLdPath()
    if (static) {                               # static is default on Windows and OS X
        flags <- paste(RcppSharkdir, "/libRcppShark.a", sep="")
    } else {					# else for dynamic linking
        flags <- paste("-L", RcppSharkdir, " -lRcppShark", sep="") # baseline setting
        if ((.Platform$OS.type == "unix") &&    # on Linux, we can use rpath to encode path
            (length(grep("^linux",R.version$os)))) {
            flags <- paste(flags, " -Wl,-rpath,", RcppSharkdir, sep="")
        }
    }
    invisible(flags)
}

CxxFlags <- function() {
    path <- RcppShark.system.file("include")
    paste("-I", path, sep="")
}

LdFlags <- function(static=staticLinking()) {
    cat(RcppSharkLdFlags(static=static))
}
