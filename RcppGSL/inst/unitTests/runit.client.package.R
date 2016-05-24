#!/usr/bin/r -t
#
# Copyright (C) 2010         Romain Francois and Dirk Eddelbuettel
# Copyright (C) 2014 - 2015  Dirk Eddelbuettel
#
# This file is part of RcppGSL.
#
# RcppGSL is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# RcppGSL is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppGSL.  If not, see <http://www.gnu.org/licenses/>.

test.client.package <- function(){
    pkg <- "RcppGSLExample"
    td <- tempfile()
    cwd <- getwd()
    dir.create(td)
    file.copy(system.file("examples", pkg, package = "RcppGSL" ) , td, recursive = TRUE)
    Sys.chmod(file.path(td, "RcppGSLExample", "configure"), "0755")
    setwd(td)
    on.exit( {
        setwd(cwd)
        unlink(td, recursive = TRUE)
    } )
    R <- shQuote(file.path(R.home(component = "bin" ), "R"))
    cmd <- paste(R , "CMD build", pkg)
    system(cmd)
    dir.create("templib")
    install.packages(paste(pkg, "_0.0.3.tar.gz", sep = ""), "templib", repos = NULL, type = "source")
    require(Rcpp)
    require(pkg, "templib", character.only = TRUE)
    m <- matrix(1:16, nc = 4)
    res <- colNorm(m)
    val <- apply(m, 2, function(x) sqrt(sum(x^2)))
    unlink("templib", recursive = TRUE)
    checkEquals(res, val, msg = "colNorm in client package")
}
