
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C) for this R-port:
#   1999 - 2012 Diethelm Wuertz, Zurich, <wuertz@itp.phys.ethz.ch>
#   2009 - 2012 Rmetrics Association, Zurich, www.rmetrics.org


################################################################################
# FUNCTION:             DESCRIPTION:
#  .dQuote.ascii         Creates a temporary function to quote a string
#  .xls2csv              Converts a xls file into acsv file
#  .read.xls             reads low level a xls files
################################################################################


# Package: gdata
# Title: Various R programming tools for data manipulation
# Description: Various R programming tools for data manipulation
# Version: 2.4.2
# Date: 2008-05-12
# Author: Gregory R. Warnes and Gregor Gorjanc.
# Maintainer: Gregory Warnes <gregory_warnes@urmc.rochester.edu>
# License: GPL-2


.dQuote.ascii <-
function(x)
{
    # Borrowed from contributed Package 'gdata'

    # FUNCTION:

    # Creating a temporary function to quote the string
    paste('"', x,'"', sep = '')
}


#-------------------------------------------------------------------------------


.xls2csv <-
    function(xls, sheet = 1, verbose = FALSE, ..., perl = "perl")
{
    # Borrowed from contributed Package 'gdata'

    # FUNCTION:

    # Directories:
    package.dir <- path.package('fImport')
    perl.dir <- file.path(package.dir,'perl')

    # Files:
    tf <- NULL
    if (substring(xls, 1, 7) == "http://") {
        tf <- paste(tempfile(), "xls", sep = ".")
            if(verbose)
                cat("Downloading",
                    .dQuote.ascii(xls), " to ",
                    .dQuote.ascii(tf), "...\n")
            else
          cat("Downloading...\n")
        download.file(xls, tf, mode = "wb")
        cat("Done.\n")
        xls <- tf
    }

    if(file.access(xls, 4)!=0)
        stop("Unable to read xls file '", xls, "'." )

    .xls2csv <- file.path(perl.dir,'xls2csv.pl')
    csv <- paste(tempfile(), "csv", sep = ".")

    # Execution command
    cmd <- paste(perl, .xls2csv, .dQuote.ascii(xls), .dQuote.ascii(csv),
        sheet, sep=" ")

    if(verbose)
      {
        cat("\n")
        cat("Converting xls file\n")
        cat("   ", .dQuote.ascii(xls), "\n")
        cat("to csv file \n")
        cat("   ", .dQuote.ascii(csv), "\n")
        cat("... \n\n")
      }
      else
        cat("Converting xls file to csv file... ")

    # Do the translation
    if(verbose)  cat("Executing ", cmd, "... \n\n")
    results <- system(cmd, intern=!verbose)
    if (verbose) cat("Done.\n\n")


    if(file.access(csv, 4)!=0)
        stop("Unable to read translated csv file '", csv, "'." )
    cat("Done.\n")


    # Prepare for cleanup now, in case of error reading file
    file(csv)
}


# ------------------------------------------------------------------------------


.read.xls <-
    function(xls, sheet = 1, verbose = FALSE, pattern, ..., perl = "perl")
{
    # Borrowed from contributed Package 'gdata'

    # FUNCTION:

    # Connection:
    con <- tfn <- NULL
    on.exit({
        if (inherits(con, "connection") && isOpen(con)) close(con)
        if (file.exists(tfn)) file.remove(tfn)})

    # Expand file path, translating ~ to user's home directory, etc.
    xls <- path.expand(xls)


    # Translate from xls to csv format (returns csv file name)
    con <- .xls2csv(xls, sheet, verbose=verbose, ..., perl = perl)

    # Load the csv file
    open(con)
    tfn <- summary(con)$description
    if (missing(pattern)) {
        if(verbose)
            cat("Reading csv file ", .dQuote.ascii(tfn), "...\n")
        else
            cat("Reading csv file... ")
        retval <- read.csv(con, ...)
        cat("Done.\n")
    } else {
        cat("Searching for lines containing pattern ", pattern, "... ")
        idx <- grep(pattern, readLines(con))
        if (length(idx) == 0) {
            warning("pattern not found")
            return(NULL)
        }
        cat("Done.\n")

        seek(con, 0)

        if(verbose)
            cat("Reading csv file ", .dQuote.ascii(tfn), "...\n")
        else
            cat("Reading csv file... ")
        retval <- read.csv(con, skip = idx[1]-1, ...)
        cat("Done.\n")
    }

    # Return Value:
    retval
}


################################################################################
