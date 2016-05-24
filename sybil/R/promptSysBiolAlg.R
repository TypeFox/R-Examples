#  promptSysBiolAlg.R
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


################################################
# Function: promptSysBiolAlg
#
# 
# 

promptSysBiolAlg <- function(algorithm,
                             prefix = "sysBiolAlg",
                             sep = "_",
                             suffix = "R",
                             fpath = ".",
                             ...) {

    stopifnot(is(algorithm, "character"))

    on.exit( closeAllConnections() )
    
    # classname
    cname <- paste(prefix, algorithm, sep = sep)
    
    # filename
    fname <- paste(paste(cname, "Class", sep = ""), suffix, sep = ".")

    # path to output file
    sbfile <- file.path(fpath, fname)

    sbfh <- try(file(description = sbfile, open = "wt", ...), silent = TRUE)

    if (is(sbfh, "try-error")) {
        stop("can not write to file ", sQuote(sbfh))
    }

    #--------------------------------------------------------------------------#
    # write file
    #--------------------------------------------------------------------------#

    cat("#",rep("-", 78), "#\n",
        "# definition of class ", cname, "\n",
        "#",rep("-", 78), "#\n",
        "\n",
        "setClass(Class = \"", cname, "\",\n",
        "         contains = \"sysBiolAlg\"\n",
        ")",
        "\n\n",
        "#",rep("-", 78), "#\n",
        "# default constructor\n",
        "#",rep("-", 78), "#\n",
        "\n",
        "# contructor for class ", cname, "\n",
        sep = "", file = sbfh)

    cat("setMethod(f = \"initialize\",\n",
        "          signature = \"", cname, "\",\n",
        "          definition = function(.Object,\n",
        "                                model,\n",
        "## PLACE FURTHER ARGUMENTS TO CONSTRUCTOR HERE                 ##\n",
        "## ARGUMENT '...' SHOULD BE LEFT FOR ARGUMENTS                 ##\n",
        "## 'solver', 'method', 'solverParm', 'termOut' AND 'retAlgPar' ##\n",
        "                                useNames = SYBIL_SETTINGS(\"USE_NAMES\"),\n",
        "                                cnames = NULL,\n",
        "                                rnames = NULL,\n",
        "                                pname = NULL,\n",
        "                                scaling = NULL,\n",
        "                                writeProbToFileName = NULL, ...) {\n\n",
        sep = "", file = sbfh, append = TRUE)



        
    cat("              if ( ! missing(model) ) {\n",
        "\n",
        "                  stopifnot(is(model, \"modelorg\"),\n",
        "## PLACE FURTHER ARGUMENT TESTS HERE                           ##\n",
        "                            )\n\n",
        sep = "", file = sbfh, append = TRUE)


    cat("## GENERATE ROW AND COLUMN NAMES FOR                           ##\n",
        "## THE PROBLEM OBJECT (OPTIONALLY)                             ##\n",
        sep = "", file = sbfh, append = TRUE)

    cat("                 if (isTRUE(useNames)) {\n",
        "                     # you maybe want to do something more fancy here:\n",
        "                     # e.g. use .makeLPcompatible() for a default value\n",
        "                     colNames <- cnames\n",
        "                     rowNames <- rnames\n",
        "                     probName <- pname\n",
        "                 }\n",
        "                 else {\n",
        "                     colNames <- NULL\n",
        "                     rowNames <- NULL\n",
        "                     probName <- NULL\n",
        "                 }\n\n",
        sep = "", file = sbfh, append = TRUE)

    cat("## GENERATE DATA STRUCTURES TO USE FOR                         ##\n",
        "## BUILDING THE PROBLEM OBJECT VIA THE DEFAULT                 ##\n",
        "## CONSTRUCTOR FOR CLASS sysBiolAlg                            ##\n",
        sep = "", file = sbfh, append = TRUE)


    cat("\n", sep = "", file = sbfh, append = TRUE)

    cat("                  # use the default constructor for class sysBiolAlg\n",
        "                  .Object <- callNextMethod(.Object,\n",
        "                                            sbalg      = \"", algorithm, "\",\n",
        "                                            pType      = \n",
        "                                            scaling    = scaling,\n",
        "                                            fi         = \n",
        "                                            nCols      = \n",
        "                                            nRows      = \n",
        "                                            mat        = \n",
        "                                            ub         = \n",
        "                                            lb         = \n",
        "                                            obj        = \n",
        "                                            rlb        = \n",
        "                                            rtype      = \n",
        "                                            lpdir      = \n",
        "                                            rub        = \n",
        "                                            ctype      = \n",
        "                                            cnames     = colNames,\n",
        "                                            rnames     = rowNames,\n",
        "                                            pnames     = probName,\n",
        "                                            algPar     = list(),\n",
        "                                            ...)\n\n",
        sep = "", file = sbfh, append = TRUE)

    cat("                  # write problem to lp file\n",
        "                  if (!is.null(writeToFileName)) {\n",
        "                      writeProb(problem(.Object),\n",
        "                                fname = as.character(writeProbToFileName))\n",
        "                  }\n\n",
        sep = "", file = sbfh, append = TRUE)

    cat("                  validObject(.Object)\n",
        sep = "", file = sbfh, append = TRUE)

    cat("              }\n", sep = "", file = sbfh, append = TRUE)
    cat("              return(.Object)\n", sep = "", file = sbfh, append = TRUE)
    cat("          }\n", sep = "", file = sbfh, append = TRUE)
    cat(")\n", sep = "", file = sbfh, append = TRUE)

    #--------------------------------------------------------------------------#
    # end
    #--------------------------------------------------------------------------#

    if ( (is(sbfh, "file")) && (isOpen(sbfh)) ) {
        close(sbfh)
    }

    message("created file ", sQuote(sbfile))
    
    return(invisible(NULL))

}
