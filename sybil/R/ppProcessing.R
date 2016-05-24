#  ppProcessing.R
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
# Function: ppProcessing
#
# do some kind of post or pre processing
#
#


.ppProcessing <- function(lpprob, ppCmd, loopvar = NA) {

    #----------------------------------------------------------------------#
    # run the pre/post processing command
    #----------------------------------------------------------------------#

    .runPPcommand <- function(ppCommand) {

        pAnaR <- NA

        if(!is(ppCommand, "character")) {
            warning("Each command must be passed as character vector.")
            ppC <- NA
        }

        if (nchar(ppCommand[1]) < 1) {
            warning("The first element of 'ppC' must contain a function name.")
            ppC <- NA
        }
        else {
            ppC <- ppCommand
        }
#         else {
#             if(!isTRUE(is.function(eval(parse(text = ppC[1]))))) {
#                 warning(paste("The first element of 'ppC' must contain a valid",
#                               "function name."
#                         )
#                 )
#                 ppC <- NA
#             }
#         }

        if (!is.na(ppC[1])) {
            narg <- length(ppC)
            if (narg == 1) {
                ppC <- paste(ppC, "(lpprob)", sep = "")
            }
            else {
                if (!is.na(loopvar)) {
                    ppC <- gsub("LOOP_VAR", loopvar, ppC, fixed = TRUE)
                }
                ppC <- gsub("LP_PROB", "lpprob", ppC, fixed = TRUE)
                ppC <- paste(ppC[1], "(",
                            paste(ppC[-1], collapse = ", "),
                            ")", sep = ""
                       )
            }

            # run the command
            pAnaR <- tryCatch(eval(parse(text = ppC)), error = function(x) x)
            if (is(pAnaR, "simpleError")) {
                pAnaR <- sybilError(pAnaR)
            }

            ppC <- gsub("lpprob", "LP_PROB", ppC, fixed = TRUE)
            if (!is.na(loopvar)) {
                ppC <- gsub(loopvar, "LOOP_VAR", ppC, fixed = TRUE)
            }

        }

        return(list(com = ppC, res = pAnaR))
    }


#------------------------------------------------------------------------------#


    # cmdList: list of commands
    # resList: list of results for each command
    # instance of class ppProc, if successfull, otherwise NA

    cmdList <- list(NA)
    resList <- list(NA)
    pAna    <- NA

    if (all(!is.na(ppCmd))) {
        if (!is(ppCmd, "list")) {
            warning("Argument 'ppCmd' must be a list.")
        }
        else {
            tmp_res <- NULL
            nCmd    <- length(ppCmd)
            resList <- vector(mode = "list", length = nCmd)
            cmdList <- vector(mode = "list", length = nCmd)
            for (i in seq(along = ppCmd)) {
                tmp_res <- .runPPcommand(ppCmd[[i]])
                resList[[i]] <- tmp_res[["res"]]
                cmdList[[i]] <- tmp_res[["com"]]
                tmp_res <- NULL
            }
        }

        pAna <- ppProc(cmdList)
        pa(pAna) <- resList

    }

    return(pAna)
}


