#  zzz.R
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


.SYBILenv <- new.env()

.onLoad <- function(lib, pkg) {

    # -------------------------------------------------------------- #
    # stacks and queues in sybil

    .SYBILenv$sybilStack <- vector(mode = "list", length = 0)
    

    # -------------------------------------------------------------- #
    # settings in sybil

    .SYBILenv$settings <- list(
        SOLVER           = "glpkAPI",
        METHOD           = "simplex",
        TOLERANCE        = 1E-6,
        MAXIMUM          = 1000,
        ALGORITHM        = "fba",
        OPT_DIRECTION    = "max",
        USE_NAMES        = FALSE,
        PATH_TO_MODEL    = ".",
        SOLVER_CTRL_PARM = as.data.frame(NA)
    )


    # -------------------------------------------------------------- #
    # available solvers and methods

    # The first method in each vector is the default one, the default
    # alternate method is NA.

    # solvers
    .SYBILenv$solvers <- c("glpkAPI", "clpAPI", "lpSolveAPI", "cplexAPI")

    # methods
    .SYBILenv$solverMethods <- list(
        glpkAPI    = c("simplex", "interior", "exact", "mip"),

        clpAPI     = c("general_solve", "inidual", "iniprimal", "inibarrier",
                       "inibarriernoc", "idiot", "dual", "primal"),

        lpSolveAPI = c("lp_solve"),

        cplexAPI   = c("lpopt", "primopt", "dualopt", "baropt", "hybbaropt",
                       "hybnetopt", "siftopt", "mipopt", "qpopt")
    )

    # default parameters
    .SYBILenv$solverCtrlParm <- list(

        glpkAPI    = list(simplex  = as.data.frame(NA),
                          interior = as.data.frame(NA),
                          exact    = as.data.frame(NA),
                          mip      = as.data.frame(NA)
        ),

        clpAPI     = list(general_solve = as.data.frame(NA),
                          inidual       = as.data.frame(NA),
                          iniprimal     = as.data.frame(NA),
                          inibarrier    = as.data.frame(NA),
                          inibarriernoc = as.data.frame(NA),
                          dual          = as.data.frame(NA),
                          primal        = as.data.frame(NA)
        ),

        lpSolveAPI = list(lp_solve = as.data.frame(NA)),

        cplexAPI   = list(lpopt     = as.data.frame(NA),
                          primopt   = as.data.frame(NA),
                          #primopt   = data.frame(
                          #   #CPX_PARAM_REDUCE    = as.integer(1),
                          #    CPX_PARAM_DATACHECK = as.integer(1)
                          #),
                          dualopt   = as.data.frame(NA),
                          baropt    = as.data.frame(NA),
                          hybbaropt = as.data.frame(NA),
                          hybnetopt = as.data.frame(NA),
                          siftopt   = as.data.frame(NA),
                          mipopt    = as.data.frame(NA),
                          qpopt     = as.data.frame(NA)
        )
    )


    # -------------------------------------------------------------- #
    # algorithms simulation genetic perturbations

    .SYBILenv$algorithm[["pert"]] <- c("fba", "mtf", "moma", "lmoma", "room")


    # -------------------------------------------------------------- #
    # methods problem types

    .SYBILenv$ptype[["mip"]] <- list(
        glpkAPI     = c("mip"),
        lpSolveAPI  = c("lp_solve"),
        cplexAPI    = c("mipopt")
    )

    .SYBILenv$ptype[["qp"]] <- list(
        cplexAPI     = c("qpopt", "baropt")
    )

}


