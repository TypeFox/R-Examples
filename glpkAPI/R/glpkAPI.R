#------------------------------------------------------------------------------#
#                             R interface to GLPK                              #
#------------------------------------------------------------------------------#

#  glpkAPI.R
#  R interface to GLPK.
#
#  Copyright (C) 2011-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of glpkAPI.
#
#  GlpkAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GlpkAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with glpkAPI  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                              the interface                                   #
#------------------------------------------------------------------------------#


delProbGLPK <- function(lp) {

    invisible(
        .Call("delProb", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )

}


#------------------------------------------------------------------------------#

eraseProbGLPK <- function(lp) {

    invisible(
        .Call("eraseProb", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )

}


#------------------------------------------------------------------------------#

copyProbGLPK <- function(lp, clp, name = GLP_OFF) {

    invisible(
        .Call("copyProb", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              glpkPointer(clp),
              as.integer(name)
        )
    )

}


#------------------------------------------------------------------------------#

initProbGLPK <- function(ptrtype = "glpk_prob") {

    lp <- .Call("initProb", PACKAGE = "glpkAPI",
                as.character(ptrtype)
          )

    lpP <- glpk_Pointer(lp)

    return(lpP)
}


#------------------------------------------------------------------------------#

setProbNameGLPK <- function(lp, pname = NULL) {

    if (is.null(pname)) {
        Cpname <- as.null(pname)
    }
    else {
        Cpname <- as.character(pname)
    }

    invisible(
        .Call("setProbName", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              Cpname
        )
    )

}


#------------------------------------------------------------------------------#

getProbNameGLPK <- function(lp) {

    pname <- .Call("getProbName", PACKAGE = "glpkAPI",
                   glpkPointer(lp)
             )

    return(pname)

}


#------------------------------------------------------------------------------#

setObjNameGLPK <- function(lp, oname = NULL) {

    if (is.null(oname)) {
        Coname <- as.null(oname)
    }
    else {
        Coname <- as.character(oname)
    }

    invisible(
        .Call("setObjName", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              Coname
        )
    )

}


#------------------------------------------------------------------------------#

getObjNameGLPK <- function(lp) {

    oname <- .Call("getObjName", PACKAGE = "glpkAPI",
                   glpkPointer(lp)
             )

    return(oname)

}


#------------------------------------------------------------------------------#

createIndexGLPK <- function(lp) {

    invisible(
        .Call("createIndex", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )

}


#------------------------------------------------------------------------------#

deleteIndexGLPK <- function(lp) {

    invisible(
        .Call("deleteIndex", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )

}


#------------------------------------------------------------------------------#

setDefaultSmpParmGLPK <- function() {

    invisible(.Call("setDefaultSmpParm", PACKAGE = "glpkAPI"))

}


#------------------------------------------------------------------------------#

setDefaultIptParmGLPK <- function() {

    invisible(.Call("setDefaultIptParm", PACKAGE = "glpkAPI"))

}


#------------------------------------------------------------------------------#

setDefaultMIPParmGLPK <- function() {

    invisible(.Call("setDefaultMIPParm", PACKAGE = "glpkAPI"))

}


#------------------------------------------------------------------------------#

setSimplexParmGLPK <- function(parm, val) {

    if (!identical(length(parm), length(val))) {
        stop("Arguments 'parm' and 'val' must have the same length!")
    }

    indi <- which(parm > 100 & parm < 200)
    indd <- which(parm > 200 & parm < 300)

    npari <- length(indi)
    npard <- length(indd)

    if (npari == 0) {
        parmi <- as.null(parm)
        vali  <- as.null(val)
    }
    else {
        parmi <- as.integer(parm[indi])
        vali  <- as.integer(val[indi])
    }

    if (npard == 0) {
        parmd <- as.null(parm)
        vald  <- as.null(val)
    }
    else {
        parmd <- as.integer(parm[indd])
        vald  <- as.numeric(val[indd])
    }

    invisible(
        .Call("setSimplexParm", PACKAGE = "glpkAPI",
              as.integer(npari),
              parmi,
              vali,
              as.integer(npard),
              parmd,
              vald
        )
    )

}


#------------------------------------------------------------------------------#

setInteriorParmGLPK <- function(parm, val) {

    if (!identical(length(parm), length(val))) {
        stop("Arguments 'parm' and 'val' must have the same length!")
    }

    nparm <- length(parm)

    invisible(
        .Call("setInteriorParm", PACKAGE = "glpkAPI",
              as.integer(nparm),
              as.integer(parm),
              as.integer(val)
        )
    )

}


#------------------------------------------------------------------------------#

setMIPParmGLPK <- function(parm, val) {

    if (!identical(length(parm), length(val))) {
        stop("Arguments 'parm' and 'val' must have the same length!")
    }

    indi <- which( (parm > 100 & parm < 200) | (parm > 600 & parm < 700) )
    indd <- which(parm > 700 & parm < 800)

    npari <- length(indi)
    npard <- length(indd)

    if (npari == 0) {
        parmi <- as.null(parm)
        vali  <- as.null(val)
    }
    else {
        parmi <- as.integer(parm[indi])
        vali  <- as.integer(val[indi])
    }

    if (npard == 0) {
        parmd <- as.null(parm)
        vald  <- as.null(val)
    }
    else {
        parmd <- as.integer(parm[indd])
        vald  <- as.numeric(val[indd])
    }

    invisible(
        .Call("setMIPParm", PACKAGE = "glpkAPI",
              as.integer(npari),
              parmi,
              vali,
              as.integer(npard),
              parmd,
              vald
        )
    )

}


#------------------------------------------------------------------------------#

getSimplexParmGLPK <- function() {

    parmS <- .Call("getSimplexParm", PACKAGE = "glpkAPI")

    return(parmS)

}


#------------------------------------------------------------------------------#

getInteriorParmGLPK <- function() {

    parmI <- .Call("getInteriorParm", PACKAGE = "glpkAPI")

    return(parmI)

}


#------------------------------------------------------------------------------#

getMIPParmGLPK <- function() {

    parmM <- .Call("getMIPParm", PACKAGE = "glpkAPI")

    return(parmM)

}


#------------------------------------------------------------------------------#

setObjDirGLPK <- function(lp, lpdir) {

    invisible(
        .Call("setObjDir", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(lpdir)
        )
    )

}


#------------------------------------------------------------------------------#

getObjDirGLPK <- function(lp) {

    lpdir <- .Call("getObjDir", PACKAGE = "glpkAPI",
                   glpkPointer(lp)
                  )
    return(lpdir)

}


#------------------------------------------------------------------------------#

addRowsGLPK <- function(lp, nrows) {

    frow <- .Call("addRows", PACKAGE = "glpkAPI",
                  glpkPointer(lp),
                  as.integer(nrows)
             )
    return(frow)

}


#------------------------------------------------------------------------------#

setRowNameGLPK <- function(lp, i, rname = NULL) {

    if (is.null(rname)) {
        Crname <- as.null(rname)
    }
    else {
        Crname <- as.character(rname)
    }

    invisible(
        .Call("setRowName", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i),
              Crname
        )
    )

}


#------------------------------------------------------------------------------#

setRowsNamesGLPK <- function(lp, i, rnames = NULL) {

    if (is.null(rnames)) {
        Crnames <- as.null(rnames)
    }
    else {
        Crnames <- as.character(rnames)
    }

    invisible(
        .Call("setRowsNames", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i),
              Crnames
        )
    )

}


#------------------------------------------------------------------------------#

getRowNameGLPK <- function(lp, i) {

    rname <- .Call("getRowName", PACKAGE = "glpkAPI",
                   glpkPointer(lp),
                   as.integer(i)
             )

    return(rname)

}


#------------------------------------------------------------------------------#

findRowGLPK <- function(lp, rname) {

    rind <- .Call("findRow", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(rname)
        )

    return(rind)

}


#------------------------------------------------------------------------------#

addColsGLPK <- function(lp, ncols) {

    fcol <- .Call("addCols", PACKAGE = "glpkAPI",
                  glpkPointer(lp),
                  as.integer(ncols)
             )

    return(fcol)

}


#------------------------------------------------------------------------------#

setColNameGLPK <- function(lp, j, cname = NULL) {

    if (is.null(cname)) {
        Ccname <- as.null(cname)
    }
    else {
        Ccname <- as.character(cname)
    }

    invisible(
        .Call("setColName", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              Ccname
        )
    )

}


#------------------------------------------------------------------------------#

setColsNamesGLPK <- function(lp, j, cnames = NULL) {

    if (is.null(cnames)) {
        Ccnames <- as.null(cnames)
    }
    else {
        Ccnames <- as.character(cnames)
    }

    invisible(
        .Call("setColsNames", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              Ccnames
        )
    )

}


#------------------------------------------------------------------------------#

getColNameGLPK <- function(lp, j) {

    cname <- .Call("getColName", PACKAGE = "glpkAPI",
                   glpkPointer(lp),
                   as.integer(j)
             )

    return(cname)

}


#------------------------------------------------------------------------------#

findColGLPK <- function(lp, cname) {

    cind <- .Call("findCol", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(cname)
        )

    return(cind)

}


#------------------------------------------------------------------------------#

getNumRowsGLPK <- function(lp) {

    nrows <- .Call("getNumRows", PACKAGE = "glpkAPI",
                   glpkPointer(lp)
             )
    return(nrows)

}


#------------------------------------------------------------------------------#

getNumColsGLPK <- function(lp) {

    ncols <- .Call("getNumCols", PACKAGE = "glpkAPI",
                   glpkPointer(lp)
                  )
    return(ncols)

}


#------------------------------------------------------------------------------#

setColsBndsGLPK <- function(lp, j, lb, ub, type = NULL) {

    if (is.null(type)) {
        Ctype <- as.null(type)
    }
    else {
        Ctype <- as.integer(type)
    }

    invisible(
        .Call("setColsBnds", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              Ctype,
              as.numeric(lb),
              as.numeric(ub)
        )
    )

}



#------------------------------------------------------------------------------#

setColsBndsObjCoefsGLPK <- function(lp, j, lb, ub, obj_coef, type = NULL) {

    if (is.null(type)) {
        Ctype <- as.null(type)
    }
    else {
        Ctype <- as.integer(type)
    }

    invisible(
        .Call("setColsBndsObjCoefs", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              Ctype,
              as.numeric(lb),
              as.numeric(ub),
              as.numeric(obj_coef)
        )
    )

}

#------------------------------------------------------------------------------#

setColBndGLPK <- function(lp, j, type, lb, ub) {

    invisible(
        .Call("setColBnd", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              as.integer(type),
              as.numeric(lb),
              as.numeric(ub)
        )
    )
}


#------------------------------------------------------------------------------#

getColsLowBndsGLPK <- function(lp, j) {

    lowbnd <- .Call("getColsLowBnds", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )

    return(lowbnd)

}


#------------------------------------------------------------------------------#

getColLowBndGLPK <- function(lp, j) {

    lowbnd <- .Call("getColLowBnd", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )
    return(lowbnd)

}


#------------------------------------------------------------------------------#

getColsUppBndsGLPK <- function(lp, j) {

    uppbnd <- .Call("getColsUppBnds", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )

    return(uppbnd)

}


#------------------------------------------------------------------------------#

getColUppBndGLPK <- function(lp, j) {

    uppbnd <- .Call("getColUppBnd", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )

    return(uppbnd)

}


#------------------------------------------------------------------------------#

setColKindGLPK <- function(lp, j, kind) {

    invisible(
        .Call("setColKind", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              as.integer(kind)
        )
    )

}


#------------------------------------------------------------------------------#

setColsKindGLPK <- function(lp, j, kind) {

    invisible(
        .Call("setColsKind", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              as.integer(kind)
        )
    )

}


#------------------------------------------------------------------------------#

getColKindGLPK <- function(lp, j) {

    kind <- .Call("getColKind", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )

    return(kind)

}


#------------------------------------------------------------------------------#

getColsKindGLPK <- function(lp, j) {

    kind <- .Call("getColsKind", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )

    return(kind)

}


#------------------------------------------------------------------------------#

getNumIntGLPK <- function(lp) {

    num <- .Call("getNumInt", PACKAGE = "glpkAPI",
                 glpkPointer(lp)
            )

    return(num)

}


#------------------------------------------------------------------------------#

getNumBinGLPK <- function(lp) {

    num <- .Call("getNumBin", PACKAGE = "glpkAPI",
                 glpkPointer(lp)
            )

    return(num)

}


#------------------------------------------------------------------------------#

setRowsBndsGLPK <- function(lp, i, lb, ub, type = NULL) {

    if (is.null(type)) {
        Ctype <- as.null(type)
    }
    else {
        Ctype <- as.integer(type)
    }

    invisible(
        .Call("setRowsBnds", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i),
              Ctype,
              as.numeric(lb),
              as.numeric(ub)
        )
    )

}


#------------------------------------------------------------------------------#

setRhsZeroGLPK <- function(lp) {

    invisible(
        .Call("setRhsZero", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )

}


#------------------------------------------------------------------------------#

setRowBndGLPK <- function(lp, i, type, lb, ub) {

    invisible(
        .Call("setRowBnd", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i),
              as.integer(type),
              as.numeric(lb),
              as.numeric(ub)
        )
    )
}


#------------------------------------------------------------------------------#

getRowsLowBndsGLPK <- function(lp, i) {

    lowbnd <- .Call("getRowsLowBnds", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i)
        )

    return(lowbnd)

}


#------------------------------------------------------------------------------#

getRowLowBndGLPK <- function(lp, i) {

    lowbnd <- .Call("getRowLowBnd", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i)
        )
    return(lowbnd)

}


#------------------------------------------------------------------------------#

getRowsUppBndsGLPK <- function(lp, i) {

    uppbnd <- .Call("getRowsUppBnds", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i)
        )

    return(uppbnd)

}


#------------------------------------------------------------------------------#

getRowUppBndGLPK <- function(lp, i) {

    uppbnd <- .Call("getRowUppBnd", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i)
        )
    return(uppbnd)

}


#------------------------------------------------------------------------------#

getRowTypeGLPK <- function(lp, i) {

    type <- .Call("getRowType", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i)
        )

    return(type)

}


#------------------------------------------------------------------------------#

getRowsTypesGLPK <- function(lp, i) {

    type <- .Call("getRowsTypes", PACKAGE = "glpkAPI",
                  glpkPointer(lp),
                  as.integer(i)
        )

    return(type)

}


#------------------------------------------------------------------------------#

getColTypeGLPK <- function(lp, j) {

    type <- .Call("getColType", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )

    return(type)

}


#------------------------------------------------------------------------------#

setObjCoefsGLPK <- function(lp, j, obj_coef) {

    invisible(
        .Call("setObjCoefs", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              as.numeric(obj_coef)
        )
    )

}


#------------------------------------------------------------------------------#

setObjCoefGLPK <- function(lp, j, obj_coef) {

    invisible(
        .Call("setObjCoef", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              as.numeric(obj_coef)
        )
    )

}


#------------------------------------------------------------------------------#

getObjCoefsGLPK <- function(lp, j) {

    obj_coef <- .Call("getObjCoefs", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )

    return(obj_coef)

}


#------------------------------------------------------------------------------#

getObjCoefGLPK <- function(lp, j) {

    obj_coef <- .Call("getObjCoef", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )
    return(obj_coef)

}


#------------------------------------------------------------------------------#

loadMatrixGLPK <- function(lp, ne, ia, ja, ra) {

    invisible(
        .Call("loadMatrix", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(ne),
              as.integer(ia),
              as.integer(ja),
              as.numeric(ra)
        )
    )

}


#------------------------------------------------------------------------------#

checkDupGLPK <- function(m, n, ne, ia, ja) {

    dup <- .Call("checkDup", PACKAGE = "glpkAPI",
                 as.integer(m),
                 as.integer(n),
                 as.integer(ne),
                 as.integer(ia),
                 as.integer(ja)
    )

    return(dup)

}


#------------------------------------------------------------------------------#

sortMatrixGLPK <- function(lp) {

    invisible(
        .Call("sortMatrix", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )

}


#------------------------------------------------------------------------------#

delRowsGLPK <- function(lp, nrows, i) {

    if (i[1] != 0) {
        Ci <- as.integer(append(i, 0, 0))
    }
    else {
        Ci <- as.integer(i)
    }

    invisible(
        .Call("delRows", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(nrows),
              Ci
        )
    )

}


#------------------------------------------------------------------------------#

delColsGLPK <- function(lp, ncols, j) {


    if (j[1] != 0) {
        Cj <- as.integer(append(j, 0, 0))
    }
    else {
        Cj <- as.integer(j)
   }

    invisible(
        .Call("delCols", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(ncols),
              Cj
        )
    )

}


#------------------------------------------------------------------------------#

setRiiGLPK <- function(lp, i, rii) {

    invisible(
        .Call("setRii", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i),
              as.numeric(rii)
        )
    )

}


#------------------------------------------------------------------------------#

setSjjGLPK <- function(lp, j, sjj) {

    invisible(
        .Call("setSjj", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              as.numeric(sjj)
        )
    )

}


#------------------------------------------------------------------------------#

getRiiGLPK <- function(lp, i) {

    rii <- .Call("getRii", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i)
        )
    return(rii)

}


#------------------------------------------------------------------------------#

getSjjGLPK <- function(lp, j) {

    sjj <- .Call("getSjj", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j)
        )
    return(sjj)

}


#------------------------------------------------------------------------------#

scaleProbGLPK <- function(lp, opt) {

    invisible(
        .Call("scaleProb", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(opt)
        )
    )

}


#------------------------------------------------------------------------------#

unscaleProbGLPK <- function(lp) {

    invisible(
        .Call("unscaleProb", PACKAGE = "glpkAPI",
            glpkPointer(lp)
        )
    )

}


#------------------------------------------------------------------------------#

setRowStatGLPK <- function(lp, i, stat) {

    invisible(
        .Call("setRowStat", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i),
              as.integer(stat)
        )
    )

}


#------------------------------------------------------------------------------#

setColStatGLPK <- function(lp, j, stat) {

    invisible(
        .Call("setColStat", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              as.integer(stat)
        )
    )

}


#------------------------------------------------------------------------------#

stdBasisGLPK <- function(lp) {

    invisible(
        .Call("stdBasis", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )

}


#------------------------------------------------------------------------------#

advBasisGLPK <- function(lp) {

    invisible(
        .Call("advBasis", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )
}


#------------------------------------------------------------------------------#

cpxBasisGLPK <- function(lp) {

    invisible(
        .Call("cpxBasis", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )
    )
}


#------------------------------------------------------------------------------#

warmUpGLPK <- function(lp) {

    wup <- .Call("warmUp", PACKAGE = "glpkAPI",
                 glpkPointer(lp)
    )

    return(wup)
}


# ------------------------------------------------------------------------------

termOutGLPK <- function(flag) {

    fl <- .Call("termOut", PACKAGE = "glpkAPI",
                as.integer(flag)
    )

    return(fl)
}


#------------------------------------------------------------------------------#

solveSimplexGLPK <- function(lp) {

    ret <- .Call("solveSimplex", PACKAGE = "glpkAPI",
              glpkPointer(lp)
           )

    return(ret)
}


#------------------------------------------------------------------------------#

solveSimplexExactGLPK <- function(lp) {

    ret <- .Call("solveSimplexExact", PACKAGE = "glpkAPI",
              glpkPointer(lp)
           )

    return(ret)
}


#------------------------------------------------------------------------------#

getObjValGLPK <- function(lp) {

    obj <- .Call("getObjVal", PACKAGE = "glpkAPI",
              glpkPointer(lp)
           )

    return(obj)
}


#------------------------------------------------------------------------------#

getSolStatGLPK <- function(lp) {

    stat <- .Call("getSolStat", PACKAGE = "glpkAPI",
              glpkPointer(lp)
            )

    return(stat)
}


#------------------------------------------------------------------------------#

getColsPrimGLPK <- function(lp) {

    col_prim <- .Call("getColsPrim", PACKAGE = "glpkAPI",
                      glpkPointer(lp)
                )

    return(col_prim)
}


#------------------------------------------------------------------------------#

getColPrimGLPK <- function(lp, j) {

    col_prim <- .Call("getColPrim", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(j)
                )

    return(col_prim)
}


#------------------------------------------------------------------------------#

getPrimStatGLPK <- function(lp) {

    prim_stat <- .Call("getPrimStat", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                  )

    return(prim_stat)
}


#------------------------------------------------------------------------------#

getDualStatGLPK <- function(lp) {

    dual_stat <- .Call("getDualStat", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                 )

    return(dual_stat)
}


#------------------------------------------------------------------------------#

getRowStatGLPK <- function(lp, i) {

    row_stat <- .Call("getRowStat", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(i)
                )

    return(row_stat)
}


#------------------------------------------------------------------------------#

getRowsStatGLPK <- function(lp) {

    rows_stat <- .Call("getRowsStat", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                 )

    return(rows_stat)
}


#------------------------------------------------------------------------------#

getRowPrimGLPK <- function(lp, i) {

    row_prim <- .Call("getRowPrim", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(i)
                )

    return(row_prim)
}


#------------------------------------------------------------------------------#

getRowsPrimGLPK <- function(lp) {

    rows_prim <- .Call("getRowsPrim", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                 )

    return(rows_prim)

}


#------------------------------------------------------------------------------#

getRowDualGLPK <- function(lp, i) {

    row_dual <- .Call("getRowDual", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(i)
                )

    return(row_dual)
}


#------------------------------------------------------------------------------#

getRowsDualGLPK <- function(lp) {

    rows_dual <- .Call("getRowsDual", PACKAGE = "glpkAPI",
              glpkPointer(lp)
        )

    return(rows_dual)
}

#------------------------------------------------------------------------------#

getColStatGLPK <- function(lp, j) {

    col_stat <- .Call("getColStat", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(j)
                )

    return(col_stat)
}


#------------------------------------------------------------------------------#

getColsStatGLPK <- function(lp) {

    cols_stat <- .Call("getColsStat", PACKAGE = "glpkAPI",
                      glpkPointer(lp) 
                 )

    return(cols_stat)
}


#------------------------------------------------------------------------------#

getColDualGLPK <- function(lp, j) {

    col_dual <- .Call("getColDual", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(j)
                )

    return(col_dual)

}


#------------------------------------------------------------------------------#

getColsDualGLPK <- function(lp) {

    cols_dual <- .Call("getColsDual", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                 )

    return(cols_dual)
}

#------------------------------------------------------------------------------#

getUnbndRayGLPK <- function(lp) {

    unbnd <- .Call("getUnbndRay", PACKAGE = "glpkAPI",
                   glpkPointer(lp)
             )

    return(unbnd)
}


#------------------------------------------------------------------------------#

solveInteriorGLPK <- function(lp) {

    ret <- .Call("solveInterior", PACKAGE = "glpkAPI",
              glpkPointer(lp)
           )

    return(ret)
}


#------------------------------------------------------------------------------#

getObjValIptGLPK <- function(lp) {

    obj <- .Call("getObjValIpt", PACKAGE = "glpkAPI",
                 glpkPointer(lp)
           )

    return(obj)
}


#------------------------------------------------------------------------------#

getSolStatIptGLPK <- function(lp) {

    stat <- .Call("getSolStatIpt", PACKAGE = "glpkAPI",
                  glpkPointer(lp)
            )

    return(stat)
}


#------------------------------------------------------------------------------#

getColsPrimIptGLPK <- function(lp) {

    cols_prim <- .Call("getColsPrimIpt", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                 )

    return(cols_prim)
}


#------------------------------------------------------------------------------#

getColPrimIptGLPK <- function(lp, j) {

    col_prim <- .Call("getColPrimIpt", PACKAGE = "glpkAPI",
                       glpkPointer(lp),
                       as.integer(j)
                 )

    return(col_prim)
}


#------------------------------------------------------------------------------#

getRowPrimIptGLPK <- function(lp, i) {

    row_prim <- .Call("getRowPrimIpt", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(i)
                )

    return(row_prim)

}


#------------------------------------------------------------------------------#

getRowsPrimIptGLPK <- function(lp) {

    rows_prim <- .Call("getRowsPrimIpt", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                 )

    return(rows_prim)

}

#------------------------------------------------------------------------------#

getRowDualIptGLPK <- function(lp, i) {

    row_dual <- .Call("getRowDualIpt", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(i)
                )

    return(row_dual)

}


#------------------------------------------------------------------------------#

getRowsDualIptGLPK <- function(lp) {

    rows_dual <- .Call("getRowsDualIpt", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                 )

    return(rows_dual)

}


#------------------------------------------------------------------------------#

getColDualIptGLPK <- function(lp, j) {

    col_dual <- .Call("getColDualIpt", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(j)
                )

    return(col_dual)
}


#------------------------------------------------------------------------------#

getColsDualIptGLPK <- function(lp) {

    cols_dual <- .Call("getColsDualIpt", PACKAGE = "glpkAPI",
                       glpkPointer(lp)
                 )

    return(cols_dual)
}


#------------------------------------------------------------------------------#

solveMIPGLPK <- function(lp) {

    ret <- .Call("solveMIP", PACKAGE = "glpkAPI",
                 glpkPointer(lp)
           )

    return(ret)
}


#------------------------------------------------------------------------------#

mipStatusGLPK <- function(lp) {

    stat <- .Call("mipStatus", PACKAGE = "glpkAPI",
                  glpkPointer(lp)
            )

    return(stat)
}


#------------------------------------------------------------------------------#

mipObjValGLPK <- function(lp) {

    obj_val <- .Call("mipObjVal", PACKAGE = "glpkAPI",
                     glpkPointer(lp)
               )

    return(obj_val)
}


#------------------------------------------------------------------------------#

mipRowValGLPK <- function(lp, i) {

    row_val <- .Call("mipRowVal", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(i)
                )

    return(row_val)
}


#------------------------------------------------------------------------------#

mipRowsValGLPK <- function(lp) {

    row_val <- .Call("mipRowsVal", PACKAGE = "glpkAPI",
                     glpkPointer(lp)
               )

    return(row_val)
}


#------------------------------------------------------------------------------#

mipColValGLPK <- function(lp, j) {

    col_val <- .Call("mipColVal", PACKAGE = "glpkAPI",
                      glpkPointer(lp),
                      as.integer(j)
                )

    return(col_val)
}


#------------------------------------------------------------------------------#

mipColsValGLPK <- function(lp) {

    col_val <- .Call("mipColsVal", PACKAGE = "glpkAPI",
                     glpkPointer(lp)
               )

    return(col_val)
}


#------------------------------------------------------------------------------#

getNumNnzGLPK <- function(lp) {

    nnz <- .Call("getNumNnz", PACKAGE = "glpkAPI",
                 glpkPointer(lp)
           )

    return(nnz)
}


#------------------------------------------------------------------------------#

getMatRowGLPK <- function(lp, i) {

    row_val <- .Call("getMatRow", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i)
        )

    return(row_val)
}


#------------------------------------------------------------------------------#

setMatRowGLPK <- function(lp, i, len, ind, val) {

    if ( (len == 0) && (is.null(ind)) ) {
        Cind <- as.null(ind)
    }
    else {
        if (ind[1] != 0) {
            Cind <- as.integer(append(ind, 0, 0))
        }
        else {
            Cind <- as.integer(ind)
        }
    }

    if ( (len == 0) && (is.null(val)) ) {
        Cval <- as.null(val)
    }
    else {
        if (val[1] != 0) {
            Cval <- as.numeric(append(val, 0, 0))
        }
        else {
            Cval <- as.numeric(val)
        }
    }

    check <- .Call("setMatRow", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(i),
              as.integer(len),
              Cind,
              Cval
    )

    return(check)
}


#------------------------------------------------------------------------------#

getMatColGLPK <- function(lp, j) {

    row_val <- .Call("getMatCol", PACKAGE = "glpkAPI",
                     glpkPointer(lp),
                     as.integer(j)
               )

    return(row_val)
}


#------------------------------------------------------------------------------#

setMatColGLPK <- function(lp, j, len, ind, val) {

    if ( (len == 0) && (is.null(ind)) ) {
        Cind <- as.null(ind)
    }
    else {
        if (ind[1] != 0) {
            Cind <- as.integer(append(ind, 0, 0))
        }
        else {
            Cind <- as.integer(ind)
        }
    }

    if ( (len == 0) && (is.null(val)) ) {
        Cval <- as.null(val)
    }
    else {
        if (val[1] != 0) {
            Cval <- as.numeric(append(val, 0, 0))
        }
        else {
            Cval <- as.numeric(val)
        }
    }

    check <- .Call("setMatCol", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(j),
              as.integer(len),
              Cind,
              Cval
    )

    return(check)

}


#------------------------------------------------------------------------------#

readMPSGLPK <- function(lp, fmt, fname) {

    check <- .Call("readMPS", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(fmt),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

readLPGLPK <- function(lp, fname) {

    check <- .Call("readLP", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

readProbGLPK <- function(lp, fname) {

    check <- .Call("readProb", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

writeMPSGLPK <- function(lp, fmt, fname) {

    check <- .Call("writeMPS", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(fmt),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

writeLPGLPK <- function(lp, fname) {

    check <- .Call("writeLP", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

writeProbGLPK <- function(lp, fname) {

    check <- .Call("writeProb", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

printSolGLPK <- function(lp, fname) {

    check <- .Call("printSol", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

readSolGLPK <- function(lp, fname) {

    check <- .Call("readSol", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

writeSolGLPK <- function(lp, fname) {

    check <- .Call("writeSol", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

printIptGLPK <- function(lp, fname) {

    check <- .Call("printIpt", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

readIptGLPK <- function(lp, fname) {

    check <- .Call("readIpt", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

writeIptGLPK <- function(lp, fname) {

    check <- .Call("writeIpt", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

printMIPGLPK <- function(lp, fname) {

    check <- .Call("printMIP", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

readMIPGLPK <- function(lp, fname) {

    check <- .Call("readMIP", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

writeMIPGLPK <- function(lp, fname) {

    check <- .Call("writeMIP", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.character(fname)
        )
    return(check)

}


#------------------------------------------------------------------------------#

versionGLPK <- function() {

    version <- .Call("version", PACKAGE = "glpkAPI")
    return(version)

}


#------------------------------------------------------------------------------#

bfExistsGLPK <- function(lp) {

    check <- .Call("bfExists", PACKAGE = "glpkAPI",
                   glpkPointer(lp)
             )

    return(check)
}


#------------------------------------------------------------------------------#

factorizeGLPK <- function(lp) {

    check <- .Call("factorize", PACKAGE = "glpkAPI",
                   glpkPointer(lp)    
             )

    return(check)
}


#------------------------------------------------------------------------------#

bfUpdatedGLPK <- function(lp) {

    check <- .Call("bfUpdated", PACKAGE = "glpkAPI",
                   glpkPointer(lp)    
             )

    return(check)
}


#------------------------------------------------------------------------------#

setBfcpGLPK <- function(lp, parm, val) {

    if (!identical(length(parm), length(val))) {
        stop("Arguments 'parm' and 'val' must have the same length!")
    }

    indi <- which(parm > 400 & parm < 500)
    indd <- which(parm > 500 & parm < 600)

    npari <- length(indi)
    npard <- length(indd)

    if (npari == 0) {
        parmi <- as.null(parm)
        vali  <- as.null(val)
    }
    else {
        parmi <- as.integer(parm[indi])
        vali  <- as.integer(val[indi])
    }

    if (npard == 0) {
        parmd <- as.null(parm)
        vald  <- as.null(val)
    }
    else {
        parmd <- as.integer(parm[indd])
        vald  <- as.numeric(val[indd])
    }

    invisible(
        .Call("setBfcp", PACKAGE = "glpkAPI",
              glpkPointer(lp),
              as.integer(npari),
              parmi,
              vali,
              as.integer(npard),
              parmd,
              vald
        )
    )
}


#------------------------------------------------------------------------------#

getBfcpGLPK <- function(lp) {

    parmB <- .Call("getBfcp", PACKAGE = "glpkAPI",
                   glpkPointer(lp)
             )

    return(parmB)
}


#------------------------------------------------------------------------------#

getBheadGLPK <- function(lp, k) {

    bh <- .Call("getBhead", PACKAGE = "glpkAPI",
                glpkPointer(lp),
                as.integer(k)
    )

    return(bh)
}


#------------------------------------------------------------------------------#

getRbindGLPK <- function(lp, i) {

    rh <- .Call("getRbind", PACKAGE = "glpkAPI",
                glpkPointer(lp),
                as.integer(i)
    )

    return(rh)
}


#------------------------------------------------------------------------------#

getCbindGLPK <- function(lp, j) {

    ch <- .Call("getCbind", PACKAGE = "glpkAPI",
                glpkPointer(lp),
                as.integer(j)
    )

    return(ch)
}


#------------------------------------------------------------------------------#

printRangesGLPK <- function(lp, numrc = 0, rowcol = NULL, fname = "sar.txt") {

    if ( (numrc == 0) || (is.null(rowcol)) ) {
        Crowcol <- as.null(rowcol)
    }
    else {
        Crowcol <- as.integer(c(0, rowcol))
    }

    sensit <- .Call("printRanges", PACKAGE = "glpkAPI",
                    glpkPointer(lp),
                    as.integer(numrc),
                    Crowcol,
                    as.character(fname)
              )

    return(sensit)
}


#------------------------------------------------------------------------------#

mplAllocWkspGLPK <- function(ptrtype = "tr_wksp") {

    wk <- .Call("mplAllocWksp", PACKAGE = "glpkAPI",
                as.character(ptrtype)
          )

    wkP <- trwks_Pointer(wk)

    return(wkP)
}


#------------------------------------------------------------------------------#

mplFreeWkspGLPK <- function(wk) {

    invisible(
        .Call("mplFreeWksp", PACKAGE = "glpkAPI",
              glpkPointer(wk)
        )
    )

}


#------------------------------------------------------------------------------#

mplReadModelGLPK <- function(wk, fname, skip) {

    check <- .Call("mplReadModel", PACKAGE = "glpkAPI",
              glpkPointer(wk),
              as.character(fname),
              as.integer(skip)
        )

    return(check)
}


#------------------------------------------------------------------------------#

mplReadDataGLPK <- function(wk, fname) {

    check <- .Call("mplReadData", PACKAGE = "glpkAPI",
              glpkPointer(wk),
              as.character(fname)
        )

    return(check)
}


#------------------------------------------------------------------------------#

mplGenerateGLPK <- function(wk, fname = NULL) {

    if (is.null(fname)) {
        Cfname <- as.null(fname)
    }
    else {
        Cfname <- as.character(fname)
    }

    check <- .Call("mplGenerate", PACKAGE = "glpkAPI",
              glpkPointer(wk),
              Cfname
        )

    return(check)
}


#------------------------------------------------------------------------------#

mplBuildProbGLPK <- function(wk, lp) {

    invisible(
        .Call("mplBuildProb", PACKAGE = "glpkAPI",
              glpkPointer(wk),
              glpkPointer(lp)
        )
    )
}


#------------------------------------------------------------------------------#

mplPostsolveGLPK <- function(wk, lp, sol) {

    check <- .Call("mplPostsolve", PACKAGE = "glpkAPI",
              glpkPointer(wk),
              glpkPointer(lp),
              as.integer(sol)
        )

    return(check)
}


