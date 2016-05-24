#  summaryOptsolClass.R
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


# summaryOptsolClass


#------------------------------------------------------------------------------#
#                            class definitions                                 #
#------------------------------------------------------------------------------#

setClass("summaryOptsol",
    representation(
        mod_id        = "character",     # model id
        mod_key       = "character",     # model key
        nzeros        = "integer",       # number of zeros in fld
        nnonzero      = "integer",       # number of non-zeros in fld
        mod_obj       = "numeric",       # objective coefficients of org. model %*% flux
        ex_met        = "character",     # metabolite id of exchange reaction
        ex_val        = "Matrix",        # flux value
        react_id      = "list",          # id's of limiting reactions
        chksol        = "checksol"       # solution status values
    )
)


#------------------------------------------------------------------------------#
#                              user constructors                               #
#------------------------------------------------------------------------------#

summaryOptsol <- function(opt,
                          mod,
                          perc = 1,
                          tol = SYBIL_SETTINGS("TOLERANCE")) {

    stopifnot(is(opt, "optsol"),
              is(mod, "modelorg"),
              identical(mod_id(opt), mod_id(mod)),
              length(perc) == 1,
              (perc <= 1 || perc >= 0))

    if (any(is.na(fluxes(opt)))) {
        #cat("\noptimal values of model objective function:\n")
        os <- summary(mod_obj(opt))
    }
    else {
        # zero and non-zero elements
        zeros <- nnzero(fluxdist(opt))
        nzero <- num_of_fluxes(fluxdist(opt)) - nnzero(fluxdist(opt))
        
        # exchange fluxes
        exch <- findExchReact(mod)

        # limiting reactions
        fl <- fluxes(opt)[fldind(opt), , drop = FALSE]
        
        lt <- apply(fl, 2, function(x) which( ( (x <= (lowbnd(mod)*perc)) |
                                                (x >= (uppbnd(mod)*perc)) )
                                               & abs(x) > tol ) )

        ltstr <- lapply(lt, function(x) react_id(mod)[x])

        # summary object
        os <- new("summaryOptsol",
                  mod_id      = mod_id(opt),
                  mod_key     = mod_key(opt),
                  mod_obj     = mod_obj(opt),
                  nnonzero    = nzero,
                  nzeros      = zeros,
                  ex_met      = met_id(exch),
                  ex_val      = fluxes(opt)[fldind(opt)[react_pos(exch)], , drop = FALSE],
                  react_id    = ltstr,
                  chksol      = checkOptSol(opt))
    }        

    return(os)
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# mod_id
setMethod("mod_id", signature(object = "summaryOptsol"),
          function(object) {
              return(object@mod_id)
          }
)

setReplaceMethod("mod_id", signature = (object = "summaryOptsol"),
                 function(object, value) {
                     object@mod_id <- value
                     return(object)
                 }
)


# mod_key
setMethod("mod_key", signature(object = "summaryOptsol"),
          function(object) {
              return(object@mod_key)
          }
)

setReplaceMethod("mod_key", signature = (object = "summaryOptsol"),
                 function(object, value) {
                     object@mod_key <- value
                     return(object)
                 }
)


# mod_obj
setMethod("mod_obj", signature(object = "summaryOptsol"),
          function(object) {
              return(object@mod_obj)
          }
)

setReplaceMethod("mod_obj", signature = (object = "summaryOptsol"),
                 function(object, value) {
                     object@mod_obj <- value
                     return(object)
                 }
)


# nzeros
setMethod("nzeros", signature(object = "summaryOptsol"),
          function(object) {
              return(object@nzeros)
          }
)

# nnzero
setMethod("nnzero", signature(x = "summaryOptsol"),
          function(x) {
              return(x@nnonzero)
          }
)

# ex_met
setMethod("ex_met", signature(object = "summaryOptsol"),
          function(object) {
              return(object@ex_met)
          }
)

# ex_val
setMethod("ex_val", signature(object = "summaryOptsol"),
          function(object) {
              return(object@ex_val)
          }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

# show
setMethod("show", signature(object = "summaryOptsol"),
    function(object) {
        
        nf <- ncol(ex_val(object))
        
        cat("flux distribution:\n")
        cat(" number of elements:         ", nnzero(object) + nzeros(object), "\n")
        cat(" number of zeros:            ", nnzero(object), "\n")
        cat(" number of non zero elements:", nzeros(object),   "\n\n")
        cat("exchange metabolites:\n")
        print(ex_met(object), quote = FALSE, justify = "right")
        cat("\nsubstrates (-) and products (+):\n")
        if (nf < 21) {
            printExchange(object, dense = FALSE)
        }
        else if (nf < 66) {
            printExchange(object, dense = TRUE)
        }
        else {
            cat(" use method ", sQuote("printExchange()"), "\n", sep = "")
        }
        cat("\nlimiting reactions:\n")
        avglt <- floor(sum(sapply(object@react_id, length))/length(object@react_id))
        msg <- sprintf(ngettext(avglt,
           " there is about one limiting reaction per optimization;\n",
           " there are about %s limiting reactions per optimization;\n"), avglt)
        cat(msg)
        cat(" use method", sQuote("printReaction()"), "to see more details\n")
        cat("\noptimal values of model objective function:\n")
        show(summary(mod_obj(object)))
        cat("\nsummary of optimization process:\n")
        show(object@chksol)
    }
)


# draw a histogramm (package lattice)
setMethod("plot", signature(x = "summaryOptsol", y = "missing"),
          function(x, y,
                   col = "grey",
                   main = "",
                   xlab = "optimal values of model objective function",
                   ...) {

              histogram(mod_obj(x), col = col, main = main, xlab = xlab, ...)
              
          }
)


# printMetabolite
setMethod("printReaction", signature(object = "summaryOptsol",
                                     mod = "modelorg"),
          function(object, mod, j, ...) {

              stopifnot(identical(mod_id(object), mod_id(mod)))

              if (missing(j)) {
                  j <- 1:length(object@react_id)
              }
              
              ltstr <- lapply(object@react_id[j],
                              function(x) {
                                  msg <- sprintf(ngettext(length(x),
                                     "\none reaction:\n",
                                     "\n%s reactions:\n"), length(x))
                                  cat(msg)
                                  printReaction(mod, react = x, ...)
                              })
              names(ltstr) <- j

              return(invisible(ltstr))
          }
)


# printExchange
setMethod("printExchange", signature(object = "summaryOptsol"),
          function(object, i, j, dense = FALSE) {

              if(missing(i)) {
                  i <- 1:nrow(object@ex_val)
              }
              if(missing(j)) {
                  j <- 1:ncol(object@ex_val)
              }

              mat <- .recodeMatrix(object@ex_val[i, j, drop = FALSE])
              rownames(mat) <- sprintf("%10s", object@ex_met[i])
              cw <- paste("%-", floor(log10(ncol(mat)))+1, "i", sep = "")
              colnames(mat) <- sprintf(cw, j)

              if (isTRUE(dense)) {
                  vec <- apply(mat, 1, paste, collapse = "", sep = "")
                  dm  <- sprintf("%10s %s", rownames(mat), vec)
                  cat(dm, sep = "\n")
              }
              else {
                  print(mat, quote = FALSE)
              }

              return(invisible(mat))
          }
)


# image
setMethod("image", signature(x = "summaryOptsol"),
          function(x,
                   cuts = 2,
                   xlab = "optimization no.",
                   ylab = "exchange reaction no.",
                   useAbs = TRUE,
                   sub = NULL,
                   printOut = TRUE, ...) {

              mat <- .recodeMatrix(x@ex_val, signs = c(-1, 0, 1))

              imgmat <- image(Matrix(mat),
                              cuts = 2,
                              xlab = xlab,
                              ylab = ylab,
                              useAbs = useAbs,
                              sub = sub, ...)

              if (isTRUE(printOut)) {
                  print(imgmat)
              }

              return(invisible(imgmat))
          }
)
