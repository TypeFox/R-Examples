#  fluxDistributionClass.R
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


# fluxDistributionClass

#------------------------------------------------------------------------------#
#                  definition of the class fluxDistribution                    #
#------------------------------------------------------------------------------#

setClass("fluxDistribution",
         representation(
              fluxes = "Matrix",
              num_of_fluxes = "integer"
         )
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

fluxDistribution <- function(fluxes = NA, nrow = 1, ncol = 1) {

    if (is(fluxes, "matrix")) {
        Mfluxes <- as(fluxes, "Matrix")
    }
    else if (is(fluxes, "Matrix")) {
        Mfluxes <- fluxes
    }
    else {
        Mfluxes <- Matrix::Matrix(fluxes, nrow = nrow, ncol = ncol)
    }

    num_of_fluxes <- length(Mfluxes)

    obj <- new("fluxDistribution",
               fluxes = Mfluxes,
               num_of_fluxes = as.integer(num_of_fluxes))
    return(obj)
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# fluxes
setMethod("fluxes", signature(object = "fluxDistribution"),
          function(object) {
              return(object@fluxes)
          }
)

setReplaceMethod("fluxes", signature = (object = "fluxDistribution"),
                 function(object, value) {
                     object@fluxes <- value
                     return(object)
                 }
)


# number of fluxes
setMethod("num_of_fluxes", signature(object = "fluxDistribution"),
          function(object) {
              return(object@num_of_fluxes)
          }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

# nnzero
setMethod("nnzero", signature(x = "fluxDistribution"),
    function(x) {

        nnz <- Matrix::nnzero(fluxes(x))
        return(nnz)

    }
)


# nvar
setMethod("nvar", signature(object = "fluxDistribution"),
    function(object) {

        nr <- nrow(fluxes(object))
        return(nr)

    }
)


# [
setMethod("[", signature(x = "fluxDistribution"),
    function(x, i, j, ..., drop = FALSE) {

        if (missing(i)) {
            i <- c(1:nrow(x@fluxes))
        }
        
        if (missing(j)) {
            j <- c(1:ncol(x@fluxes))
        }


        newfld <- fluxDistribution(x@fluxes[i,j],
                                   nrow = length(i),
                                   ncol = length(j))

        return(newfld)

    }

)


setMethod("plot", signature(x = "fluxDistribution", y = "missing"),
          function(x, y, ordReact, ordMut, todo = "absdev", ...) {
              
              if (any(is.na(fluxes(x)))) {
                  stop("fluxdistribution is required.")
              }
              
              checkPackage <- require("lattice")
              if(!isTRUE(checkPackage)) {
                  stop("Package lattice not found.")
              }
              checkPackage <- require("grid")
              if(!isTRUE(checkPackage)) {
                  stop("Package grid not found.")
              }

#g <- apply(subSys(bla), 1, function(x) colnames(subSys(bla))[x])
#order(g)

              # Arguments ordReact and ordMut can be numeric vectors or lists:
              #
              # numeric vectors:
              # a vector of length n (number of reactions) if ordReact, or
              # a vector of length number of mutations (number of problems to
              # solve, minus wild type).
              # The vector contains a desired permutation of indices.
              #
              # list containing two elements 'cat' and 'bin':
              # Element 'cat' is a character vector containing category names.
              # Element 'bin' is a numeric vector connecting the elements in
              # 'cat' with the reactions or mutations.
              
              
              if (missing(ordReact)) {
                  cat <- NA
                  fld <- fluxes(x)
              }
              else {
                  if (is.list(ordReact)) {
                      cat <- ordReact$cat[ordReact$bin]
                      fld <- fluxes(x)[ordReact$bin, ]
                  }
                  else if (is.numeric(ordReact)) {
                      cat <- NA
                      fld <- fluxes(x)[ordReact, ]
                  }
                  else {
                      stop("argument 'ordReact' must be numeric or list.")
                  }
              }


              # first column is always wild type
              ref <- fld[, 1]
              bin <- fld[ ,(2 : ncol(fld)-1)]
              #print(is(mat))
              
              mat <- switch(todo,
                  "reldev" = {
                      bin <- abs(ref) < SYBIL_SETTINGS("TOLERANCE")
                  },
                  "absdev" = {
                      fld - ref
                  },
                  {
                      fld
                  }
              )

              
              #nz <- abs(mat) < SYBIL_SETTINGS("TOLERANCE")
              #mat[nz] <- 0
              
              # 20 breaks
              br <- c(-1e+03, -1e+02, -1e+01, -1e+00,
                      -1e-01, -1e-02, -1e-03, -1e-04, -1e-05, -1e-06,
                       #0e+00,
                       1e-06,  1e-05,  1e-04,  1e-03,  1e-02,  1e-01,
                       1e+00,  1e+01,  1e+02,  1e+03)
              
              
              
              ylabels = c(expression(-10^3),
                          expression(-10^2),
                          expression(-10^1),
                          expression(-10^0),
                          expression(-10^-1),
                          expression(-10^-2),
                          expression(-10^-3),
                          expression(-10^-4),
                          expression(-10^-5),
                          expression(-10^-6),
                          expression(~~10^-6),
                          expression(~~10^-5),
                          expression(~~10^-4),
                          expression(~~10^-3),
                          expression(~~10^-2),
                          expression(~~10^-1),
                          expression(~~10^0),
                          expression(~~10^1),
                          expression(~~10^2),
                          expression(~~10^3)
                          )


              #
              ## colors are from RColorBrewer:
              ## rev(brewer.pal(9, "Blues")) for negative values
              ## brewer.pal(9, "Reds")  for positive values
              
              # colors: blue negative vale, red positive value
              neg <- c("#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6",
                       "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")

              # green
              #neg <- c("#00441B", "#006D2C", "#238B45", "#41AB5D", "#74C476",
              #         "#A1D99B", "#C7E9C0", "#E5F5E0", "#F7FCF5")

              pos <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
                       "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
              
              colors <- c(neg, "#FFFFFF", pos)


              #layout(matrix(c(1, 2), ncol = 2), widths = c(1,4))
              #image(z = t(ref), axes = FALSE, breaks = br, col = colors, ...)
              #image(z = t(mat), breaks = br, col = colors, ...)
             
              levelplot(x = mat,
                        xlab = "Reaction No.",
                        ylab = "Gene/Flux No.",
                        at = br,
                        col.regions = colors,
                        colorkey = list(at = seq(-1e+03, 1e+03, length = 20),
                                        labels = ylabels
                        ),
                        panel = function(...) {
                        
#                            grid.segments(x0 = c(2, 2, 10), y0 = c(1, 5, 5),
#                                          x1 = c(2, 10, 10), y1 = c(5, 5, 1),
#                                          default.units = "native")
                            
                            grid::grid.points(5, 5, pch = 16, size=grid::unit(5, "mm"))
#                            default.units = "native")
                            
                            panel.levelplot(...)

                        },
                        
                        
                        ...)
                        
            
          }
)
