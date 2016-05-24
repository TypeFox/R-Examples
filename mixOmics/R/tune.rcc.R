# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


tune.rcc <-
function(...) UseMethod("tune.rcc")


tune.rcc <-
function(X, 
         Y, 
         grid1 = seq(0.001, 1, length = 5), 
         grid2 = seq(0.001, 1, length = 5), 
         validation = c("loo", "Mfold"), 
         folds, 
         M = 10, 
         plt = TRUE) 
{

    # validation des arguments #
	#--------------------------#
    if (length(dim(X)) != 2 || length(dim(Y)) != 2) 
        stop("'X' and/or 'Y' must be a numeric matrix.")
     
    X = as.matrix(X)
    Y = as.matrix(Y)
     
    if (!is.numeric(X) || !is.numeric(Y)) 
        stop("'X' and/or 'Y' must be a numeric matrix.")
		
    if (nrow(X) != nrow(Y)) 
        stop("unequal number of rows in 'X' and 'Y'.")
     
    validation = match.arg(validation)
    grid = expand.grid(grid1, grid2)
     
    if (validation == "loo") {
        cv.score = apply(grid, 1, function(lambda) 
                                    {loo(X, Y, lambda[1], lambda[2])})
    }
    else {
        nr = nrow(X)
        if (missing(folds)) folds = split(sample(1:nr), rep(1:M, length = nr))
        cv.score = apply(grid, 1, function(lambda) 
                                    {Mfold(X, Y, lambda[1], lambda[2], folds, M)})
    }
     
    cv.score.grid = cbind(grid, cv.score)
    mat = matrix(cv.score, nrow = length(grid1), ncol = length(grid2))
     
    if (isTRUE(plt)) {
        image.tune.rcc(list(grid1 = grid1, grid2 = grid2, mat = mat))
    }
     
    opt = cv.score.grid[cv.score.grid[, 3] == max(cv.score.grid[, 3]), ]
    cat("  lambda1 = ", opt[[1]], "\n", " lambda2 = ", opt[[2]], "\n",
	"CV-score = ", opt[[3]], "\n")
     
    out = list(opt.lambda1 = opt[[1]], opt.lambda2 = opt[[2]], 
	opt.score = opt[[3]], grid1 = grid1, grid2 = grid2, mat = mat)
     
    class(out) = "tune.rcc"	
    return(invisible(out))
}

