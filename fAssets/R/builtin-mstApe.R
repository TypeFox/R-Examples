
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA
    
    
################################################################################
# FUNCTION:                 DESCRIPTION:
#   .mst                     Minimum spanning tree
#   .sortIndexMST
#   .mstPlot
#   .nsca
################################################################################


# Rmetrics:
#   Note that covRobust is not available on Debian as of 2009-04-28. 
#   To run these functions under Debian/Rmetrics we have them    
#   implemented here as a builtin.
#   We also made modifications for tailored usage with Rmetrics. 


# Package: ape
# Version: 2.3
# Date: 2009-03-30
# Title: Analyses of Phylogenetics and Evolution
# Author: Emmanuel Paradis, Ben Bolker, Julien Claude, Hoa Sien Cuong,
#   Richard Desper, Benoit Durand, Julien Dutheil, Olivier Gascuel,
#   Gangolf Jobb, Christoph Heibl, Daniel Lawson, Vincent Lefort,
#   Pierre Legendre, Jim Lemon, Yvonnick Noel, Johan Nylander,
#   Rainer Opgen-Rhein, Korbinian Strimmer, Damien de Vienne
# Maintainer: Emmanuel Paradis <Emmanuel.Paradis@ird.fr>
# Depends: R (>= 2.6.0)
# Suggests: gee
# Imports: gee, nlme, lattice
# ZipData: no
# Description: ape provides functions for reading, writing, plotting, and
#   manipulating phylogenetic trees, analyses of comparative data
#   in a phylogenetic framework, analyses of diversification and
#   macroevolution, computing distances from allelic and nucleotide
#   data, reading nucleotide sequences, and several tools such as
#   Mantel's test, computation of minimum spanning tree, the
#   population parameter theta based on various approaches,
#   nucleotide diversity, generalized skyline plots, estimation of
#   absolute evolutionary rates and clock-like trees using mean
#   path lengths, non-parametric rate smoothing and penalized
#   likelihood, classifying genes in trees using the
#   Klastorin-Misawa-Tajima approach. Phylogeny estimation can be
#   done with the NJ, BIONJ, and ME methods.
# License: GPL (>= 2)
# URL: http://ape.mpl.ird.fr/
# Packaged: Mon Mar 30 08:46:28 2009; paradis
# Repository: CRAN
# Date/Publication: 2009-03-30 06:56:17
 

# ------------------------------------------------------------------------------


.mst <-  
function(X)
{
    # Description:
    #   The function mst finds the minimum spanning tree between  
    #   a set of observations using a matrix of pairwise distances.
    
    # Authors:
    #   Original Code: Yvonnick Noel, Julien Claude, and Emmanuel Paradis
    
    # Source:
    #   Contributed R-packe "ape".
    
    # FUNCTION:
    
    # Minimum Spanning Tree:
    if (class(X) == "dist") X = as.matrix(X)
    n = dim(X)[1]
    N = matrix(0, n, n)
    tree = NULL
    large.value = max(X) + 1
    diag(X) = large.value
    index.i = 1

    for (i in 1:(n - 1)) {
        tree = c(tree, index.i)
        # calcul les minimum par colonne
        m = apply(as.matrix(X[, tree]), 2, min)  
        a = .sortIndexMST(X[, tree])[1, ]
        b = .sortIndexMST(m)[1]
        index.j = tree[b]
        index.i = a[b]

        N[index.i, index.j] = 1
        N[index.j, index.i] = 1

        for (j in tree) {
            X[index.i, j] = large.value
            X[j, index.i] = large.value
        }
    }
    dimnames(N) = dimnames(X)
    class(N) = "mst"
    
    # Return Value:
    return(N)
}


# ------------------------------------------------------------------------------


.sortIndexMST <-  
function(X)
{
    # Function returning an index matrix for an increasing sort
    if(length(X) == 1) return(1)                  # sorting a scalar?
    if(!is.matrix(X)) X = as.matrix(X)            # force vector into matrix
    
    # n = nrow(X)
    apply(X, 2, function(v) order(rank(v)))       # find the permutation
}


# ------------------------------------------------------------------------------


.mstPlot <-  
function (x, graph = "circle", x1 = NULL, x2 = NULL, ...) 
{
    # Description:
    #   Plots the minimum spanning tree showing the links 
    #   where the observations are identified by their numbers.
    
    # FUNCTION:
    
    # Plot:
    n = nrow(x)
    if (is.null(x1) || is.null(x2)) {
        if (graph == "circle") {
            ang = seq(0, 2 * pi, length = n + 1)
            x1 = cos(ang)
            x2 = sin(ang)
            plot(x1, x2, 
                type = "n", 
                xlab = "", ylab = "", xaxt = "n", 
                yaxt = "n", bty = "n", ...)
        }
        if (graph == ".nsca") {
            XY = .nsca(x)
            x1 = XY[, 1]
            x2 = XY[, 2]
            xLim = c(min(x1) - 0.25 * diff(range(x1)), max(x1))
            plot(XY, 
                type = "n", 
                xlim = xLim,
                xlab = "", # "\".nsca\" -- axis 1", 
                ylab = "", # "\".nsca\" -- axis 2", 
                xaxt = "n", yaxt = "n", col = "red", 
                ...)
            # Legend:
            Names = colnames(x)
            legendtext = paste(1:length(Names), Names, sep = "-")
            legendtext = substr(legendtext, 1, 8)
            legend("topleft", legend = legendtext, bty = "n", cex = 0.8)
        }
    } else {
        plot(x1, x2, type = "n",  
            xlab = deparse(substitute(x1)), 
            ylab = deparse(substitute(x2)), ...)
    }
    
    for (i in 1:n) {
        w1 = which(x[i, ] == 1)
        segments(x1[i], x2[i], x1[w1], x2[w1], lwd = 2)
    }
    
    points(x1, x2, pch = 21, col = "red", bg = "black", cex = 4)
    text(x1, x2, 1:n, col = "white", cex = 0.7)
}


# ------------------------------------------------------------------------------


.nsca <-  
function(A)
{
    # FUNCTION:
    
    Dr = apply(A, 1, sum)
    Dc = apply(A, 2, sum)

    eig.res = eigen(diag(1 / sqrt(Dr)) %*% A %*% diag(1 / sqrt(Dc)))
    r = diag(1 / Dr) %*% (eig.res$vectors)[, 2:4]
    
    # The next line has been changed by EP (20-02-2003), since
    # it does not work if 'r' has no dimnames already defined
    # dimnames(r)[[1]] = dimnames(A)[[1]]
    rownames(r) = rownames(A)
    
    # Return Value:
    r
}


################################################################################

