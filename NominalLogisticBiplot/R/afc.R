# file NominalLogisticBiplot/R/afc.R
# copyright (C) 2012-2013 J.L. Vicente-Villardon and J.C. Hernandez
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#Function that calculates simple correspondence analysis for a data set.
#----------------------Parameters--------------
  #x: data set with nominal variables
  #dim: number of dimensions for the solution
  #alpha: this parameter is used to calculate row and column markers
afc <- function(x, dim = 2, alpha = 1) {

    if (is.data.frame(x))
      x = as.matrix(x)
    if (is.null(rownames(x)))
      rownames(x) <- rownames(x, do.NULL = FALSE, prefix = "I")
    RowNames = rownames(x)
    if (is.null(colnames(x)))
      colnames(x) <- colnames(x, do.NULL = FALSE, prefix = "V")
    VarNames = colnames(x)
    DimNames = "Dim 1"
    if(dim > 1){
      for (i in 2:dim) DimNames = c(DimNames, paste("Dim", i))
    }
    afc=list()
    afc$Title = "Correspondence Analysis"
    afc$Non_Scaled_Data = x
    afc$Minima = apply(x, 2, min)
    afc$Maxima = apply(x, 2, max)
    afc$Initial_Transformation = "Correspondence Analysis"
    n = dim(x)[1]
    p = dim(x)[2]
    nt = sum(x)
    x = x/nt
    dr = matrix(rowSums(x), n, 1)
    dc = matrix(colSums(x), p, 1)
    x = x - dr %*% t(dc)
    afc$Scaled_Data = x
    Dr = diag(1, n, n)
    diag(Dr) <- 1/sqrt(dr)
    Dc = diag(1, p, p)
    diag(Dc) <- 1/sqrt(dc)
    x = Dr %*% x %*% Dc
    UDV = svd(x)
    r = min(c(n, p))
    d = UDV$d[1:r]
    iner = ((d^2)/sum((d^2))) * 100
    U = diagonal(sqrt(1/dr)) %*% UDV$u[, 1:r]
    V = diagonal(sqrt(1/dc)) %*% UDV$v[, 1:r]
    D = diagonal(d)
    A = U %*% D
    B = V %*% D
    sf = rowSums(A^2)
    cf = solve(diagonal(sf)) %*% (A^2)
    sc = rowSums(B^2)
    cc = solve(diagonal(sc)) %*% (B^2)
    if (alpha <= 1) {
      A = U %*% diagonal(d ^ alpha)
      B = V %*% diagonal(d^(1 - alpha))
      afc$Expected=A[, 1:dim] %*% t(B[, 1:dim])
    }
    afc$nrows = n
    afc$ncols = p
    afc$dim = dim
    afc$EigenValues = d^2
    afc$Inertia = iner
    afc$CumInertia = cumsum(iner)
    afc$RowCoordinates = as.matrix(A[, 1:dim])
    rownames(afc$RowCoordinates)=RowNames
    colnames(afc$RowCoordinates)=DimNames
    afc$ColCoordinates = as.matrix(B[, 1:dim])
    rownames(afc$ColCoordinates)=VarNames
    colnames(afc$ColCoordinates)=DimNames
    afc$RowContributions = round(as.matrix(cf[, 1:dim]) * 100, digits = 2)
    rownames(afc$RowContributions)=RowNames
    colnames(afc$RowContributions)=DimNames
    afc$ColContributions = round(as.matrix(cc[, 1:dim]) * 100, digits = 2)
    rownames(afc$ColContributions)=VarNames
    colnames(afc$ColContributions)=DimNames
    afc$Scale_Factor = 1
    class(afc) = "afc.sol"
    return(afc)
}