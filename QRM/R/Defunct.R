## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


mk.returns <- function(tsdata, type = "log"){
  .Defunct(new = "returns()", package = "timeSeries",
           msg = "The function mk.returns() previously contained in the package QRMlib has been removed. Use returns() contained in the package timeSeries instead.")
}
plotMultiTS <- function(tS, colvec = 1:ncol(tS), type = "l", ltypvec = 1, lwdvec = 1, yrange, format, at, reference.grid, ...){
  .Defunct(new = "plot()", package = "timeSeries",
           msg = "The function plotMultiTS() previously contained in the package QRMlib has been removed. Use plot-method contained in the package timeSeries instead.")
}
signalSeries <- function(data, positions., units, units.position, from = 1, by = 1){
  .Defunct(new = "series()", package = "timeSeries",
           msg = "The function signalSeries() previously contained in the package QRMlib has been removed. Use series() contained in the package timeSeries instead.")
}
aggregateQuarterlySeries <- function(timeseries, FUNC = colSums){
  .Defunct(new = "aggregate()", package = "timeSeries",
           msg = "The function aggregateQuarterlySeries() previously contained in the package QRMlib has been removed. Use aggregate() contained in the package timeSeries instead.")
}
aggregateMonthlySeries <- function(timeseries, FUNC = colSums){
  .Defunct(new = "aggregate()", package = "timeSeries",
           msg = "The function aggregateMonthlySeries() previously contained in the package QRMlib has been removed. Use aggregate() contained in the package timeSeries instead.")
}
aggregateWeeklySeries <- function(timeseries, FUNC = colSums){
  .Defunct(new = "aggregate()", package = "timeSeries",
           msg = "The function aggregateWeeklySeries() previously contained in the package QRMlib has been removed. Use aggregate() contained in the package timeSeries instead.")
}
aggregateSignalSeries <- function(x, pos, AGGFUNC, together = FALSE, drop.empty = TRUE, include.ends = FALSE, adj,offset, colnames, by){
  .Defunct(new = "aggregate()", package = "timeSeries",
           msg = "The function aggregateSignalSeries() previously contained in the package QRMlib has been removed. Use aggregate() contained in the package timeSeries instead.")
}
ConvertDFToTimeSeries <- function(dataframe){
  .Defunct(new = "timeSeries()", package = "timeSeries",
           msg = "The function ConvertDFToTimeSeries() previously contained in the package QRMlib has been removed. Use aggregate() contained in the package timeSeries instead.")
}
CovToCor <- function(mat){
  .Defunct(new = "cov2cor()", package = "stats",
           msg = "The function CovToCor() previously contained in the package QRMlib has been removed. Use cov2cor() contained in the package stats instead.")
}
symmetrize <- function(matrix){
  .Defunct(new = "forceSymmetric()", package = "Matrix",
           msg = "The function symmetrize() previously contained in the package QRMlib has been removed. Use forecSymmetric() contained in the package Matrix instead.")
}
hessb <- function(f, x, ep = 0.0001, ...){
  .Defunct(new = "hessb()", package = "numDeriv",
           msg = "The function hessb() previously contained in the package QRMlib has been removed. Use hessian() contained in the package numDeriv instead.")
}
fit.Archcopula2d <- function(Udata, name){
  .Defunct(new = "fit.AC()", package = "QRM",
           msg = "The function fit.Archcopula2d() previously contained in the package QRMlib has been removed. Use fit.AC() instead.")
}
besselM3 <- function(lambda = 9/2, x = 2, log = FALSE){
  .Defunct(new = "besselK()", package = "base",
           msg = "The function besselM3() previously contained in the package QRMlib has been removed. Use besselK() instead.")
}
psifunc <- function(x = 2, logvalue = FALSE){
  .Defunct(new = "psi()", package = "gsl",
           msg = "The function psifunc() previously contained in the package QRMlib has been removed. Use psi() instead.")
}
kurtosisSPlus <- function(x, na.rm, method = "fisher"){
  .Defunct(new = "kurtosis()", package = "timeDate",
           msg = "The function kurtosisSPlus() previously contained in the package QRMlib has been removed. Use kurtosis() instead.")
}
fit.tcopula.rank <- function(Udata, method = "Kendall"){
  .Defunct(new = "fit.tcopula()", package = "QRM",
           msg = "The function fit.tcopula.rank() previously contained in the package QRMlib has been removed. Use fit.tcopula() with appropriate method selection instead.")
}
fit.GPDb <- function(data, threshold = NA, nextremes = NA, method = "ml", information = "observed"){
  .Defunct(new = "fit.GPD()", package = "QRM",
           msg = "The function fit.GPDb() previously contained in the package QRMlib has been removed. Use fit.GPD() with appropriate method selection instead.")
}
lbeta <- function(a, b){
  .Defunct(new = "lbeta()", package = "base",
           msg = "The function lbeta() previously contained in the package QRMlib has been removed. Use lbeta() of the base package instead.")
}
