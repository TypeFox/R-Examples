### ----------------------------------------------------------------------------
### This file is part of boolean3
###
### Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
###
### boolean3 represents a substantial re-write of the original boolean package
### developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
### was developed under the direction of Bear Braumoeller and with support from
### The Ohio State University's College of Social and Behavioral Sciences.
###
### boolean3 and is free software: you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the Free
### Software Foundation, either version 3 of the License, or (at your option)
### any later version.
###
### This program is distributed in the hope that it will be useful, but WITHOUT
### ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
### FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
### more details.
###
### You should have received a copy of the GNU General Public License along with
### this program.  If not, see <http://www.gnu.org/licenses/>.
###
### ----------------------------------------------------------------------------


##' Boolean Model Class
##'
##' An object of class \code{boolean} as returned from \code{\link{boolprep}}.
##' @name boolean-class
##' @docType class
##' @return A boolean class contains the following items:
##' \describe{
##' \item{call}{call to estimate the model.}
##' \item{links}{link functions used to estimate the model.}
##' \item{model}{list, specifying the structure of the boolean model (used
##'              internally).}
##' \item{N}{number of observations.}
##' \item{k}{number of covariates.}
##' \item{coef.labels}{preformatted labels for the estimated coefficients.}
##' \item{coef.idx}{indicies for the coefficients in the model matrix (used
##'                 internally).}
##' \item{response}{response vector.}
##' \item{frame}{model frame used to estimate the model.}
##' }
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
NULL

##' A boolfit class
##'
##' An object of class boolfit.
##'
##' @name boolfit-class
##' @docType class
##' @return A boolfit object is a component added to a
##' \code{\link{boolean-class}} object (as \code{model.fit}) once the model is
##' fit with \code{\link{boolean}}. It consists of a list of results returned by
##' the particular optimization routines used to fit the model.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
NULL

##' A boolsum class
##'
##' An object of class boolsum.
##'
##' @name boolsum-class
##' @docType class
##' @return The object returned by \code{\link{summary.boolean}}.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
NULL

##' A boolprof class
##'
##' An object of class boolprof.
##'
##' @name boolprof-class
##' @docType class
##' @return The object returned by \code{\link{boolprof}}. It's comprised of the
##' following items:
##' \describe{
##' \item{est}{a data.frame of three columns: \code{llik}, the calculated
##'            log-likelihood; \code{x}, the coefficient value; \code{var}, the
##'            preformatted variable name.}
##' \item{coef.labels}{preformatted coefficient labels.}
##' \item{default.plot}{a lattice plot produced by \code{xyplot}.}
##' }
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
NULL

##' A boolprob class
##'
##' An object of class boolprob.
##'
##' @name boolprob-class
##' @docType class
##' @return The object returned by \code{\link{boolprob}}. It's comprised of the
##' following items:
##' \describe{
##' \item{est}{a data.frame of five columns: \code{llik}, the calculated range
##'            for the predicted probability
##'            \code{x}, the covariate value; \code{coef}, the
##'            preformatted variable name.}
##' \item{coef.labels}{preformatted coefficient labels.}
##' \item{default.plot}{a lattice plot produced by \code{xyplot}.}
##' } 
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
NULL
