## rmatio, a R interface to the C library matio, MAT File I/O Library.
## Copyright (C) 2013-2014  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## rmatio is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' \pkg{rmatio}: reading and writing Matlab MAT files from R
##'
##' Reading and writing Matlab MAT files from R
##'
##' \code{rmatio} supports reading MAT version 4, MAT version 5 and
##' MAT compressed version 5.
##'
##' \code{rmatio} can write version 5 MAT files and version 5 files
##' with variable compression.
##' @import Matrix
##' @import methods
##' @import lattice
##' @name rmatio
##' @references \itemize{
##'   \item Christopher C. Hulbert, MATIO User Manual for version 1.5.2.\cr
##'   \url{http://sourceforge.net/projects/matio/files/matio/1.5.2/matio_user_guide.pdf}
##'
##'   \item The MathWorks Inc., MATLAB - MAT-File Format, version R2013b, September 2013.\cr
##'   \url{http://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf}
##' }
##' @docType package
NULL
