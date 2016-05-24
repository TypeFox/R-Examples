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

library(rmatio)

## If everything works as expected, no data are written to file...
filename <- tempfile(fileext = ".mat")

##
## Argument checking
##
## Check that the function write.mat stop if not the expected
## arguments are given
##

##
## 'filename' must be a character vector of length one
##
tools::assertError(write.mat(list(a=1:5), filename=NULL))
tools::assertError(write.mat(list(a=1:5), filename=5))
tools::assertError(write.mat(list(a=1:5), filename=c('a', 'b')))
tools::assertError(write.mat(list(a=1:5), filename=''))

##
## 'compression' must be a logical vector of length one
##
tools::assertError(write.mat(list(a=1:5), filename=filename, compression=NULL))
tools::assertError(write.mat(list(a=1:5), filename=filename, compression=5))
tools::assertError(write.mat(list(a=1:5), filename=filename, compression=c(TRUE, TRUE)))
tools::assertError(write.mat(list(a=1:5), filename=filename, compression=logical(0)))

##
## All values in the list must have a unique name
##
tools::assertError(write.mat(list(1:5), filename=filename, compression=FALSE))
tools::assertError(write.mat(list(a=1:5, 6:10), filename=filename, compression=FALSE))
tools::assertError(write.mat(list(a=1:5, a=6:10), filename=filename, compression=FALSE))

## Make sure the file is removed in case test failure and data are written...
unlink(filename)
