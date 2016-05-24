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

##
## Check write and read of vector in MAT5 format:
## 1) without compression
## 2) with compression
##

##
## vector: case-1
##
a1.exp <- 1:5
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a1.exp), filename=filename, compression=FALSE, version='MAT5')
a1.obs <- read.mat(filename)[['a']]
unlink(filename)
storage.mode(a1.obs) <- 'integer'
str(a1.obs)
stopifnot(identical(a1.obs, a1.exp))

# Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a1.exp), filename=filename, compression=TRUE, version='MAT5')
    a1.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    storage.mode(a1.zlib.obs) <- 'integer'
    str(a1.zlib.obs)
    stopifnot(identical(a1.zlib.obs, a1.exp))
}

##
## vector: case-2
##
a2.exp <- c(1,2,3,4,5)
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a2.exp), filename=filename, compression=FALSE, version='MAT5')
a2.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a2.obs)
stopifnot(identical(a2.obs, a2.exp))

# Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a2.exp), filename=filename, compression=TRUE, version='MAT5')
    a2.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a2.zlib.obs)
    stopifnot(identical(a2.zlib.obs, a2.exp))
}

##
## vector: case-3
##
a3.exp <- 1
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a3.exp), filename=filename, compression=FALSE, version='MAT5')
a3.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a3.obs)
stopifnot(identical(a3.obs, a3.exp))

# Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a3.exp), filename=filename, compression=TRUE, version='MAT5')
    a3.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a3.zlib.obs)
    stopifnot(identical(a3.zlib.obs, a3.exp))
}
