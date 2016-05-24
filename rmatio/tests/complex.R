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
## Check write and read of complex in MAT5 format:
## 1) without compression
## 2) with compression
##

##
## complex: case-1
##
a1.exp <- array(complex(real=1:20,imaginary=21:40), c(4,5))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a1.exp), filename=filename, compression=FALSE, version='MAT5')
a1.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a1.obs)
stopifnot(identical(a1.obs, a1.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a1.exp), filename=filename, compression=TRUE, version='MAT5')
    a1.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a1.zlib.obs)
    stopifnot(identical(a1.zlib.obs, a1.exp))
}
