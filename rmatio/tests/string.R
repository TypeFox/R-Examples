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
## Check write and read of character strings in MAT5 format:
## 1) without compression
## 2) with compression
##

##
## Note:
## If the character vector contains elements of differents lengths i.e.
## a2.in <- list(c("a", "bb"), c("c", "dd"))
## then the expected result is not identical, since each element
## is saved in a cell and the expected result of a2.in is therefore
## a2.exp <- list("a", "bb")
##

##
## string: case-1
##
a1.exp <- c("abcdefghijklmnopqrstuvwxyz",
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
          "1234567890!@#$%^&*()-_=+`~", #
          "[{]}\\|;:'\",<.>/?          ")
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

##
## string: case-2
##
a2.in <- c("a", "bb")
a2.exp <- list("a", "bb")
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a2.in), filename=filename, compression=FALSE, version='MAT5')
a2.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a2.obs)
stopifnot(identical(a2.obs, a2.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a2.in), filename=filename, compression=TRUE, version='MAT5')
    a2.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a2.zlib.obs)
    stopifnot(identical(a2.zlib.obs, a2.exp))
}

##
## string: case-3
##
a3.in <- list(y=c("a", "bb"))
a3.exp <- list(y=list("a", "bb"))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a3.in), filename=filename, compression=FALSE, version='MAT5')
a3.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a3.obs)
stopifnot(identical(a3.obs, a3.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a3.in), filename=filename, compression=TRUE, version='MAT5')
    a3.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a3.zlib.obs)
    stopifnot(identical(a3.zlib.obs, a3.exp))
}

##
## string: case-4
##
a4.in <- list(c("a", "bb"))
a4.exp <- list(list("a", "bb"))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a4.in), filename=filename, compression=FALSE, version='MAT5')
a4.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a4.obs)
stopifnot(identical(a4.obs, a4.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a4.in), filename=filename, compression=TRUE, version='MAT5')
    a4.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a4.zlib.obs)
    stopifnot(identical(a4.zlib.obs, a4.exp))
}
